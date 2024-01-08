#! /usr/bin/python

import sys
import os
import subprocess
import argparse
import shutil
import logging
from utils import parse_fasta


def fasta_to_stockholm(inf):
    nseq = 0
    with open(inf, "r") as seqf:
        for line in seqf:
            nseq += 1
    if nseq == 2:  # 1 sequence
        print("Single sequence, going straight to HMM")
        shutil.copy(inf, inf + ".sto")
        # name .sto for compatibility, but actually .fa
    else:
        print("Aligning baits with FSA")
        cmd = ["fsa", "--fast", "--stockholm", inf]
        print(subprocess.list2cmdline(cmd))
        out = subprocess.run(cmd, shell=False, capture_output=True,
                             text=True).stdout
        with open(inf + ".sto", "w") as outf:
            outf.write(out)


def stockholm_to_profile(sto):
    print("Building hmm profile")
    cmd = ["hmmbuild", sto[:-3] + "hmm", sto]
    print(subprocess.list2cmdline(cmd))
    subprocess.run(cmd, shell=False)


def hmmsearch_db(hmm, dbf):
    print("Searching database with hmm")
    cmd = ["hmmsearch", "--noali", "--tblout", dbf + ".out", hmm, dbf]
    print(subprocess.list2cmdline(cmd))
    subprocess.run(cmd, shell=False)


def make_blast_db(dbf):
    print("Making blastdb")
    cmd = ["makeblastdb", "-in", dbf, "-dbtype", "prot", "-parse_seqids"]
    print(subprocess.list2cmdline(cmd))
    subprocess.run(cmd, shell=False)


def blast_db(bait, dbf, nt):
    print("Searching database with blastp")
    blastout = dbf + "." + bait + ".blastp.outfmt6"
    cmd = ["blastp", "-query", bait, "-db", dbf, "-num_threads", str(nt),
           "-evalue", "10", "-out", blastout, "-outfmt",
           "6 qseqid sseqid evalue bitscore"]
    print(subprocess.list2cmdline(cmd))
    subprocess.run(cmd, shell=False)
    return blastout


def parse_hmmsearch_out(outf):
    """"Parse hmmsearch output with -tblout"""
    print("Parsing hmmsearch output")
    hits = open(outf[:-3] + "hits", "w")
    with open(outf, "r") as search_out:
        for line in search_out:
            if not line.startswith("#"):
                hits.write(line.split(" ")[0].strip() + "\n")
    hits.close()


def parse_blastp_out(outf, min_bitscore=30.0, thresh=0.1):
    """Parse blastp outfmt 6 with columns qseqid sseqid evalue bitscore
    By default ignores hits with bitscore < 30.0 and < 0.1 * max bitscore of
    that query"""
    sys.stderr.write("Parsing blastp output\n")
    hitsfilename = outf[:-7] + "hits"
    hitsfile = open(hitsfilename, "w")
    hits = []
    with open(outf, "r") as f:
        for line in f:
            # out is four columns, qid, sid, evalue, bitscore
            line = line.strip().split("\t")
            q, h, b = line[0], line[1], float(line[3])
            if q != h and b > min_bitscore:  # filter self and lower than min
                hits.append((h, b))
    highest = 0.0
    currentq = ""
    hits_filtered = []  # contains hits after bitscore filtering
    for h in hits:
        # here we filter hits with respect to highest score of that query
        if h[0] == currentq:
            continue
        else:
            currentq = h[0]
            highest = h[1]
        if h[1] > thresh * highest:
            hits_filtered.append(h)
    hits_filtered = [x[0] for x in sorted(hits_filtered, key=lambda x: x[1],
                                          reverse=True)]
    hits_unique = list(set(hits_filtered))
    for h in hits_unique:
        hitsfile.write(h + "\n")
    hitsfile.close()
    return hitsfilename


def gather_sequences(hitsf, dbf, outf, nhits: int = None):
    """Takes a per db hits file and agglomerates sequences for hits,
    optionally with a numerical cutoff of hits to keep"""
    print("Compiling sequences")
    seqout = open(outf, "a")
    if nhits is not None:
        seql = [x.rstrip("\n") for x in open(hitsf, "r").readlines()[:nhits]]
    else:
        seql = [x.rstrip("\n") for x in open(hitsf, "r").readlines()]
    dbdict = dict([x for x in parse_fasta(dbf)])
    for s in seql:
        seqout.write(">" + s + "\n")
        seqout.write(dbdict[s] + "\n")
    seqout.close()


def search_proteomes(bait, database_dir, output_dir, blast=False, nhits=None,
                     nt=1, min_bitscore=30.0, thresh=0.1, dbkeep=None):
    # file name
    if "/" in bait:
        name = bait.split("/")[-1].split(".")[0]
    else:
        name = bait.split(".")[0]

    # collect individual sequence sets
    dblist = []
    outf = open("database.txt", "w")
    for dirpath, _, filenames in os.walk(database_dir):
        for f in filenames:
            if f.endswith(".pep.fa") or f.endswith(".cdhit"):  # add suffix
                if dbkeep is not None:
                    if f in dbkeep:
                        dblist.append(os.path.abspath(os.path.join(dirpath,
                                                                    f)))
                        outf.write(f"{f}\n")
                else:
                    dblist.append(os.path.abspath(os.path.join(dirpath, f)))
                    outf.write(f"{f}\n")
    outf.close()

    if dbkeep is not None:
        logging.info(f"searching only {len(dblist)} sequence DBs "
                                 f"in DB list")
    else:
        logging.info(f"searching {len(dblist)} sequence DBs "
                                 f"in {database_dir}")
        
    if blast:
        logging.info("using blastp")
        outfile = os.path.abspath(output_dir) + "/" + name + ".blastp.fa"
        if os.path.isfile(outfile):
            os.remove(outfile)  # prevent appending partial file
        for db in dblist:
            blastdbsuf = [".pdb", ".phr", ".pin", ".pog", ".pos", ".pot",
                          ".psq", ".ptf", ".pto"]
            for s in blastdbsuf:
                if not os.path.isfile(db + s):
                    make_blast_db(db)
            blastout = blast_db(bait, db, nt)
            hits = parse_blastp_out(blastout, min_bitscore, thresh)
            gather_sequences(hits, db, outfile, nhits)
            os.remove(blastout)
            os.remove(hits)
    else:
        logging.info("using hmmsearch")
        outfile = os.path.abspath(output_dir) + "/" + name + ".hmmsearch.fa"
        if os.path.isfile(outfile):
            os.remove(outfile)  # prevent appending partial file
        stockholm = bait + ".sto"
        if stockholm not in os.listdir(os.getcwd()):
            fasta_to_stockholm(bait)
            stockholm_to_profile(stockholm)
        else:
            try:
                stockholm_to_profile(stockholm)
            except subprocess.CalledProcessError:
                fasta_to_stockholm(bait)
                stockholm_to_profile(stockholm)
        for db in dblist:
            hmmsearch_db(bait + ".hmm", db)
            parse_hmmsearch_out(db + ".out")
            gather_sequences(db + ".hits", db,
                             outfile,
                             nhits=nhits)
            os.remove(db + ".out")
            os.remove(db + ".hits")

    baitdict = dict([x for x in parse_fasta(bait)])
    with open(outfile, "a") as outf:  # append baits back to search results
        for key, value in baitdict.items():
            outf.write(">" + key + "\n")
            outf.write(value + "\n")
    return outfile


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("bait", help="FASTA file of query sequences")
    parser.add_argument("database_dir", help="Directory containing sequence \
                        databases to be searched")
    parser.add_argument("output_dir", help="Directory to put output")
    parser.add_argument("-k", "--keep", help="Number of hits to retain \
                        (default all)", type=int)
    parser.add_argument("-b", "--blast", help="Use blastp for similarity \
                        searches, instead of hmmer", action="store_true")
    parser.add_argument("--min_bitscore", help="For blast. Filter hits with \
                        bitscores lower than min_bitscore (default 30.0)",
                        type=float, default=30.0)
    parser.add_argument("--threshold", help="For blast. Filter hits with \
                        bitscores lower than threshold * max bitscore of \
                        query (default 0.1)", type=float, default=0.1)
    parser.add_argument("-t", "--threads", help="Number of threads", type=int,
                        default=1)
    parser.add_argument("-dbl", "--dblist", help="Text file with files in \
                        database_dir to to search, if not all, one per line")
    # parser.add_argument("")
    args = parser.parse_args()

    # _ = search_proteomes(args.bait, args.database_dir, args.output_dir,
    #                      blast=False, nhits=args.keep)
    if args.dblist is not None:
        DBLIST = []
        with open(args.dblist, "r") as f:
            for line in f:
                DBLIST.append(line.strip())
        _ = search_proteomes(args.bait, args.database_dir, args.output_dir,
                             args.blast, args.keep, args.threads,
                             args.min_bitscore, args.threshold, DBLIST)
    else:
        _ = search_proteomes(args.bait, args.database_dir, args.output_dir,
                             args.blast, args.keep, args.threads,
                             args.min_bitscore, args.threshold)
