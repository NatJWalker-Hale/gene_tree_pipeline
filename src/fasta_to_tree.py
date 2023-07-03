#! /usr/bin/python3

import sys
import subprocess
import argparse


def fasta_to_aln(inf, thread=2, alner="mafft", accurate=False):
    if alner == "mafft":
        if accurate:
            cmd = ["mafft", "--genafpair", "--maxiterate", "1000", "--amino",
                   "--thread", str(thread), inf]
            out = inf + ".mafftgp.aln"
            outf = open(out, "w")
        else:
            cmd = ["mafft", "--auto", "--amino", "--thread", str(thread), inf]
            out = inf + ".mafft.aln"
            outf = open(out, "w")
    elif alner == "fsa":
        cmd = ["fsa", "--fast", inf]
        out = inf + ".fsa.aln"
        outf = open(out, "w")
    print("Aligning sequences with "+alner)
    print(subprocess.list2cmdline(cmd))
    proc = subprocess.Popen(cmd, shell=False, stdout=outf,
                            stderr=subprocess.PIPE)
    while proc.poll() is None:
        print(proc.stderr.readline().decode().strip())
    return out


def clean(aln):
    cleaned = aln + "-cln"
    cmd = ["pxclsq", "-s", aln, "-o", cleaned, "-p", "0.1"]
    print("Cleaning alignment")
    print(subprocess.list2cmdline(cmd))
    subprocess.run(cmd, shell=False)
    return cleaned


def aln_to_tree(aln, treeblder="fasttree", thread=2):
    if treeblder == "fasttree":
        out = aln + ".fasttree.tre"
        cmd = ["fasttree", "-wag", "-out", out, aln]
    elif treeblder == "raxml-ng":
        # cmd = ["iqtree", "-s", aln, "-m", "WAG+G", "-nt", str(thread)]
        cmd = ["raxml-ng", "--search", "--msa", aln, "--model", "WAG+G",
               "--threads", str(thread)]
    print("Inferring tree with "+treeblder)
    print(subprocess.list2cmdline(cmd))
    proc = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE)
    proc.wait()
    if treeblder == "raxml-ng":
        out = aln + ".raxml.bestTree"
    return out


def fasta_to_tree(inf, thread=2, alner="mafft", treeblder="fasttree",
                  accurate=False):
    if accurate:
        aln = fasta_to_aln(inf, thread, alner, accurate)
    else:
        aln = fasta_to_aln(inf, thread, alner)
    cleaned = clean(aln)
    out = aln_to_tree(cleaned, treeblder, thread)
    return cleaned, out


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--aligner", help="The software used to infer \
                        alignments. Options are mafft (auto, default), fsa",
                        default="mafft")
    parser.add_argument("-t", "--tree_builder", help="The software used to \
                        infer trees. Options are fasttree (default), iqtree",
                        default="fasttree")
    parser.add_argument("-nt", "--threads", help="The number of threads to \
                        use for alignment and tree inference (default 2)",
                        default=2)
    parser.add_argument("fasta", help="FASTA of sequences to align and infer \
                         a tree")
    args = parser.parse_args()

    _, _ = fasta_to_tree(args.fasta, args.threads, args.aligner,
                         args.tree_builder)
