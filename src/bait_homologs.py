#! /usr/bin/python3

import os
import sys
import argparse
import newick3
import logging
import datetime
from utils import parse_fasta
from search_proteomes import search_proteomes
from fasta_to_tree import fasta_to_tree
from trim_tree import trim_tree
from mask_monophyly import mask_monophyly, get_names_to_exclude
from cut_internal_branches import cut_internal_branches
from write_fasta_from_tree import write_fasta_from_tree
from tree_utils import get_front_labels


def check_bait_presence(seqlist, tree):
    with open(tree, "r") as inf:
        tre = newick3.parse(inf.readline())
    taxa = get_front_labels(tre)
    if set(seqlist).intersection(set(taxa)):
        return True
    else:
        return False


def fasta_to_subtree(fa, alner, treblder, relcut, abscut, intcut, mintaxa,
                     nt, ignore=[], mask=True, para=True, accurate=False):
    cln, t = fasta_to_tree(fa, nt, alner, treblder, accurate)
    tt = trim_tree(t, relcut, abscut)
    if mask:
        ttmm = mask_monophyly(tt, cln, para, ignore)
        subtrees = cut_internal_branches(ttmm, intcut, mintaxa)
    else:
        subtrees = cut_internal_branches(tt, intcut, mintaxa)
    return subtrees


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--blast", help="Use blastp for similarity \
                        search instead of default hmmsearch",
                        action="store_true")
    parser.add_argument("--min_bitscore", help="Filter blastp hits with \
                        bitscore lower than min_bitscore (default 30.0)",
                        type=float, default=30.0)
    parser.add_argument("--threshold", help="Filter blastp hits with bitscore \
                        lower than threshold * max bitscore of query \
                        (default 0.1)", type=float, default=0.1)
    parser.add_argument("-a", "--aligner", help="Alignment software to use: \
                        mafft (auto, default), fsa (defaults to --fast)",
                        default="mafft")
    parser.add_argument("-t", "--tree_builder", help="Tree building software \
                        to use: fasttree (wag, default), raxml-ng (defaults to\
                        WAG+G)", default="fasttree")
    parser.add_argument("-tc", "--tip_abs_cutoff", help="Absolute branch \
                        length cutoff for trimming tips. Tips longer than \
                        this will be trimmed. Defaults to 1.5", default=1.5)
    parser.add_argument("-tcf", "--tip_abs_cutoff_final", help="Absolute \
                        branch length cutoff for trimming tips in the final \
                        round. Tips longer than this will be trimmed. \
                        Defaults to 1.0", type=float, default=1.0)
    parser.add_argument("-rc", "--tip_rel_cutoff", help="Relative branch \
                        length cutoff for trimming tips. Tips longer than \
                        this and at least 10x longer than sister will be \
                        trimmed. Defaults to 1.0", default=1.0)
    parser.add_argument("-rcf", "--tip_rel_cutoff_final", help="Relative \
                        branch length cutoff for trimming tips in the final \
                        round. Tips longer than this and at least 10x longer \
                        than sister will be trimmed. Defaults to 0.5",
                        type=float, default=0.5)
    parser.add_argument("-ic", "--internal_cutoff", help="Branch length \
                        cutoff for internal branches. Subtrees subtended by \
                        branches longer than this will be trimmed. Defaults \
                        to 1.0", default=1.0)
    parser.add_argument("-icf", "--internal_cutoff_final", help="Branch \
                        length cutoff for internal branches in the final \
                        round. Subtrees subtended by branches longer than \
                        this will be trimmed. Defaults to 0.8", type=float,
                        default=0.8)
    parser.add_argument("-mt", "--min_taxa", help="Minimum taxa in a subtree \
                        to conserve and check for bait presence. Defaults to \
                        4", default=4)
    parser.add_argument("-nt", "--threads", help="Number of threads to use. \
                        Defaults to 2", default=2)
    parser.add_argument("-m", "--mask", help="If this flag is selected, \
                        monophyletic masking will be conducted",
                        action="store_true")
    parser.add_argument("-mp", "--mask_paraphyly", help="Whether to mask \
                        paraphyletic sequences while doing monophyletic \
                        masking", action="store_true")
    parser.add_argument("-if", "--ignore_file", help="File containing taxon \
                        names (one per line) to ignore while masking \
                        monophyletic tips. Defaults masks all taxa",
                        default=None)
    parser.add_argument("-po", "--prune_og", help="Whether to extract rooted \
                        ingroup clades containing baits after first round of \
                        tree inference (requires OG file, default off) - TBI",
                        action="store_true")
    parser.add_argument("-og", "--outgroups", help="File containing outgroup \
                        taxon labels, one per line - TBI")
    parser.add_argument("-it", "--iterate", help="how many times to iterate \
                        tree building and cleaning", type=int, default=3)
    parser.add_argument("-o", "--output_dir", help="Directory to put output. \
                        Defaults to current directory", default=os.getcwd())
    parser.add_argument("-k", "--keep", help="Number of hits to keep (default \
                        all)", type=int, default=None)
    parser.add_argument("-dbl", "--dblist", help="Text file with files in \
                        database_dir to to search, if not all, one per line")
    parser.add_argument("bait", help="FASTA file of baits to search")
    parser.add_argument("database_dir", help="Path to the database containing \
                        proteomes to search. Expects file endings of .pep.fa \
                        or .cdhit")
    args = parser.parse_args()

    logging.basicConfig(filename='log.txt', encoding='utf-8',
                        level=logging.DEBUG, filemode="w",
                        format="%(levelname)s|%(message)s")
    logging.info(datetime.datetime.now().strftime("%d.%b %Y %H:%M:%S"))
    logging.info(" ".join([sys.executable] + sys.argv))

    if "/" in args.bait:
        name = args.bait.split("/")[-1].split(".")[0]
    else:
        name = args.bait.split(".")[0]

    BAITS = [key for key in dict(parse_fasta(args.bait))]
    logging.info(f"querying with {len(BAITS)} sequences in {args.bait}")

    if args.ignore_file is not None:
        IGNORE = get_names_to_exclude(args.ignore_file)
        with open("ignored.txt", "w") as outf:
            for i in IGNORE:
                outf.write(f"{i}\n")
        logging.info(f"during masking, ignoring {len(IGNORE)} taxa in " 
                     f"{args.ignore_file}")
        logging.info(f"writing to ignored.txt")
    else:
        IGNORE = []

    if args.prune_og:
        if not args.outgroups:
            sys.stderr.write("Must specify outgroups with -og/--outgroups")
            sys.exit()
        OUTGROUPS = []
        with open(args.outgroups, "r", encoding="utf-8") as f:
            for line in f:
                OUTGROUPS.append(line.strip())

    if args.dblist is not None:
        DBLIST = []
        with open(args.dblist, "r", encoding="utf-8") as f:
            for line in f:
                DBLIST.append(line.strip())
        hits = search_proteomes(args.bait, args.database_dir, args.output_dir,
                                args.blast, args.keep, args.threads,
                                args.min_bitscore, args.threshold, DBLIST)
    else:
        hits = search_proteomes(args.bait, args.database_dir, args.output_dir,
                                args.blast, args.keep, args.threads,
                                args.min_bitscore, args.threshold)

    iters = args.iterate
    # first round
    subtrees = fasta_to_subtree(hits, args.aligner, args.tree_builder,
                                args.tip_rel_cutoff, args.tip_abs_cutoff,
                                args.internal_cutoff, args.min_taxa,
                                args.threads, IGNORE, args.mask,
                                args.mask_paraphyly)
    fas = []
    for t in subtrees:
        if check_bait_presence(BAITS, t):
            newf = write_fasta_from_tree(hits, t)
            logging.info(f"found baits in {t}, writing fasta to {newf}")
            fas.append(newf)
        else:
            os.remove(t)
            # comment out this if you want to keep
    iters -= 1  # now 2 in default
    while iters > 0:
        subtrees = []
        if iters == 1:  # last round, can change
            # to do mafft --genafpair and iqtree
            for f in fas:
                subtrees += fasta_to_subtree(f, args.aligner,
                                             args.tree_builder,
                                             args.tip_rel_cutoff_final,
                                             args.tip_abs_cutoff_final,
                                             args.internal_cutoff_final,
                                             args.min_taxa,
                                             args.threads, IGNORE,
                                             args.mask,
                                             args.mask_paraphyly,
                                             accurate=False)
            for t in subtrees:
                if check_bait_presence(BAITS, t):
                    newf = write_fasta_from_tree(hits, t)
                    logging.info(f"found baits in {t}, writing fasta to {newf}")
                    fas.append(newf)
                else:
                    os.remove(t)
                    # comment out this if you want to keep
        else:
            for f in fas:
                subtrees += fasta_to_subtree(f, args.aligner,
                                             args.tree_builder,
                                             args.tip_rel_cutoff,
                                             args.tip_abs_cutoff,
                                             args.internal_cutoff,
                                             args.min_taxa,
                                             args.threads, IGNORE,
                                             args.mask,
                                             args.mask_paraphyly)
            fas = []
            for t in subtrees:
                if check_bait_presence(BAITS, t):
                    newf = write_fasta_from_tree(hits, t)
                    logging.info(f"found baits in {t}, writing fasta to {newf}")
                    fas.append(newf)
                else:
                    os.remove(t)
                    # comment out this if you want to keep
        iters -= 1
    print("Done!")
