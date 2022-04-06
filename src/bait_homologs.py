#! /usr/bin/python3

import os
import sys
import argparse
import shutil
import newick3
import phylo3
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


def check_same_tree(tree1, tree2):
    with open(tree1, "r") as inf1:
        tre1 = newick3.parse(inf1.readline())
    taxa1 = get_front_labels(tre1)
    with open(tree2, "r") as inf2:
        tre2 = newick3.parse(inf2.readline())
    taxa2 = get_front_labels(tre2)
    if taxa1 == taxa2:
        return True
    else:
        return False


def sub_loop(fa, alner, treeblder, abscut, relcut, nt, ignore=[],
             mask=True, para=True):
    fasta_to_tree(fa, nt, alner, treeblder)
    trim_tree(fa + "." + alner + ".aln-cln." + treeblder + ".tre", relcut,
              abscut)
    if mask:
        mask_monophyly(fa + "." + alner + ".aln-cln." + treeblder + ".tre.tt",
                       fa + "." + alner + ".aln-cln", para, ignore)


def main_loop(bait, fa, alner, treeblder, abscut, relcut, intcut, mintaxa,
              nt, ignore=[], mask=True, para=True, iterate=False):
    bait_seqs = [key for key in dict([x for x in parse_fasta(bait)]).keys()]
    fasta_to_tree(fa, nt, alner, treeblder)
    trim_tree(fa + "." + alner + ".aln-cln." + treeblder + ".tre", relcut,
              abscut)
    if mask:
        mask_monophyly(fa + "." + alner + ".aln-cln." + treeblder + ".tre.tt",
                       fa + "." + alner + ".aln-cln", para, ignore)
        subtrees = cut_internal_branches(fa + "." + alner + ".aln-cln." +
                                         treeblder + ".tre.tt.mm", intcut,
                                         mintaxa)
    else:
        subtrees = cut_internal_branches(fa + "." + alner + ".aln-cln." +
                                         treeblder + ".tre.tt", intcut,
                                         mintaxa)
    print(subtrees)
    if subtrees is not None and len(subtrees) > 1:
        for t in subtrees:
            if check_bait_presence(bait_seqs, t):
                print(bait + " in " + t)
                newfa = write_fasta_from_tree(fa, t)
                main_loop(bait, newfa, alner, treeblder, abscut, relcut,
                          intcut, mintaxa, nt, ignore, mask, para, iterate)
            else:
                os.remove(t)
    else:
        if iterate:
            print("No more subtrees to cut. Iterating tip trimming and \
                   monophyletic masking until topology stabilises.")
            newfa = write_fasta_from_tree(fa, subtrees[0])
            tree1 = subtrees[0]
            name = subtrees[0].split(".")[0]+"_m_1"
            print(ignore)
            sub_loop(newfa, alner, treeblder, abscut, relcut, nt, ignore,
                     mask, para)
            if mask:
                shutil.copyfile(newfa + "." + alner + ".aln-cln." + treeblder +
                                ".tre.tt.mm", name+".tre")
            else:
                shutil.copyfile(newfa + "." + alner + ".aln-cln." + treeblder +
                                ".tre.tt", name + ".tre")
            tree2 = name + ".tre"
            going = True
            if check_same_tree(tree1, tree2):
                going = False
            counter = 2
            while going:
                tree1 = tree2
                name = tree1.split(".")[0].rsplit("_", 1)[0]+"_"+str(counter)
                newfa = write_fasta_from_tree(fa, tree1)
                print(ignore)
                sub_loop(newfa, alner, treeblder, abscut, relcut, nt, ignore,
                         mask, para)
                if mask:
                    shutil.copyfile(newfa + "." + alner + ".aln-cln." +
                                    treeblder + ".tre.tt.mm", name+".tre")
                else:
                    shutil.copyfile(newfa + "." + alner + ".aln-cln." +
                                    treeblder + ".tre.tt", name+".tre")
                tree2 = name+".tre"
                counter += 1
                if check_same_tree(tree1, tree2):
                    going = False
            print("Finished iteration")
            _ = write_fasta_from_tree(fa, tree2)
        else:
            _ = write_fasta_from_tree(fa, subtrees[0])
            print("No iteration requested. Finished")


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--aligner", help="Alignment software to use: \
                        mafft (auto, default), fsa (defaults to --fast)",
                        default="mafft")
    parser.add_argument("-t", "--tree_builder", help="Tree building software \
                        to use: fasttree (wag, default), iqtree (defaults to \
                        WAG+G)", default="fasttree")
    parser.add_argument("-tc", "--tip_abs_cutoff", help="Absolute branch \
                        length cutoff for trimming tips. Tips longer than \
                        this will be trimmed. Defaults to 1.5", default=1.5)
    parser.add_argument("-rc", "--tip_rel_cutoff", help="Relative branch \
                        length cutoff for trimming tips. Tips longer than \
                        this and at least 10x longer than sister will be \
                        trimmed. Defaults to 1.0", default=1.0)
    parser.add_argument("-ic", "--internal_cutoff", help="Branch length \
                        cutoff for internal branches. Subtrees subtended by \
                        branches longer than this will be trimmed. Defaults \
                        to 1.0", default=1.0)
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
    parser.add_argument("-it", "--iterate", help="If this flag is selected, \
                        tip trimming and monophyletic masking will be \
                        iterated until topology has stabilised following \
                        final subtree cut", action="store_true")
    parser.add_argument("-o", "--output_dir", help="Directory to put output. \
                        Defaults to current directory", default=os.getcwd())
    parser.add_argument("-k", "--keep", help="Number of hits to keep (default \
                        all)", type=int, default=None)
    parser.add_argument("bait", help="FASTA file of baits to search")
    parser.add_argument("database_dir", help="Path to the database containing \
                        proteomes to search. Expects file endings of .pep.fa \
                        or .cdhit")
    args = parser.parse_args()

    if "/" in args.bait:
        name = args.bait.split("/")[-1].split(".")[0]
    else:
        name = args.bait.split(".")[0]
    if args.ignore_file is not None:
        IGNORE = get_names_to_exclude(args.ignore_file)
    else:
        IGNORE = []
    search_proteomes(args.bait, args.database_dir, args.output_dir,
                     blast=False, nhits=args.keep)
    main_loop(args.bait, name+".hmmsearch.fa", args.aligner, args.tree_builder,
              args.tip_abs_cutoff, args.tip_rel_cutoff, args.internal_cutoff,
              args.min_taxa, args.threads, IGNORE, args.mask,
              args.mask_paraphyly, args.iterate)
