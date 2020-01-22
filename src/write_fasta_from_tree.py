#! /usr/bin/python3

import sys 
import os
import argparse
import newick3
import phylo3
from utils import parse_fasta
import tree_utils


def write_fasta_from_tree(allfa,tree):
    tree = os.path.abspath(tree)
    name = tree.split("/")[-1].split(".")[0]
    fa_dict = dict([x for x in parse_fasta(allfa)])
    with open(tree,"r") as inf:
        t = newick3.parse(inf.readline())
    seqids = set(tree_utils.get_front_labels(t))
    with open(name+".pep.fa","w") as outf:
        for s in seqids:
            outf.write(">"+s+"\n")
            outf.write(fa_dict[s]+"\n")
    return name+".pep.fa"

if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("intree",help="Tree file in newick format to write fasta from tips")
    parser.add_argument("allfasta",help="Master FASTA containing sequences to be written")
    args = parser.parse_args()

    write_fasta_from_tree(args.allfasta,args.intree)


    