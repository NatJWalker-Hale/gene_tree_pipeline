#! /usr/bin/python3

import sys
import argparse
import newick3
import logging
from tree_utils import get_name, remove_kink
from utils import parse_fasta


def get_names_to_exclude(ignoref):
    ignore = [line.strip() for line in open(ignoref, "r").readlines()]
    return ignore


def mask_monophyletic_tips(curroot, unamb_chrDICT, ignore=[]):
    going = True
    while going and curroot is not None and len(curroot.leaves()) >= 4:
        going = False
        for node in curroot.iternodes():   # walk through nodes
            if not node.istip:
                continue   # only look at tips
            name = get_name(node.label)
            if name in ignore:
                continue   # do not mask the genomes
            for sister in node.get_sisters():
                if sister.istip and name == get_name(sister.label):  # mask
                    if unamb_chrDICT[node.label] > unamb_chrDICT[sister.label]:
                        node = sister.prune()
                    else:
                        node = node.prune()
                    if len(curroot.leaves()) >= 4:
                        if (
                            node == curroot
                            and node.nchildren == 2
                            ) or (
                                node != curroot
                                and node.nchildren == 1):
                            node, curroot = remove_kink(node, curroot)
                    going = True
                    break
    return curroot


def mask_paraphyletic_tips(curroot, unamb_chrDICT, ignore=[]):
    going = True
    while going and curroot is not None and len(curroot.leaves()) >= 4:
        going = False
        for node in curroot.iternodes():  # walk through nodes
            if not node.istip:
                continue  # only look at tips
            name = get_name(node.label)
            if name in ignore:
                continue  # do not mask the genomes
            parent = node.parent
            if node == curroot or parent == curroot or parent is None:
                continue  # no paraphyletic tips for the root
            for para in parent.get_sisters():
                if para.istip and name == get_name(para.label):  # mask
                    if unamb_chrDICT[node.label] > unamb_chrDICT[para.label]:
                        node = para.prune()
                    else:
                        node = node.prune()
                    if len(curroot.leaves()) >= 4:
                        if (
                            node == curroot
                            and node.nchildren == 2
                            ) or (
                                node != curroot
                                and node.nchildren == 1):
                            node, curroot = remove_kink(node, curroot)
                    going = True
                    break
    return curroot


def mask(curroot, clnfile, para=True, ignore=[]):
    chrDICT = {}  # key is seqid, value is number of unambiguous chrs
    for key, value in dict([x for x in parse_fasta(clnfile)]).items():
        for ch in ['-', 'X', "x", "?", "*"]:
            value = value.replace(ch, "")  # ignore gaps, xs and Xs
        chrDICT[key] = len(value)
    curroot = mask_monophyletic_tips(curroot, chrDICT, ignore)
    if para:
        curroot = mask_paraphyletic_tips(curroot, chrDICT, ignore)
    return curroot


def mask_monophyly(tre, clnaln, para=True, ignore=[]):
    with open(tre, "r") as inf:
        intree = newick3.parse(inf.readline())
    in_tips = len(intree.leaves())
    curroot = mask(intree, clnaln, para, ignore)
    out_tips = len(curroot.leaves())
    masked = tre + ".mm"
    if para:
        logging.info(f"masking monophyletic and paraphyletic tips in {tre}")
    else:
        logging.info(f"masking monophyletic tips in {tre}")
    logging.info(f"masked {in_tips - out_tips} tips, writing to {masked}")
    with open(masked, "w") as outf:
        outf.write(newick3.tostring(curroot)+";\n")
    return masked


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--mask_paraphyly",
                        help="Mask tips that are paraphyletic (default True)",
                        default=True)
    parser.add_argument("-f", "--exclude_file",
                        help="File containing taxon names (one per line) \
                        to ignore while masking",
                        default=None)
    parser.add_argument("tree",
                        help="The tree file in newick format to mask tips")
    parser.add_argument("aln_cln",
                        help="The cleaned alignment in FASTA format \
                        corresponding to the tips of the tree")
    args = parser.parse_args()
    if args.mask_paraphyly:
        para = True
    else:
        para = False
    if args.exclude_file is not None:
        ignore = get_names_to_exclude(args.exclude_file)
    else:
        ignore = []
    _ = mask_monophyly(args.tree, args.aln_cln, para, ignore)
