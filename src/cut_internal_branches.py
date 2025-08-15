#! /usr/bin/python3

import sys
import os
import argparse
import newick3
import logging
from tree_utils import get_front_names, remove_kink
from shutil import copy


def count_taxa(node):
    """given a node count how many taxa it has in frount"""
    return len(set(get_front_names(node)))


def cut_long_internal_branches(curroot, cutoff, mintaxa=4):
    """cut long branches and output all subtrees >= mintaxa"""
    going = True
    subtrees = []  # store all subtrees after cutting
    while going:
        going = False
        # only keep going if long branches were found during last round
        for node in curroot.iternodes():  # Walk through nodes
            if node == curroot:  # skip root
                continue
            if node.istip:  # skip tips
                continue
            if node.nchildren == 1:
                node,curroot = remove_kink(node, curroot)
                going = True
                break
            child0, child1 = node.children[0], node.children[1]
            if node.length > cutoff:
                print(node.length)
                if not child0.istip and not child1.istip and child0.length + child1.length > cutoff:
                    # resulting subtree would have two clades separated
                    # by length longer than cutoff, hence prune already
                    print(child0.length + child1.length)
                    if count_taxa(child0) >= int(mintaxa):
                        subtrees.append(child0)
                    if count_taxa(child1) >= int(mintaxa):
                        subtrees.append(child1)
                else:
                    subtrees.append(node)
                node = node.prune()
                if len(curroot.leaves()) > 2:  # no kink if only two left
                    node, curroot = remove_kink(node, curroot)
                    going = True
                break
    if count_taxa(curroot) >= int(mintaxa):
        subtrees.append(curroot)  # write out the residue after cutting
    return subtrees


def cut_internal_branches(tre, brlencutoff=1.0, mintaxa=4):
    mintaxa = int(mintaxa)
    brlencutoff = float(brlencutoff)
    print("Cutting at branches longer than "+str(brlencutoff))
    logging.info(f"cutting at branches longer than {brlencutoff} in {tre}")
    subtree_names = []
    with open(tre, "r") as inf:
        intree = newick3.parse(inf.readline())
        subtrees = sorted([x for x in
                           cut_long_internal_branches(intree, brlencutoff,
                                                      mintaxa)],
                          reverse=True, key=lambda x: count_taxa(x))
        print(subtrees)
        if len(subtrees) == 0:
            print("No branches to cut in "+tre)
            logging.info(f"no branches to cut in {tre}, writing input to "
                         f"{tre.split('.')[0]}_1.subtree")
            subtree_names.append(tre.split(".")[0]+"_1.subtree")
        else:
            count = 0
            for t in subtrees:
                if count_taxa(t) >= mintaxa:
                    if t.nchildren == 2:  # fix bifurcating roots from cutting
                        _, t = remove_kink(t, t)
                    count += 1
                    sub_name = "".join([tre.split(".")[0],
                                        "_",
                                        str(count),
                                        ".subtree"])
                    print(sub_name)
                    subtree_names.append(sub_name)
                    with open(sub_name, "w") as outfile:
                        outfile.write(newick3.tostring(t)+";\n")
    logging.info(f"writing {len(subtree_names)} subtree(s)")
    return subtree_names


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-bc", "--brlencut",
                        help="Cut internal branches longer than this value \
                        (default 1.0)", default=1.0, type=float)
    parser.add_argument("-m", "--mintaxa", help="Minimum taxa in subtree to \
                        retain (default 4)", default=4, type=int)
    parser.add_argument("intree", help="Tree file in newick format to cut \
                        internal branches")
    args = parser.parse_args()

    tref = os.path.abspath(args.intree)
    # print(tref)
    cut_internal_branches(tref, args.brlencut, args.mintaxa)
