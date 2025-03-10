#! /usr/bin/python3

import sys
import os
import argparse
import newick3
import logging
from tree_utils import *


# from Y. Yang https://bitbucket.org/yangya/adh_2016/src/master/trim_tips.py
# return the outlier tip, with abnormal high contrast and long branch
def check_countrast_outlier(node0, node1, above0, above1, relative_cutoff):
    if node0.istip and above0 > relative_cutoff:
        if above1 == 0.0 or above0/above1 > 10:
            return node0
    if node1.istip and above1 > relative_cutoff:
        if above0 == 0.0 or above1/above0 > 10:
            return node1
    return None

# from Y. Yang https://bitbucket.org/yangya/adh_2016/src/master/trim_tips.py
def remove_a_tip(root, tip_node):
    node = tip_node.prune()
    if len(root.leaves()) > 3:
        node, root = remove_kink(node, root)
        return root
    else:
        print("Fewer than four tips left")
        return None

# from Y. Yang https://bitbucket.org/yangya/adh_2016/src/master/trim_tips.py
def trim(curroot, relative_cutoff, absolute_cutoff):
    if curroot.nchildren == 2:
        temp, root = remove_kink(curroot, curroot)
    going = True
    while going and curroot is not None and len(curroot.leaves()) > 3:
        going = False
        for i in curroot.iternodes(order=1):  # POSTORDER
            if i.nchildren == 0:  # at the tip
                i.data['len'] = i.length
                if i.length > absolute_cutoff:
                    curroot = remove_a_tip(curroot, i)
                    going = True
                    break
            elif i.nchildren == 1:  # kink in tree
                remove_kink(i, curroot)
                going = True
                break
            elif i.nchildren == 2:  # normal bifurcating internal nodes
                child0, child1 = i.children[0], i.children[1]
                above0, above1 = child0.data['len'], child1.data['len']
                i.data['len'] = ((above0+above1)/2.) + i.length
                # stepwise average
                outlier = check_countrast_outlier(child0, child1, above0,
                                                  above1, relative_cutoff)
                if outlier is not None:
                    curroot = remove_a_tip(curroot, outlier)
                    going = True  # need to keep checking
                    break
            else:  # 3 or more branches from this node. Pair-wise comparison
                total_len = 0
                nchild = i.nchildren
                for child in i.children:
                    total_len += child.data['len']
                i.data['len'] = total_len / float(i.nchildren)
                keep_checking = True
                for index1 in range(nchild):  # do all the pairwise comparison
                    for index2 in range(nchild):
                        if index2 <= index1:
                            continue  # avoid repeatedly checking a pair
                        child1, child2 = i.children[index1], i.children[index2]
                        above1, above2 = child1.data['len'], child2.data['len']
                        outlier = check_countrast_outlier(child1, child2,
                                                          above1, above2,
                                                          relative_cutoff)
                        if outlier is not None:
                            print(above1, above2)
                            curroot = remove_a_tip(curroot, outlier)
                            going = True  # need to keep checking
                            keep_checking = False  # to break the nested loop
                            break
                    if not keep_checking:
                        break
    return curroot


def trim_tree(inf, relative_cut=1.0, absolute_cut=1.5):
    with open(inf, "r") as infile:
        intree = newick3.parse(infile.readline())
    in_tips = len(intree.leaves())
    outtree = trim(intree, float(relative_cut), float(absolute_cut))
    out_tips = len(outtree.leaves())
    # if outtree is not None:
    out = inf + ".tt"
    logging.info(f"trimming tips in {inf} longer than {absolute_cut} "
                f"or >10x length of sister and longer than {relative_cut}")
    logging.info(f"trimmed {in_tips - out_tips} tip(s), writing to {out}")
    with open(out, "w") as outfile:
        outfile.write(newick3.tostring(outtree)+";\n")
    return out


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--absolute_cutoff", help="Trim terminal edges \
                        longer than absolute cutoff (default 1.5)",
                        default=1.5)
    parser.add_argument("-r", "--relative_cutoff", help="Trim terminal edges \
                        longer than relative cutoff and > 10x length of \
                        sister (default 1.0)",
                        default=1.0)
    parser.add_argument("intree", help="Tree in newick format to trim long \
                        branches")
    args = parser.parse_args()
    _ = trim_tree(args.intree, args.absolute_cutoff, args.relative_cutoff)
