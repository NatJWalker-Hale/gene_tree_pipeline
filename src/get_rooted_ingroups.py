#! /usr/bin/python3


import os
import sys
import argparse
import newick3
import tree_utils


if __name__ == "__main__":
    ingroups = []
    with open(sys.argv[1], "r") as igf:
        for line in igf:
            ingroups.append(line.strip())
    print(ingroups)

    outgroups = []
    with open(sys.argv[2], "r") as ogf:
        for line in ogf:
            outgroups.append(line.strip())
    print(outgroups)

    t = newick3.parse_from_file(sys.argv[3])
    # for n in t.iternodes():
    #     print(n.label)

    rt1 = tree_utils.extract_rooted_ingroup_clades(t, ingroups, outgroups, 1)
    print(rt1)
    print([newick3.to_string(x) for x in rt1])

    t = newick3.parse_from_file(sys.argv[3])
    rt2 = tree_utils.extract_rooted_ingroup_clades_with_outgroup(t, ingroups,
                                                                 outgroups, 1)
    print(rt2)
    print([newick3.to_string(x) for x in rt2])
