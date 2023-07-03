#! /usr/bin/python3


import sys
import argparse
import newick3


def get_taxon_code_name_corres(inf):
    taxonDict = {}
    with open(inf, "r") as taxTable:
        for line in taxTable:
            code, name = line.strip().split("\t")
        taxonDict[code] = name
    return taxonDict


# def relabel_tree(tre):
#     for n in tre.leaves():
