#! /usr/bin/python3

import sys
import argparse
import subprocess
from utils import parse_fasta

def get_names_for_tree(inf):
    seqs = dict([x for x in parse_fasta(inf)])
    orig_names = [key for key in seqs.keys()]
    iq_names = [s.replace("@","_") for s in orig_names]
    return orig_names,iq_names

def run_iqtree(aln,model,nt,bs=1000,alrt=1000):
    """Runs IQ-TREE with the given model string (can also be replaced with MFP for model finder). By default does 1000 UFB reps and SH-aLRT"""
    cmd = ["iqtree","-s",aln,"-m",model,"-nt",nt,"-bb",bs,"-alrt",alrt]
    subprocess.run(cmd,shell=False)

def convert_names_on_iqtree()