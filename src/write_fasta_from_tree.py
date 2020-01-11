#! /usr/bin/python3

import sys 
import os
import argparse
from utils import parse_fasta
from newick3 import *

def write_fasta_from_tree(allfa,tree):
    