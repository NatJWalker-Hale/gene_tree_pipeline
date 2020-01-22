#! /usr/bin/python3

import sys
import os
import subprocess
import argparse

def fasta_to_aln(inf,thread=2,alner="mafft"):
    if alner == "mafft":
        cmd = ["mafft","--auto","--amino","--thread",str(thread),inf]
        outf = open(inf+".mafft.aln","w")
    elif alner == "fsa":
        cmd = ["fsa","--fast",inf]
        outf = open(inf+".fsa.aln","w")
    print("Aligning sequences with "+alner)
    proc = subprocess.Popen(cmd,shell=False,stdout=outf,stderr=subprocess.PIPE)
    while proc.poll() is None:
        print(proc.stderr.readline().decode().strip())

def clean(aln):
    cmd = ["pxclsq","-s",aln,"-o",aln+"-cln","-p","0.1"]
    print("Cleaning alignment")
    subprocess.run(cmd,shell=False)

def aln_to_tree(aln,treeblder="fasttree",thread=2):
    if treeblder == "fasttree":
        cmd = ["fasttree","-wag","-out",aln+".fasttree.tre",aln]
    elif treeblder == "iqtree":
        cmd = ["iqtree","-s",aln,"-m","WAG+G","-nt",str(thread)]
    print("Inferring tree with "+treeblder)
    proc = subprocess.Popen(cmd,shell=False,stdout=subprocess.PIPE)
    proc.wait()
    if treeblder == "iqtree":
        os.rename(aln+".treefile",aln+".iqtree.tre")

def fasta_to_tree(inf,thread=2,alner="mafft",treeblder="fasttree"):
    if not inf+"."+alner+".aln" in os.listdir(os.getcwd()):
        fasta_to_aln(inf,thread,alner)
    if not inf+"."+alner+".aln-cln" in os.listdir(os.getcwd()):
        try:
            clean(inf+"."+alner+".aln")
        except:
            fasta_to_aln(inf,thread,alner)
            clean(inf+"."+alner+".aln")
    aln_to_tree(inf+"."+alner+".aln-cln",treeblder,thread)

if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-a","--aligner",help="The software used to infer alignments. Options are mafft (auto, default), fsa",default="mafft")
    parser.add_argument("-t","--tree_builder",help="The software used to infer trees. Options are fasttree (default), iqtree",default="fasttree")
    parser.add_argument("-nt","--threads",help="The number of threads to use for alignment and tree inference (default 2)",default=2)
    parser.add_argument("fasta", help="FASTA of sequences to align and infer a tree")
    args = parser.parse_args()

    fasta_to_tree(args.fasta,args.threads,args.aligner,args.tree_builder)