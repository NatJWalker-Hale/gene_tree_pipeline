#! /usr/bin/python3

import sys
import os
import subprocess
import argparse
from utils import parse_fasta

def fasta_to_stockholm(inf):
    print("Aligning baits with FSA")
    cmd = ["fsa","--stockholm",inf]
    out = subprocess.run(cmd,shell=False,capture_output=True,text=True).stdout
    with open(inf+".sto","w") as outf:
        outf.write(out) 

def stockholm_to_profile(sto):
    print("Building HMM profile")
    cmd = ["hmmbuild",sto[:-3]+"hmm",sto]
    subprocess.run(cmd,shell=False)

def search_db(hmm,dbf):
    print("Searching Database")
    cmd = ["hmmsearch","--noali","--tblout",dbf+".out",hmm,dbf]
    subprocess.run(cmd,shell=False)

def parse_search_out(outf):
    print("Parsing search output")
    hits = open(outf[:-3]+"hits","w")
    with open(outf,"r") as search_out:
        for line in search_out:
            if not line.startswith("#"):
                hits.write(line.split(" ")[0].strip()+"\n")
    hits.close()

def gather_sequences(inf,hitsf,dbf,outf):
    print("Compiling sequences")
    seqout = open(outf,"a")
    seql = open(hitsf,"r").readlines()
    dbdict = dict([x for x in parse_fasta(dbf)])
    for s in seql:
        seqout.write(">"+s)
        seqout.write(dbdict[s.rstrip("\n")]+"\n")
    seqout.close()

def search_proteomes():
    if "/" in args.bait:
        name = args.bait.split("/")[-1].split(".")[0]
    else:
        name = args.bait.split(".")[0]
    fasta_to_stockholm(args.bait)
    stockholm_to_profile(args.bait+".sto")
    dblist = []
    for dirpath,_,filenames in os.walk(args.database_dir):
        for f in filenames:
            if f.endswith(".pep.fa") or f.endswith(".cdhit"):
                dblist.append(os.path.abspath(os.path.join(dirpath,f)))
    print(dblist)
    for db in dblist:
        search_db(args.bait+".hmm",db)
        parse_search_out(db+".out")
        gather_sequences(args.bait,db+".hits",db,os.path.abspath(args.output_dir)+"/"+name+".hmmsearch.fa")
        os.remove(db+".out")
        os.remove(db+".hits")
    baitdict = dict([x for x in parse_fasta(args.bait)])
    with open(os.path.abspath(args.output_dir)+"/"+name+".hmmsearch.fa","a") as outf:
        for key,value in baitdict.items():
            outf.write(">"+key+"\n")
            outf.write(value+"\n")

    
if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")
    
    parser = argparse.ArgumentParser()
    parser.add_argument("bait",help="FASTA file of query sequences")
    parser.add_argument("database_dir",help="Directory containing sequence databases to be searched")
    parser.add_argument("output_dir",help="Directory to put output")
    args = parser.parse_args()  
    search_proteomes()  





