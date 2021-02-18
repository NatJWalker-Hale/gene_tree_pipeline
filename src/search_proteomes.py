#! /usr/bin/python3

import sys
import os
import subprocess
import argparse
import shutil
from utils import parse_fasta

def fasta_to_stockholm(inf):
    nseq = 0
    with open(inf, "r") as seqf:
        for line in seqf:
            nseq += 1
    if nseq == 2:  # 1 sequence
        print("Single sequence, going straight to HMM")
        shutil.copy(inf, inf+".sto")
        # name .sto for compatibility, but actually .fa
    else:
        print("Aligning baits with FSA")
        cmd = ["fsa","--fast","--stockholm",inf]
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

def search_proteomes(bait,database_dir,output_dir):
    if "/" in bait:
        name = bait.split("/")[-1].split(".")[0]
    else:
        name = bait.split(".")[0]
    if not bait+".sto" in os.listdir(os.getcwd()):
        fasta_to_stockholm(bait)
        stockholm_to_profile(bait+".sto")
    else:
        try:
            stockholm_to_profile(bait+".sto")
        except subprocess.CalledProcessError:
            fasta_to_stockholm(bait)
            stockholm_to_profile(bait+".sto")
    dblist = []
    for dirpath,_,filenames in os.walk(database_dir):
        for f in filenames:
            if f.endswith(".pep.fa") or f.endswith(".cdhit"):
                dblist.append(os.path.abspath(os.path.join(dirpath,f)))
    for db in dblist:
        search_db(bait+".hmm",db)
        parse_search_out(db+".out")
        gather_sequences(bait,db+".hits",db,os.path.abspath(output_dir)+"/"+name+".hmmsearch.fa")
        os.remove(db+".out")
        os.remove(db+".hits")
    baitdict = dict([x for x in parse_fasta(bait)])
    with open(os.path.abspath(output_dir)+"/"+name+".hmmsearch.fa","a") as outf:
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
    search_proteomes(args.bait,args.database_dir,args.output_dir)  





