#! /usr/bin/python3


import os
import re
import sys
import argparse
from utils import parse_fasta


def get_seqs(q_ids: list, s_seqdict: dict) -> dict:
    out_dict = {}
    for s in q_ids:
        try:
            out_dict[s] = s_seqdict[s]
        except KeyError:
            try:
                out_dict[s] = s_seqdict[re.sub(".p$", "", s)]
            except KeyError:
                sys.stderr.write("unable to match %s, check\n" % s)
    return out_dict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("sequences", help="FASTA-formatted query sequences")
    parser.add_argument("database_dir", help="Directory containing coding \
                        sequence DBs to extract from. Names should match \
                        query IDs")
    args = parser.parse_args()

    q_dict = dict([x for x in parse_fasta(args.sequences)])
    q_seq_ids = [x for x in q_dict.keys()]
    q_taxa = set([x.split("@")[0] for x in q_seq_ids if "@" in x])
    # no "@" is likely bait, so we skip

    db_dict = {}  # key is file, value is path
    for dirpath, _, filenames in os.walk(args.database_dir):
        for f in filenames:
            if f.endswith(".cds.fa"):  # add suffix
                db_dict[f] = os.path.abspath(os.path.join(dirpath, f))
    # debug
    # print(db_dict)
    db_names = [re.sub(".cds.fa$", "", x) for x in db_dict.keys()]
    # debug
    # print(db_names)

    out_seqs = {}
    for t in q_taxa:  # iterate over taxa we have proteins for
        # get seq ids for that taxon
        q_taxon_ids = [x for x in q_seq_ids if x.split("@")[0] == t]
        # first, try to match file name directly
        if t in db_names:  # match
            taxon_db_path = db_dict[t + ".cds.fa"]  # get path
            cds_dict = dict([x for x in parse_fasta(taxon_db_path)])
            s_dict = get_seqs(q_taxon_ids, cds_dict)
            out_seqs = out_seqs | s_dict
        else:  # try to match by first header
            sys.stderr.write("searching headers for %s\n" % t)
            match = False
            for v in db_dict.values():
                with open(v, "r") as checkf:
                    sid = checkf.readline().split("@")[0].lstrip(">")
                    # debug
                    # print(sid)
                    # can't have redundant taxon labels...
                    if sid == t:
                        if not match:
                            match = True
                            cds_dict = dict([x for x in parse_fasta(v)])
                            s_dict = get_seqs(q_taxon_ids, cds_dict)
                            out_seqs = out_seqs | s_dict
                        else:
                            sys.stderr.write("taxon %s already matched, \
                                             check that labels are \
                                             non-redundant\n" % t)
                            sys.exit()

    for k, v in out_seqs.items():
        print(">" + k)
        print(v)
