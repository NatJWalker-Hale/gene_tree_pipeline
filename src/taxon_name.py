#! /usr/bin/python3


import sys
import argparse
import phylo3
import newick3
from utils import parse_fasta


def get_taxon_code_name_corres(inf: str) -> dict:
    taxon_dict = {}
    with open(inf, "r") as tax_table:
        for line in tax_table:
            code, name = line.strip().split("\t")
            taxon_dict[code] = name
    return taxon_dict


def relabel_tree(root: phylo3.Node, taxon_dict: dict) -> phylo3.Node:
    for n in root.leaves():
        code, seq = n.label.split("@")
        try:
            taxon = taxon_dict[code]
            print(n.label + "\t" + "@".join([taxon, seq]))
            n.label = "@".join([taxon, seq])
        except KeyError:
            sys.stderr.write(
                "code %s not in taxon dictionary\n" % code
            )
            continue
    return root


def relabel_fasta(seq_dict: dict, taxon_dict: dict) -> dict:
    new_dict = {}
    for k, v in seq_dict.items():
        code, seq = k.split("@")
        try:
            taxon = taxon_dict[code]
            print(k + "\t" + "@".join([taxon, seq]))
            new_dict["@".join([taxon, seq])] = v
        except KeyError:
            sys.stderr.write(
                "code %s not in taxon dictionary\n" % code
            )
            new_dict[k] = v
    return new_dict


def detect_input(inf) -> bool:
    with open(inf, "r") as f:
        first = f.readline()
        if first.startswith(">"):
            return True
        elif first.startswith("("):
            return False
        else:
            sys.stderr.write("input not recognised\n")
            sys.exit()


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input tree (newick) or sequence \
                        (FASTA) to rename. Will autodetect")
    parser.add_argument("taxon_table", help="code-to-taxon correspondences, \
                        tab-separated, one-per-line")
    args = parser.parse_args()

    taxon_corres = get_taxon_code_name_corres(args.taxon_table)
    # print(taxon_corres)

    if detect_input(args.input):
        fasta = dict([x for x in parse_fasta(args.input)])
        relabelled_fasta = relabel_fasta(fasta, taxon_corres)
        with open(args.input + ".name", "w") as outf:
            for k, v in relabelled_fasta.items():
                outf.write(">" + k + "\n")
                outf.write(v + "\n")
    else:
        t = newick3.parse_from_file(args.input)
        relabel_tree(t, taxon_corres)  # modifies in place
        with open(args.input + ".name", "w") as outf:
            outf.write(newick3.to_string(t) + "\n")
