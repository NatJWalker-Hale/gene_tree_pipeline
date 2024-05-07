#! /usr/bin/python3


import sys
import argparse
import phylo3
import newick3
from utils import parse_fasta


def get_taxon_code_name_corres(inf: str) -> dict:
    taxon_dict = {}
    with open(inf, "r", encoding='utf-8') as tax_table:
        for line in tax_table:
            code, name = line.strip().split("\t")
            taxon_dict[code] = name
    return taxon_dict


def relabel_tree(root: phylo3.Node, taxon_dict: dict) -> phylo3.Node:
    for n in root.leaves():
        try:
            code, seq = n.label.split("@")
        except ValueError:
            continue  # don't rename e.g. baits
        try:
            taxon = taxon_dict[code]
            print(f"{n.label}\t{'@'.join([taxon, seq])}")
            n.label = "@".join([taxon, seq])
        except KeyError:
            sys.stderr.write(
                f"code {code} not in taxon dictionary\n"
            )
            continue
    return root


def relabel_fasta(seq_dict: dict, taxon_dict: dict) -> dict:
    new_dict = {}
    for header, seq in seq_dict.items():
        try:
            code, sid = header.split("@")
            try:
                taxon = taxon_dict[code]
                print(f"{header}\t{'@'.join([taxon, sid])}")
                new_dict["@".join([taxon, sid])] = seq
            except KeyError:
                sys.stderr.write(
                    f"code {code} not in taxon dictionary\n"
                )
                new_dict[header] = seq
        except ValueError:
            new_dict[header] = seq
    return new_dict


def detect_input(inf) -> bool:
    with open(inf, "r", encoding="utf-8") as f:
        first = f.readline()
        if first.startswith(">"):
            return True
        if first.startswith("("):
            return False
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
        fasta = dict(parse_fasta(args.input))
        relabelled_fasta = relabel_fasta(fasta, taxon_corres)
        # print(relabelled_fasta)
        with open(args.input + ".name", "w", encoding="utf-8") as outf:
            for k, v in relabelled_fasta.items():
                outf.write(">" + k + "\n")
                outf.write(v + "\n")
    else:
        t = newick3.parse_from_file(args.input)
        relabel_tree(t, taxon_corres)  # modifies in place
        with open(args.input + ".name", "w", encoding="utf-8") as outf:
            outf.write(newick3.to_string(t) + "\n")
