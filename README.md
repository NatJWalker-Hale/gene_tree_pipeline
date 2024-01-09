# Gene Tree Pipeline 

Python3 reimplementation of Ya Yang's bait homolog search pipeline from [Lopez-Nieves et al 2018](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.14822). Many of the scripts come courtesy of Ya Yang and Stephen Smith and are sourced from [here](https://bitbucket.org/yangya/adh_2016/src/master/). The core logic remains the same but various changes have been implemented.

This pipeline is actively developed. Feel free to share and use (at your own risk), but please cite the above paper, and the relevant papers for any dependencies when you do.

## dependencies

- Python3 (tested on 3.10.6)
- HMMER (http://hmmer.org/, tested on v3.2.1)
- phyx (https://github.com/FePhyFoFum/phyx, all utilities in path as `px*`)
- FSA (https://fsa.sourceforge.net/, in path as `fsa`)
- mafft (https://mafft.cbrc.jp/alignment/software/, tested on v7.490, in path as `mafft`)
- FastTree2 (http://www.microbesonline.org/fasttree/, in path as `fasttree`)
- one of...
  - raxml-ng (https://github.com/amkozlov/raxml-ng, in path as `raxml-ng`, default)
  - IQ-TREE2 (https://github.com/iqtree/iqtree2, in path as `iqtree2`, needs changes)

## description

This pipeline is designed to use a set of query homologs ("baits") to search a set of proteomes (incl. some transcriptome-derived) for similar sequences, and then conduct multiple rounds of alignment, tree inference, long branch cleaning, monophyletic masking, and subtree pruning. The output is (a set of) subtree(s) containing the input baits and their homologs from the queried proteomes.

By default, the pipeline uses FSA (`--fast`) to align the bait sequences and produce a stockholm-formatted alignment which is used to create an HMM with HMMER and search each proteome with hmmsearch. Blastp can also be called for searching, in which case the full set of queries are used to query each proteome, hits with bitscores less than min_bitscore are filtered, and hits with bitscores less than a percentage threshold (default 10%) of the query's best hit are filtered. In both cases, the top k hits can optionally be retained (directly from the hmmsearch output, or post sorting the filtered hits, for blastp). The resulting hits are combined with the bait sequences and used to infer an alignment with mafft (`--auto`). The alignment is trimmed of sites containing less than 10% non-ambiguities (`pxclsq -p 0.1`) and a tree is inferred with FastTree (`-wag`). Long terminal branches (1.5, or 1.0 and > 10x sister length) are removed. If `-m` and `-mp` are specified, monophyletic masking will be used to compress monophyletic (and single-node paraphyletic) sets of transcripts from the same taxon to the longest transcript present in the cleaned alignment (useful for transcriptome sequences). `-if` can be used to give a file specifying taxa to ignore while masking (useful if a mix of genomes and transcriptomes). Finally, long internal branches (> 1.0) are cut to separate subtrees, and subtrees containing baits are preserved. If `-it` is > 1, this process is repeated on the resulting subtree(s) FASTA(s) (default 3 rounds).

TBA: rooted-ingroup subtree extraction (options `-po`, `-og`).

## required inputs

A FASTA-formatted file of amino acid sequences of the queries of interest, and a directory containing the proteomes to search, each with suffix `.pep.fa` or `.pep.fa.cdhit`. All sequence names should be formatted taxon_name@sequence code, e.g. `Beta@BVRB_01G0098`. Bait labels may be arbitrary.

## options

Run `python3 bait_homologs.py` to see a full list of command line options:

```
usage: bait_homologs.py [-h] [-b] [--min_bitscore MIN_BITSCORE] [--threshold THRESHOLD] [-a ALIGNER] [-t TREE_BUILDER] [-tc TIP_ABS_CUTOFF] [-rc TIP_REL_CUTOFF] [-ic INTERNAL_CUTOFF] [-mt MIN_TAXA]
                        [-nt THREADS] [-m] [-mp] [-if IGNORE_FILE] [-po] [-og OUTGROUPS] [-it ITERATE] [-o OUTPUT_DIR] [-k KEEP]
                        bait database_dir

positional arguments:
  bait                  FASTA file of baits to search
  database_dir          Path to the database containing proteomes to search. Expects file endings of .pep.fa or .cdhit

options:
  -h, --help            show this help message and exit
  -b, --blast           Use blastp for similarity search instead of default hmmsearch
  --min_bitscore MIN_BITSCORE
                        Filter blastp hits with bitscore lower than min_bitscore (default 30.0)
  --threshold THRESHOLD
                        Filter blastp hits with bitscore lower than threshold * max bitscore of query (default 0.1)
  -a ALIGNER, --aligner ALIGNER
                        Alignment software to use: mafft (auto, default), fsa (defaults to --fast)
  -t TREE_BUILDER, --tree_builder TREE_BUILDER
                        Tree building software to use: fasttree (wag, default), raxml-ng (defaults to WAG+G)
  -tc TIP_ABS_CUTOFF, --tip_abs_cutoff TIP_ABS_CUTOFF
                        Absolute branch length cutoff for trimming tips. Tips longer than this will be trimmed. Defaults to 1.5
  -rc TIP_REL_CUTOFF, --tip_rel_cutoff TIP_REL_CUTOFF
                        Relative branch length cutoff for trimming tips. Tips longer than this and at least 10x longer than sister will be trimmed. Defaults to 1.0
  -ic INTERNAL_CUTOFF, --internal_cutoff INTERNAL_CUTOFF
                        Branch length cutoff for internal branches. Subtrees subtended by branches longer than this will be trimmed. Defaults to 1.0
  -mt MIN_TAXA, --min_taxa MIN_TAXA
                        Minimum taxa in a subtree to conserve and check for bait presence. Defaults to 4
  -nt THREADS, --threads THREADS
                        Number of threads to use. Defaults to 2
  -m, --mask            If this flag is selected, monophyletic masking will be conducted
  -mp, --mask_paraphyly
                        Whether to mask paraphyletic sequences while doing monophyletic masking
  -if IGNORE_FILE, --ignore_file IGNORE_FILE
                        File containing taxon names (one per line) to ignore while masking monophyletic tips. Defaults masks all taxa
  -po, --prune_og       Whether to extract rooted ingroup clades containing baits after first round of tree inference (requires OG file, default off)
  -og OUTGROUPS, --outgroups OUTGROUPS
                        File containing outgroup taxon labels, one per line
  -it ITERATE, --iterate ITERATE
                        how many times to iterate tree building and cleaning
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to put output. Defaults to current directory
  -k KEEP, --keep KEEP  Number of hits to keep (default all)
```

## example

Files for an example run can be found in [/example/](https://github.com/NatJWalker-Hale/gene_tree_pipeline/tree/master/example). The search directory is `/example/pep/`. For transcriptome-derived sequences, redundant transcripts are first clustered using CD-HIT (https://sites.google.com/view/cd-hit, `-n 5 -c 0.995`), leaving files labelled `.pep.fa.cdhit`, but this is optional. Genome-derived proteomes are not clustered, therefore these remain `.pep.fa`. `baits.pep.fa` contains two characterised DODA homologs from _Beta vulgaris_. Because two of our proteomes come from genome annotations, we include them in `ignore_file.txt` as `Athaliana` and `Beta`.

We conduct a run, using defaults and activating monophyletic masking, with the following command:

```
python3 ../src/bait_homologs.py -m -mp -if ignore_file.txt bait.pep.fa pep/
```

### output

This generates a lot of files, but we'll explore the main ones:

- `log.txt` contains a full record of the run, including input command.
- `bait.pep.fa.sto` contains the stockholm-formatted pairwise alignment of the two baits. You can view this with e.g. belvu.
- `bait.pep.fa.hmm` contains the HMM from hmmbuild.
- `bait.hmmsearch.fa` contains the amalgamated sequences from searching each of the protein databases with hmmsearch, plus the baits.
- `bait.hmmsearch.fa.mafft.aln` contains the mafft alignment of `bait.hmmsearch.fa`.
- `bait.hmmsearch.fa.mafft.aln-cln` is the mafft alignment cleaned of columns with less than 10% data with `pxclsq`.
- `bait.hmmsearch.fa.mafft.aln-cln.fasttree.tre` is the FastTree inference on the cleaned alignment, with `-wag`.
- `bait.hmmsearch.fa.mafft.aln-cln.fasttree.tre.tt` is the same tree after trimming tips.
- `bait.hmmsearch.fa.mafft.aln-cln.fasttree.tre.tt.mm` is the tip-trimmed tree after monophyletic masking.
- `bait_1.subtree` is the subtree containing the baits resulting from internal branch cutting on tip-trimmed, monophyletically-masked tree, and `bait_1.pep.fa` is the corresponding FASTA.

From here, the pipeline proceeds through two (by default) further iterations, appending _1, _2, etc. to the subtrees each time. In our example, there is only one subtree containing the baits, so we end up with bait_1_1_1.subtree and bait_1_1_1.pep.fa, containing the final subtree and final FASTA after three iterations.

## miscellaneous options

Most scripts can also be used standalone - for example `search_proteomes.py` can be used as a general wrapper for hmmsearch or blastp searching of a specified proteome(s). For any subscript, see the available options by running e.g. `python3 search_proteomes.py`.

Resulting trees and FASTAs can be automatically renamed from codes to any other name by including a `taxon_table` file where each line is tab-separated code and corresponding name. The script `taxon_name.py` can be used as follows:

```
python ../src/taxon_name.py bait_1_1_1.pep.fa taxon_table
```

An example `taxon_table` file is included in `/example/`.
