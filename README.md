# Gene Tree Pipeline 

Python3 reimplementation of Ya Yang's bait homolog search pipeline from [Lopez-Nieves et al 2017](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.14822). Many of the scripts come courtesy of Ya Yang and Stephen Smith and are sourced from [here](https://bitbucket.org/yangya/adh_2016/src/master/). The core logic remains the same but various changes have been implemented: hmmer3 searches, no specification of a number of hits to preserve, continuous iteration until tip number has stabilised. 

I intend this pipeline to be the platform for some extensions in homology searching using various different tools. Usage (at your own risk) is totally fine, but please cite the above paper and the associated papers of any tools used, e.g. MAFFT, FSA, IQ-TREE, etc. a list of which is forthcoming.  
