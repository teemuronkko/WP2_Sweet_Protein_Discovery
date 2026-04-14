# WP2_SweetProteinDiscovery

GitHub repository for WP2 on Sweet Proteins. 

## 1. Data sources
We used two complementary sources for building our databases of plant proteins. Firstly, all unreviewed and reviewed proteins under taxonomy Viridiplantae were retrieved from UniProt (date) [1]. Secondly, we accessed the data from One Thousand Plant Transcriptomes Initiative (OneKP) [2]. The data from OneKP was extracted from [CyVerse](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019) and the public [Google Drive](https://drive.google.com/drive/folders/18AOvneP_1l5uzE7tVWPKVR9MAkhrtA2N) containing some additional data. All proteins that were longer than 1000 amino acids from both sources were removed from the input data. This was done to speed up ProstT5 3Di sequence generation.

## 2. Foldseek database creation
To search for structurally similar proteins to known plant proteins, we created a custom Foldseek database using ProstT5 3Di structural sequences. These 3Di structural sequence use a custom alphabet of letters to encode for specific structural features of every given amino acid in the input sequence. Therefore, by transforming our input sequences to 3Di sequences, we can use the created sequences to search for structurally similar proteins using an input protein (e.g. Thaumatin-I) in our database of 3Di sequences.

The Foldseek database was created using the following command:
```
foldseek createdb $INPUT_DB $DB_OUTPUT --prostt5-model $MODEL_WEIGHTS --gpu 1 --threads $SLURM_CPUS_PER_TASK --write-lookup 1
```

The above database creation was done in batches of 500,000 sequences using a mixture of NVIDIA A100 and Tesla V100 GPUs. The batched databases were concatenated using `foldseek concatdbs`, using the script [./foldseek_scripts/foldseek_concatdbs.sh](./foldseek_scripts/foldseek_concatdbs.sh).

The Foldseek query database was created using the experimental (Brazzein: 4HE7, MNEI: 1IV7, Thaumatin-I: 1LR2, and Thaumatin-II: 3AOK), or predicted (Miraculin: https://www.uniprot.org/uniprotkb/P13087/entry and Honey truffle sweetener: [3]) structures of the known sweet-tasting or taste-modifying proteins. The following command was used to create the query database:
```
foldseek createdb $INPUT_DB $DB_OUTPUT
```

## 3. Foldseek search 
The Foldseek search was performed using the following command:
```
foldseek easy-search "$INPUT_DB" "$TARGET_DB" $OUTPUT $TMP_DIR \
	--exhaustive-search 1 \
	--min-ungapped-score 30 \
	--threads $SLURM_CPUS_PER_TASK \
	-e 0.05 \
	--alignment-type 2 \
	-c 0.50 \
	--cov-mode 2 \
	--format-output "query,target,evalue,bits,pident,fident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,qseq,tseq,q3di,t3di,qheader,theader,qaln,taln,q3dialn,t3dialn,qcov,tcov" \
	--format-mode 4
```
The search was performed separately on the OneKP and UniProt plant protein data. 

## 4. Download AFDB structures for UniProt hits
As the protein on UniProt have been modelled using AF2 and the structures are publicly available on AFDB, we downloaded the structures from AFDB. In total, XXXX structures were downloaded. XXXX proteins did not have modelled structures on AFDB, and these strucutures were modelled using XXXX.

## 5. Model OneKP hits using OpenFold 3
The OneKP proteins do not have publicly available protein structures, and therefore, we had to model these using OpenFold 3 using Gefion. The proteins that did not have an AFDB entry from UniProt were also modelled using OpenFold 3. 

## 6. Structural analysis and comparison of `foldseek-search` results

## 7. Toxicity and allergenity of identified hits
To ensure that the identified hits are not toxic or allergenic, we used several tools to computationally assess the toxicity and allergenity of the identified hits. The following tools were used.
- TODO decide which tools to use!

## 8. Hit selection
After filtering out the hits with toxic and allergenic features, we selected the final set of proteins to study experimentally.

# References:
[1] UniProt: the universal protein knowledgebase in 2025. Nucleic acids research 53, no. D1 (2025): D609-D617.
[2] One Thousand Plant Transcriptomes Initiative, Gane Ka-Shu Wong, James H. Leebens-Mack (2019). Data packages for One Thousand Plant transcriptomes and phylogenomics of green plants. CyVerse Data Commons. DOI 10.25739/8m7t-4e85
[3] McFarland, Chase, Bassam Alkotaini, Chloe P. Cowen, Michael G. Edwards, Elizabeth Grein, Alan D. Hahn, Justine C. Jennings et al. "Discovery, Expression, and In Silico Safety Evaluation of Honey Truffle Sweetener, a Sweet Protein Derived from Mattirolomyces terfezioides and Produced by Heterologous Expression in Komagataella phaffii." Journal of agricultural and food chemistry 72, no. 35 (2024): 19470-19479.