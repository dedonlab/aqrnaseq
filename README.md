# AQRNA-seq
This repo contains python processing scripts for AQRNA-seq data analysis. The scripts were used in conjuction with open-source tools/algorithms such as fastxtoolkit and blast as described below.

# Description of the AQRNA-seq data processing pipeline
Forward and reverse reads were trimmed of their respective adapter sequences using fastxtoolkit/0.013. A minimum adapter alignment length of 10bp was required, and unknown (N) nucleotides were kept: `fastx_clipper -Q33 -n -M 10`

Trimmed sequencing files from the BCG and sequencing control samples  were were blasted against a reference library using blast/2.6.0 with the following parameters:
```
blast -perc_identity 90 -word_size 9 -dust no -soft_masking false
```
For sequencing libraries prepared from BCG tRNA, a reference library was created based on the 48 entries in the gtRNA database for Myocbacterium bovis BCG str. Pasteur_1173P2 (http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/bacteria/Myco_bovi_BCG_Pasteur_1173P2_BCG_Pasteur_1173P2/Myco_bovi_BCG_Pasteur_1173P2_BCG_Pasteur_1173P2-gene-list.html). Sequences corresponding to duplicate tRNA genes (e.g. tRNA-Ala-TGC-1-2 and tRNA-Ile-GAT-1-2) and tRNA pseudogenes (tRNA-Ser-CGA-2-1) were removed to eliminate redundant entries and reduce the incidence of ambiguous or false positive matches. The terminal (3') CCA sequence was added to the tRNA sequence when it.
sequence for the 80-nt RNA oligo internal standard was added to the reference library. For the sequencing control samples, the sequences of the 5 synthetic RNA oligos was used to create a reference library, along with the 80-nt RNA oligo internal standard.

For each sample forward and reverse reads were merged by integrating their start and end positions to generate new start and end positions that reflect their combined coverage (blastpair.py). Multiple alignments were reduced by ranking all the alignments for a given read by their e-value and retaining only the alignment with the lowest e-value. Forward and reverse reads were required to match the same target. Paired reads that did not match the same target were discarded. Paired reads that mapped to multiple targets were stored in a separate file and not analyzed (cull.py) Uniquely mapped reads were tabulated and counted.

For the microRNA samples, the set of 963 miRNA sequences contained within the miRXplore Universal Reference product was combined with the sequence for the 80-nt RNA oligo standard was added to generate the reference library. In that analysis, 
the number of exact sequencing reads that matched to the reference miRNA sequence in each trimmed sequencing file was determined with `fgrep`:
```
numberReadsPerFile=$(fgrep $miRNA_sequence $trimmed_sequencing_file | wc -l)
```
The read counts of the miRNAs were normalized to the summed counts for all detected miRNAs to obtain a “normalized read count”. The summed counts of all detected miRNAs was also divided by 963, the total number of detected miRNAs, to obtain the “expected read count” assuming all species were equimolar. The read ratio was calculated by dividing the normalized read count by the expected read count.
