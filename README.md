# ALIdetection

The following are the data files and sample code used in the ALI detection project.

(1) npbin.r: R code of NPBin method

(2) demo_npbin.r: demo using real ChIP-seq data.

(3) data atac.txt, data ctcf.txt, data dnase.txt: the pre-processed ChIP-seq data analyzed in this paper. 
The following are the column names and their meanings of these data files.
1. chr: chromosome
2. location: genomic location based on hg19
3. m: total number of reads covering the SNP
4. xm: total number of reads at the SNP from the maternal allele
5. winning.chip: the allele with more ChIP-seq reads. "P" if xm < m/2 and "M" otherwise.
6. motif: the ID and the transcription factor name of the motif in JASPAR database (Mathelier and others, 2013).
7. pval.mat.atSNP: the p-value of the motif on the maternal allele from R package atSNP (Zuo and others, 2015).
8. pval.pat.atSNP: the p-value of the motif on the paternal allele from atSNP.
9. pval.rank.atSNP: the p-value of the rank test of the allelic motif strength difference from atSNP.
10. winnig.motif: the allele with stronger motif, e.g., it is "M" if pval.mat.atSNP < pval.pat.atSNP.
11. potential TP: whether it is a potential TP based on our criteria described in the Supplementary notes. The users can define it differently using different thresholds on the various p-values from atSNP.
12. potential FP: whether it is a potential FP based on our criteria described in the Supplementary notes. The users can define it differently using different thresholds on the various p-values from atSNP.

Note that ff the motif related items of a SNP are all "NA", it means that there is no known
motif at this SNP with either pval.mat.atSNP or pval.pat.atSNP less than 0.01.
