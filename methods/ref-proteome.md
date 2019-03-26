

## Transcriptome

A TruSeq Stranded mRNA library (Illumina) was created using geoduck larvae (5 dpf) and sequenced on NovaSeq Platform (2 lanes). Sequence reads were assembled using Trinity.


```
Trinity \
--seqType fq \
--SS_lib_type RF \
--left NR021_S8_L001_R1_001.fastq,NR021_S8_L002_R1_001.fastq \
--right NR021_S8_L001_R2_001.fastq,NR021_S8_L002_R2_001.fastq \
--CPU 50 --trimmomatic --max_memory 500G
```

[output](http://owl.fish.washington.edu/halfshell/bu-mox/analyses/0804_1818/trinity_out_dir/0804_Pgen_larvae.fasta)

```
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	145131
Total trinity transcripts:	219698
Percent GC: 36.30

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 2651
	Contig N20: 1765
	Contig N30: 1263
	Contig N40: 921
	Contig N50: 676

	Median contig length: 331
	Average contig: 533.78
	Total assembled bases: 117271312


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 2384
	Contig N20: 1541
	Contig N30: 1090
	Contig N40: 791
	Contig N50: 592

	Median contig length: 324
	Average contig: 500.45
	Total assembled bases: 72630419
```

  ### Deduced Proteome

  Transdecoder was run on the assembled transcriptome


```
  TransDecoder.LongOrfs \
-t trinity_out_dir/0804_Pgen_larvae.fasta
```

The longest putative open reading frames were identified

```
/TransDecoder-4.0.0/TransDecoder.LongOrfs \
-t $fa \
-S
```

and BlastP performed

```
blastp \
-query PgenlarS-longest_orfs.pep \
-db uniprot_sprot_080917 \
-out blastout \
-num_threads 28 \
-evalue 1E-05 \
-max_target_seqs 1 \
-outfmt 6
```

with likely ORFs determined

```
TransDecoder-4.0.0/TransDecoder.Predict \
-t $fa \
--retain_blastp_hits blastout
```


Resulting peptide sequences for candidate ORFs with all shorter candidates within longer ORFs removed:

[0804_Pgen_larvae.fasta.transdecoder.pep](http://owl.fish.washington.edu/halfshell/bu-alanine-wd/17-08-10/0804_Pgen_larvae.fasta.transdecoder.pep)

Number of sequences: 36891
