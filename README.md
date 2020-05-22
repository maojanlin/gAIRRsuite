_Updated: May. 21, 2020_
## BLASTn-based pipeline

BLASTn performs local alignment, which effectively compares local regions of a sequence with a database.
Here we use reads as queries and build two databases, one for TCR genes and one for BCR genes.
Because our reads (300-bp, paired-end) are usually longer than BCR/TCR alleles, and we haven't had a deep understanding to the genetic diversity of the flanking regions, local-alignment-based methods is useful.

### Run BLASTn

Bash script is in `scripts/blastn.sh`. 
First, we build the BLASTn database for the target gene (e.g. `TRBV`).
Second, because BLASTn doesn't support FASTQ inputs, we must convert the reads to a FASTA file.
This step discards the sequencing quality information, so an additional quality check step might be helpful (we haven't done this).
Finally, we perform a BLASTn local alignment for all converted reads. This step is quite fast because currently our data is from capture-based targeted sequencing, which includes much less reads compared to a WGS dataset.

### Process BLASTn results and compare them with annotation

We use `scripts/parse_blastn_output.py` to process BLASTn alignment results.
By default, we use the _bit_score_ reported by the program and take the top hits (alleles).
If there is a tie, we split the score evenly between each hit.
(This may not be statistically true because bit score is in log space, but empirically the results look good.)
Top `n` (which is an argument `scripts/parse_blastn_output.py`) hits are reported.

We then use `scripts/compare_with_annotation.py` to compare BLASTn calls with annotation which uses the information from personalized reference genomes.

An example for internal use:

```
cd /home/naechyunchen/NAS/Yuchun/naechyun/blast/experiments
python ../../immunogenomics/scripts/parse_blastn_output.py -f NA12878_S46_full_name.blastn.out.txt -n 600 -o all_top600_allelelen.txt --fn_allele_len ../../../20200317_filtered_V_alleles_for_probe_design/allele_len.tsv --identity_thrsd 100 -c clustering/TR_all.cluster -ocp all_top600_allelelen_tr.cluster.pickle
python ../../immunogenomics/scripts/compare_with_annotation.py -fc all_top600_allelelen.txt -fa NA12878_annotated_all.txt -g TR
```

The above `parse_blastn_output.py` example:
- uses blastn log `NA12878_S46_full_name.blastn.out.txt`
- takes top-600 ranked alleles
- writes top-ranked alleles and associated scores as `all_top600_allelelen.txt`
- normalizes scores with allele length, extracted from file `allele_len.tsv`
- set an identity threshold at 100, it is highly recommended to set the identity threshold at 100.
  - Since in the second phase analysis, 100% identity is also considered, identity threshold 100 in this phase can group reads in different clusters more precisely.
- writes reads associated with each allele group for second-step high-confidence allele calling
  - cluster labelling is from `TR_all.cluster`
  - cluster-allele-read information is shown as a dictionary, serialized by pickle, as `all_top600_allelelen_tr.cluster.pickle`

## Detailed analysis of the read-coverage on candidate alleles

After `parse_blastn_output.py` and the clusters information group up reads and alleles into different clusters, an all read-to-allele comparison is done in each cluster.

In comparison, the coverage of a read on allele should more than a threshold (100), and no Hamming distance is allowed for a match. In other word, all the matched reads possess perfect identity to the allele with more than (100) nucleotide coverage.

After the comparison, we analyze the read-depth on every alleles, and another threshold (10) is posed on the minimun read-depth on alleles. Only the alleles with minimum read-depth more than (10) are classified as high-confidence allele. It turns out that most unannotated alleles' minimum read-depths are below (10), and annotated ones are above the threshold.

```
python ../../immunogenomics/scripts/coverage_analysis.py -fr NA12878_S46.fasta -fa clustering/IG_all.fasta -ans NA12878_annotated_all.txt -md 10 -fnp ./NA12878/all_top600_allelelen_ig_idthrd100.cluster.pickle -fop ./NA12878/BCR_top600_reads_alleles.pickle -fsp ./NA12878/NA12878_bcrv_support_reads.pickle > read_depth_NA12878_bcrv.log
```
The above `coverage_analysis.py` example:
- fetches the sequence data from read fasta file `NA12878_S46.fasta` and allele fasta file `TR_all.fasta`
- reference the clustering in `all_top600_allelelen_tr.cluster.pickle` to build read/allele pickle file `TCR_top600_reads_alleles.pickle`
- If `TCR_top600_reads_alleles.pickle` is already existed, the data will be directly loaded instead of re-calculating
- the minimum coverage threshold between read and allele is set to 100
- the minimum read-depth threshold is set to 10 and can be adjusted with -md command
- the -ans argument read in the annotation file for different sample
- the -fsp argument output a pickle file showing all the supporting reads of all alleles in the coverage analysis
- `coverage_analysis.py` also shows which allele and its min/average/max coverage on terminal when processing. And can be parsed to `read_depth_NA12878_bcrv.log`

### Checking the annotation

This part can only be done with annotation from whole genome assembly. After the coverage analysis, we can see whats going on in each allele from the `NA12878_bcrv_support_reads.pickle` file. Below two command generates the supporting reads fasta file `sup_reads-f1-H1.fasta` and `chromosome-f1-H1.fasta`
```
python ../../../immunogenomics/scripts/fetch_support_reads.py -fr NA12878_S46.fasta -fsup NA12878_bcrv_support_reads.pickle -fa IGHV1/OR16-3*01 -fo sup_reads-f1-H1.fasta
python ../../../immunogenomics/scripts/fetch_chromosome.py -fg ../../../../asm/NA12878/NA12878-H1.fa -fc NA12878-S1544-H1-000005F -fo chromosome-f1-H1.fasta
```
We can also check the whole reads alignment on the region, and check the number of supporting reads of different haplotype.

```
samtools view ../../../../align_to_asm/20200514_NA12878_aligntoasm_H1/NA12878_align_to_H1asm.sorted.bam "NA12878-S1544-H1-000005F:21185-21186" > interval.sam
python parse_sam_haplotyping.py -fs interval.sam -fp 21186 21198 21238 -fo haplotype.txt
```


## Cluster alleles using Clustal-omega

Clustal-omega pre-built binary files can be downloaded from: http://www.clustal.org/omega/#Documentation

The following command clusters TCR alleles, and limit each cluster to include less than 20 alleles:

```
~/bin/clustal-omega-1.2.3-macosx -i TR_all.fasta -o TR_msa.fasta --clustering-out=TR_all.cluster --threads 4 --cluster-size 20
```

# Worki-in-progress methods that consider RSS structures

## Usage
**Generate search patterns for haptamer and nonamer**

Process heptamer/nonamer sequences from database and generate pattern for `grep`.
The search pattern includes both forward and reverse sequences.

`python process_rss_for_grep.py -f <seq.txt>` 

**Find reads carrying either haptamer or nonamer using grep**

Find reads carrying either heptamer or nonamer (both forward and reverse).
Using `grep` for fast screening.

`sh exact_filter.sh <input.fastq> <filtered.fastq>`

**Find potential RSS by looking at the distance between heptamer and nonamer**

Find reads carrying both heptamer and nonmaer, and their distance is close to either 12-bp or 23-bp.
The tolerance of distance can be set by argument `-t` , which we recommend to be `0` since further analysis doesn't support `t > 0`.
The output of this script is a list of read names that potentially carrying a recombination signal sequences (RSS).

`python find_rss.py [-t INT] -f <filtered.fastq> -o <rss.names>`
