_Updated: September. 7, 2020_
## AIRRCall pipeline

Usage:
```
./AIRRCall.sh
```

The `AIRRCall.sh` pipeline uses the capture based short reads and the alleles downloaded from IMGT database to
- Call alleles
- Call novel alleles
- Call flanking sequences

### Call alleles

Shell script:
```
./allele_calling.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
```

The pipeline finds if the known allele is possessed by the person or not.
`./allele_calling.sh` first aligns the capture-based short reads to the IMGT alleles. Afterward, it analyzed the alignment, filtering out those mismatched reads and those reads with coverage length below a threshold (100 or allele length in the code). 
For each allele, a histogram on all positions of the allele is built, the coverage area of all filtered alleles are accumulated in the histogram. The minimum value in the histogram (the mimnum filtered read coverage of the allele) is the calling score of the allele.
The scoring of positive alleles are way larger than those of negative alleles.

### Call novel alleles

Shell script:
```
./novel_allele.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
```

### Call flanking sequences

Shell script:
```
./flanking_sequence.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
```


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

The above commands are integrated into `haplotype_analysis.sh`
```
../../immunogenomics/scripts/haplotype_analysis.sh IGKV2-29*01
```
- the bash file first use `fetch_sam.py` searching the allele sam file to find the target allele `IGKV2-29*01`'s alignment location to the assembly genome. 
- `samtools view` is then used to fetch all the alignment reads in the `IGKV2-29*01` alignment region
- the contig is also fetched by `fetch_chromosome.py`
- finally with the allele_SEQ, contig_SEQ and all the reads in the region, a haplotype analysis if performed with `parse_sam_haplotyping.py`
- the same analysis is applied on the other strand of the assembly genome

## Cluster alleles using Clustal-omega

Clustal-omega pre-built binary files can be downloaded from: http://www.clustal.org/omega/#Documentation

The following command clusters TCR alleles, and limit each cluster to include less than 20 alleles:

```
~/bin/clustal-omega-1.2.3-macosx -i TR_all.fasta -o TR_msa.fasta --clustering-out=TR_all.cluster --threads 4 --cluster-size 20
```

## Assembly approach
```
./assembly_analysis.sh
```
The command runs the assembly pipeline.
First `output_reads.py` is called and produce three files each cluster according to `all_top600_allelelen_ig.cluster.pickle`. Three files are:
- reads_cluster_H1.fasta
- reads_cluster_H2.fasta
- alleles_cluster.fasta
all the pair-end counterpart of reads in the cluster are also fetched and seperated in the H1 and H2 fasta files.
`reads_cluster_H1.fasta` and `reads_cluster_H2.fasta` are then made into contigs by `spades.py`.
We align `alleles_cluster.fasta` to the contigs with BWA.
`parse_bwa_sam.py` can filter the mismatch and soft-clip alignment results of BWA, the remaining alignments are all perfect matches.
All the matches from different clusters are incorporated into assembly_call.txt

### Checking the SPAdes coverage on annotations
```
python3 annotation_locus_parser.py -fa $file -fo annotation_contigs.txt
python3 locus.py -fs1 ./flanking_region_analysis/flanking_NA12878_tcrv_sup_size_H1.sam -fs2 ./flanking_region_analysis/flanking_NA12878_tcrv_sup_size_H2.sam -fo1 flanking_region_analysis/locus_flanking_NA12878_tcrv_sup_size_H1.pickle -fo2 flanking_region_analysis/locus_flanking_NA12878_tcrv_sup_size_H2.pickle > flanking_region_analysis/NA12878_tcrv_size_edit_dis.txt -td 10
python3 flanking_coverage.py -fna ./flanking_region_analysis/annotation_contigs_H2.txt -fpa ./flanking_region_analysis/locus_annotated_H2.pickle -fpf ./flanking_region_analysis/locus_flanking_H2.pickle
```
`annotation_locus_parser.py` is used to parse the annotation files in ./20200527_NA12878_BCRV_annotated_alleles 
Since the file name of the annotation indicates the contig of the annotated allele belong. `annotation_locus_parser.py`
incorporate the information of contig_name and allele position on the contig.
`locus.py` parses the H1, H2 sam file together and indicate the regions on asm_contig that SPAdes_contig covered. If the same SPAdes_contig align to H1 with significantly less mismatches than align to H2 (with difference larger than threshold `td`), locus.py keep only SPAdes_contig to H1 region and vice versa.
`flanking_coverage.py` compares the annotated allele positions and the SPAdes_contig covered regions. 

```
./alternative_contig.sh
python3 parse_contig_realign.py -fs NA12878_tcrv_support_asm/asm_realign/TCRV_realign_225.sam -fo NA12878_tcrv_support_asm/asm_realign/TCRV_remain_225.fasta > TCRV_realign_225.rpt
```
`alternative_contig.sh` script automatically run:
- the realignment of supporting reads to first round contig.
- parse_contig_realign.py to prepare the reads for alternative assembly.
- SPAdes that assemble the alternative contig.
- cat first round contig and the alternative contig into a fasta file.
- align the fasta file back to H1 and H2 assembly to check the correctness of the flanking region.

the `parse_contig_realign.py` parse the reads-to-SPAdes_contig realignment file:
- analyze the sam file to mark the potential variant (hot spot region).
- pop out the reads covered the hot spot region and support all variant favors the first round contig.
- produce the pair-end reads fasta file `TCRV_remain_225_P1.fasta` and `TCRV_remain_225_P2.fasta` that can be assembled into alternative contig.


## Novel alleles calling pipeline
```
./novel_alleles.sh
```
The input path of `novel_allele.sh` is list below
- outer_dir="NA12878_tcrv_novel_alleles/"   # the pipeline will make a directory to store all the files
- allele_path="TCRV_alleles.fasta"          # the known allele list download from IMGT database
- read_path_1="NA12878_S46_R1.fasta"        # captured R1 short read fasta file
- read_path_2="NA12878_S46_R2.fasta"        # captured R2 short read fasta file
- asm_path_H1="../NA12878/NA12878-H1.fa"    # FOR VERIFICATION ONLY -- assembly H1 fasta file
- asm_path_H2="../NA12878/NA12878-H2.fa"    # FOR VERIFICATION ONLY -- assembly H2 fasta file

The pipeline start with BWA alignment using `TCRV_alleles.fasta` as reference and short reads `NA12878_S46_R1/R2.fasta` as query to produce `bwa_read_to_allele.sam`.

The python file `parse_cluster_realign.py` analyze the variants in `bwa_read_to_allele.sam` and call the haplotype with reference to `TCRV_alleles.fasta`. Without the loss of generality we take a reference sequence TCRV allele TRAV21\*01 as an example. There is a G/A variant in TRAV21\*01, which means that there are several short reads support the haplotype with variant A other than G in original alleles in `TCRV_alleles.fasta`. We put all the haplotypes as candidates of novel alleles. The output fasta file is `corrected_alleles_raw.fasta`.

Since the correction is done on alleles separately, duplications may happened. For example, the corrected TRAV21\*01_corrected is in fact TRAV21\*02. The duplications may happened between `corrected_alleles_raw.fasta` and `TCRV_alleles.fasta` or between `corrected_alleles_raw.fasta` itself. `filter_corrected_alleles.py` filters the `corrected_alleles_raw.fasta` to `corrected_alleles_filtered.fasta`. In some cases, the duplicated novel alleles are longer than its original allele like the case TRAV21\*01_corrected is longer than TRAV21\*02, in these cases, we rename TRAV21\*01_corrected to TRAV21\*02_txtend. 

### Verification of novel alleles
This part of the pipeline is only done for verification if the long reads assembly exists. We first align the novel alleles to whole genome `NA12878-H1.fa` and `NA12878-H2.fa`. The novel allele that matched perfectly to the whole genome assembly can be seem as high confident novel allele. On the other hand, if an allele cannot matched to both `NA12878-H1.fa` and `NA12878-H2.fa`. It is probable that the novel alleles called had some problems.

`verify_novel_alleles.py` parsed the sam file from above BWA alignment to print the results.



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
