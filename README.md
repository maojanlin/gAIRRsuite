_Updated: Nov 16, 2025_
# gAIRR-suite: Profiling genes encoding the adaptive immune receptor repertoire
gAIRR-annotate provides annotations of IG and TR genes on personal assembly data

gAIRR-call can genotype IG and TR genes for short read data sequenced from germline cell. It is designed for probe captured gAIRR-seq data, but it can also be applied to other type of short read enriched from IG and TR regions.

- The database materials are updated to the latest (by 2025/11/02) IMGT/GENE-DB.
- Single-end read input is allowed in this version
- Beta option: use the novel alleles called from HPRC year 1 released assembles by gAIRR-annotate.  Since in this callset, some allele's core region is substring of another allele's, gAIRR suite may call both while only one of them is real. 


## Prerequisite programs:
- BWA aligner (0.7.17) 
- SPAdes assembler (v3.13.0)


## Installation
pip install gAIRR-suite
- [pip](https://pypi.org/project/gAIRR-suite/)
```
pip install gAIRR-suite
```
- [Github](https://github.com/maojanlin/gAIRRsuite.git)
```
git clone https://github.com/maojanlin/gAIRRsuite.git
cd gAIRRsuite
```
Though optional, it is a good practice to install in a virtual environment to manage the dependancies:

```
python -m venv venv
source venv/bin/activate
```
Now a virtual environment (named venv) is activated:

```
python setup.py install
```



## gAIRR-annotate usage
```
$ gAIRR_annotate -wd <work_dir> -id <sample_id> -a1 <assembly_h1.fa> -a2 <assembly_h2.fa> 
```

The `-a2` argument is optional for diploid personal assemblies.

The final annotation report is `work_dir/sample_id/group_genes.1.bed`, and `work_dir/sample_id/group_genes.2.bed` if the second haplotype of the assembly is provided.

The novel alleles in the annotation will be maked with parentheses, indicating the edit distance of the novel allele to the documented genes. For example, `TRAV8-3*01(i:1)` means the annotated gene has $1$ edit-distance to TRAV8-3*01, `TRAV19*01(i:0,h:17)` has no edit-distance to original allele, but has $17$ bases being clipped.

The novel allele sequence can be found in `work_dir/sample_id/novel_sample_ID_geneLocus.fasta.1/2.bed`

If only IG or TR genes are prefered, option `-lc IG` or `-lc TR` can be specified.

Use `-hprc` argument to use the additional novel alleles from HPRC year 1 release samples.


## gAIRR-call usage
```
$ gAIRR_call -wd <work_dir> -id <sample_ID> -rd1 <read.R1.fastq.gz> -rd2 <read.R2.fastq.gz> -lc <TRV TRJ TRD IGV IGJ IGD>
```

The prefered IG/TR locus and V/D/J genes can be specified with `-lc` option.

The final calling report is `work_dir/sample_ID/gene-locus/gAIRR-call_report.rpt`, `gene-locus` is the targeted gene set, like TCRV, TCRJ, etc.

The novel allele sequence from the `gAIRR-call_report.rpt` can be found in `work_dir/sample_ID/gene-locus_novel/gene-locus_with_novel.fasta`

`--flanking` option can be specified to run the flanking sequence assembly algorithm. The flanking sequnece information can be found in `work_dir/sample_ID/gene-locus_flanking/flanking_result/flanking_haplotypes.fasta`

Note that the seriel numbers provided by the `gAIRR_call` and `gAIRR_annotate` are not necessary corresponding. When there are multiple novel alleles, the seril number starting from $0$ can be in different order.

Use `-hprc` argument to use the additional novel alleles from HPRC year 1 release samples.


### checking RSS

Usage:
```
./scripts/check_RSS.sh
```

The `check_RSS.sh` pipeline uses the RSS and separated heptamer and nonamer sequences downloaded from IMGT database to check if there are proper RSS pattern in the flanking sequences

To run the `check_RSS.sh` pipeline, BWA aligner should be installed.

The `check_RSS.sh` pipeline first align all the known IMGT RSS to the flanking sequences to check if there are identical or near-identical RSS pattern. The flanking sequences missing RSS are then recorded in `RSS_checking/first_scan/missing_RSS_HG002-part_TCRJ_first_scan.fasta` and passed to second scanning. The second scanning aligned heptamer and nonamer sequences separately to the flanking sequences and try to identify heptamer-nonamer pairs that resemble proper RSSs.

Generated files:
`RSS_checking/first_scan/database_HG002-part_TCRJ_first_scan.csv` is the RSS report file of `HG002-part`. It indicate if the RSS are known, novel or could not be found after first scanning.
`RSS_checking/second_scan/database_HG002-part_TCRJ_second_scan.csv` is the RSS report file of the flanking sequences that missed RSS in the first scanning.


### Data collection pipeline

Usage:
```
./scripts/database_collect.sh
./scripts/allele_consensus.sh
```

The `database_collect.sh` pipeline collects the novel and flanking sequence database into database files. The duplicated novel or flanking sequences will be collapsed into one. Taking TRV novel allele as an example, generated file `database_novel_TRV.tsv` indicates which samples possess which novel allele, and `database_novel_TRV.fasta` recorded the novel allele sequence.

For samples with multiple assembly. Consensus allele result can be get from `allele_consensus.sh` pipeline. Taking `database_novel_TRV.tsv` and `database_novel_TRV.fasta` as input, `allele_consensus.sh` will generate `database_novel_TRV_consensus.tsv` and `database_novel_TRV_consensus.fasta` as output according to `./example/samples/consensus_name_HGSVC.log`.

In `./example/samples/consensus_name_HGSVC.log`, terms are separated by space. The first term is the consensus name while the later terms indicate the samples' different assembly id.


## Example

The `gAIRR_suite/material/` directory contains IMGT allele sequences and RSS information.
The `example/samples/` containts two miniature samples. `HG002_part_gAIRR-seq_R1.fasta` and `HG002_part_gAIRR-seq_R2.fasta` are a small part of the pair-end gAIRR-seq reads sequenced from HG002. `HG002-S22-H1-000000F_1900000-2900000.fasta` is a genome assembly sequence extracted from (Garg, S. *et al*, 2021). The genome sequence is the 1900000:2900000 segment from the contig HG002-S22-H1 of HG002's maternal haplotype assembly.

In the example settings. Running 
```
$ gAIRR_call -wd target_call -id HG002-part -rd1 gAIRR_suite/example/samples/HG002_part_gAIRR-seq_R1.fastq.gz -rd2 gAIRR_suite/example/samples/HG002_part_gAIRR-seq_R2.fastq.gz
```
or ```./scripts/AIRRCall.sh```
will gAIRR-call the HG002's AIRR alleles based on `HG002_part_gAIRR-seq_R1.fasta` and `HG002_part_gAIRR-seq_R2.fasta`.

Running 
```
$ gAIRR_annotate -wd target_annotate -id HG002-part -a1 gAIRR_suite/example/samples/HG002-S22-H1-000000F_1900000-2900000.fasta
```
or ```./scripts/AIRRAnnotate.sh``` 
will gAIRR-annotate part of the HG002's genome assembly `HG002-S22-H1-000000F_1900000-2900000.fasta`. In `./scripts/AIRRAnnotate.sh` , several shell script commands are commented. The commented commands are the settings to gAIRR-annotate two phased assemblies while in the example is to gAIRR-annotate single strend genome assembly.




<!--
## Verification pipeline

`AIRRCall.sh` and `AIRRAnnotate.sh` can be operated independently; however, `AIRRVerify.sh` can only be operated after both `AIRRCall.sh` and `AIRRAnnotate.sh` be performed successfully.

Usage:
```
./scripts/AIRRVerify.sh
```

`AIRRVerify.sh` check if the calling results of `AIRRCall.sh` and `AIRRAnnotate.sh` agree with each others.

To run the `AIRRVerify.sh` pipeline, BWA aligner should be installed.

The path parameters in `AIRRVerify.sh` should agree with those specified in `AIRRCall.sh`:
- workspace: the working directory of `AIRRCall.sh` (e.g. "./target_call/").
- allele_name: allele type (e.g. "TCRV").
- allele_path: where IMGT allele fasta file store (e.g. "../IMGT_alleles/TCRV_alleles.fasta").
- person_name: person's id (e.g. "NA12878").
And `AIRRAnnotate.sh`:
- asm_path_H1: personal assembly H1 contig fasta file (e.g. "../asm_NA12878/NA12878-H1.fa").
- asm_path_H2: personal assembly H2 contig fasta file (e.g. "../asm_NA12878/NA12878-H2.fa").
- annotation_dir: the working directory of `AIRRAnnotate.sh` (e.g. "./target_annotation/").

Generated files:

`target_call/NA12878_TCRV_flanking/verification/annotation_report.txt` shows the comparison details of each AIRRCall flanking sequence to AIRRAnnotate alleles including the position on personlized assembly and the extending length of AIRRCall flanking sequence.

`target_call/NA12878_TCRV_flanking/verification/annotation_summary.txt` shows all the false positive (AIRRCall only) and false negative (AIRRAnnotate only) of AIRRCall called flanking sequence comparing to AIRRAnnotate. For any AIRRCall flanking sequence that cannot be aligned to both H1 and H2 assembly fasta file, the flanking sequence is reported as redundant flanking sequence.


### Verification of novel alleles
This part of the pipeline is only done for verification if the long reads assembly exists. We first align the novel alleles to whole genome `NA12878-H1.fa` and `NA12878-H2.fa`. The novel allele that matched perfectly to the whole genome assembly can be seem as high confident novel allele. On the other hand, if an allele cannot matched to both `NA12878-H1.fa` and `NA12878-H2.fa`. It is probable that the novel alleles called had some problems.

`verify_novel_alleles.py` parsed the sam file from above BWA alignment and print the results.




Old method that use two round assembly approach
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
`annotation_locus_parser.py` is used to parse the annotation files in ./20200527_NA12878_BCRV_annotated_alleles. 
Since the file name of the annotation indicates the contig of the annotated allele belong, `annotation_locus_parser.py`
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
-->

