_Updated: May 5, 2021_
## gAIRR-call

Usage:
```
./scripts/AIRRCall.sh
```

The `AIRRCall.sh` pipeline uses the capture-based short reads and the alleles downloaded from IMGT database to
- Find **novel alleles**
- **Call alleles** (including both known and novel alleles)
- Assemble and haplotype **flanking sequences**

To run the `AIRRCall.sh` pipeline, BWA aligner and SPAdes assembler should be installed.

The path parameters in `AIRRCall.sh` should be specified:
- workspace: the directory all results and intermediate data be stored (e.g. "./target_call/").
- path_SPAdes: the path to the spades.py installed (e.g. "../SPAdes-3.11.1-Linux/bin/spades.py")

The below three parameters indicate the interested allele reference (IMGT) fasta files.
- list_allele_name: target allele types (e.g. "TCRV TCRJ BCRV")
- allele_dir: the directory IMGT allele fasta file store (e.g. "../IMGT_alleles/")
- allele_suffix: the suffix of allele fasta file, should agree with the real file name ( e.g."\_alleles.fasta")

The final three parameters indicate the target sequencing fasta files.
- person_name: person's id (e.g. "NA12878").
- read_path_1: capture-based short reads R1 fasta file (e.g. "NA12878_S46_R1.fasta").
- read_path_2: capture-based short reads R2 fasta file (e.g. "NA12878_S46_R2.fasta").

Below paths in the README.md use `NA24385` and `TCRV` as examples.

### Find novel alleles

Shell script:
```
./scripts/novel_allele.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
```

The `novel_allele.sh` pipeline aligns capture-based short reads to IMGT alleles with BWA MEM. Then the program `parse_cluster_realign.py` finds the variant in each alleles. If threre are variants, the program haplotypes the allele and call the haplotypes not in the IMGT database as novel allele candidates.

Since the correction is done on each allele separately, two alleles may generate the same novel allele candidate. It is also possible that the generated novel allele is infact another known allele. For example, the corrected TRAV21\*01_corrected is the same as TRAV21\*02. To deal with duplication problems, `filter_corrected_alleles.py` filters the `corrected_alleles_raw.fasta` to `corrected_alleles_filtered.fasta`. Those duplicated sequence are elimitated. 

Generated files:

`target_call/NA24385_TCRV_novel/corrected_alleles_filtered.fasta` is all the novel allele candidates fasta file.
`target_call/NA24385_TCRV_novel/TCRV_with_novel.fasta` is the merged allele file containing IMGT alleles and haplotyped novel allele candidates.

### Call alleles

Shell script:
```
allele_path=target_call/NA24385_TCRV_novel/TCRV_with_novel.fasta
./scripts/allele_calling.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
```

The pipeline finds if the alleles in the merged fasta file (containing both IMGT and novel alleles) are positive (possessed by the person).
In `allele_calling.sh`, the capture-based short reads are aligned to the merged allele fasta file with BWA MEM again. This time the "-a" option of BWA is used to ensure that the alleles can be reached by any potential reads. Afterward, `analyze_read_depth_with_bwa.py` filters out reads with edit-distance (mismatches or indels) and reads with coverage length below a threshold (minimum between 100 or allele length in the code). 

For each allele, a histogram on all positions of the allele is built. The read-depth of all filtered alleles are accumulated in the histogram. The minimum value in the histogram (the mimnum filtered read depth of the allele) is the calling score of the allele.
Empirically, the scores (minimum read-depth) of the true alleles are way larger than those of false alleles.

Generated files:

`target_call/NA24385_TCRV/read_depth_calling_by_bwa.rpt` reports all the alleles sorted by their scores (minimum read-depth).
`target_call/NA24385_TCRV/gAIRR-call_report.rpt` reports the positive alleles with the adaptive threshold.
`target_call/NA24385_TCRV/allele_support_reads.pickle` is a pickle file containing a dictionary. The dictionary indicates the names of the read supporting each alleles. The dictionary key is the allele name and the dictionary value is a set containing all reads support (perfectly match with enough length coverage) the allele.

### Assemble and haplotype flanking sequences

Shell script:
```
./scripts/flanking_sequence.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2} ${path_SPAdes}
```

The `flanking_sequence.sh` first groups pair-end read sequences and allele sequences in the directory `target_call/NA24385_TCRV_flanking/group_allele_reads/` according to `target_call/NA24385_TCRV/allele_support_reads.pickle`. Then the sub-pipeline `denovo_backbone.sh` uses SPAdes to assemble each short reads group into an unphased flanking contig (backbone). Afterward, each allele has a backbone in the `target_call/NA24385_TCRV_flanking/asm_contigs/` directory.

The `denovo_backbone.sh` sub-pipeline also use BWA MEM to align alleles to the backbone to check the correctness of the backbone. Those backbones do not contain perfectly matched contig are discarded. The the start and end positions of the allele in qualified backbones are marked by `parse_bwa_sam.py`.

`target_call/NA24385_TCRV_flanking/flanking_result/flank_region.txt` is the start-end report and
`target_call/NA24385_TCRV_flanking/flanking_result/flanking_contigs.fasta` is the backbone collections.

Finally the `flanking_sequence.sh` align all the capture-based reads to the backbones `target_call/NA24385_TCRV_flanking/flanking_result/flanking_contigs.fasta`, the sam file is haplotyped by `shrink_sam_to_range.py`. The region extending 200 bps from two ends of the original allele is cropped and reported as the flanking sequences.

Generated file:

`target_call/NA24385_TCRV_flanking/flanking_result/flanking_haplotypes.fasta` is the final called flanking sequences.


## gAIRR-annotate

Usage:
```
./scripts/AIRRAnnotate.sh
```

The `AIRRAnnotate.sh` pipeline uses the personal assembly contigs and the alleles downloaded from IMGT database to
- **Call alleles** (novel alleles are marked)
- **Call flanking sequences**

To run the `AIRRAnnotate.sh` pipeline, BWA aligner should be installed.

The path parameters in `AIRRAnnotate.sh` should be specified:
- outer_dir: the directory all results and intermediate data be stored (e.g. "./target_annotation/")
- list_allele_name: target allele types (e.g. "TCRV TCRJ BCRV")
- allele_dir: the directory IMGT allele fasta file store (e.g. "../IMGT_alleles/")
- allele_suffix: the suffix of allele fasta file, should agree with the real file name ( e.g."\_alleles.fasta")
- person_name: person's id (e.g. "NA12878").
- asm_path_H1: personal assembly H1 contig fasta file (e.g. "../asm_NA12878/NA12878-H1.fa")
- asm_path_H2: personal assembly H2 contig fasta file (e.g. "../asm_NA12878/NA12878-H2.fa")

The `AIRRAnnotate.sh` pipeline first indexes the personalized assembly and aligns IMGT alleles to the assembly with BWT. Afterward, `annotation_with_asm.py` analyzes the alignment sam file. The perfectly matched alleles are kept and aligned alleles with edit-distance are seen as novel alleles. Finally `get_asm_flanking.py` utilize the previous alignment sam file and crop the flanking sequence from the personal assembly fasta file.

Generated files:

`target_annotation/annotation_imperfect_NA12878_TCRV.txt` is the report showing the aligned allele, aligned contig, contig position, and alignment length. If there are edit-distance in the alignment, the report shows additional tag the same as sam format.
`target_annotation/novel_NA12878_TCRV.fasta` is the collection of novel alleles.
`target_annotation/flanking_NA12878_TCRV.fasta` is the collection of flanking sequence (including novel alleles).


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

The `database_collect.sh` pipeline collects the novel and flanking sequence database into database files. The duplicated novel or flanking sequences will be collapsed into one. Take TRV novel allele as an example, generated file `database_novel_TRV.tsv` indicates which samples possess which novel allele, and `database_novel_TRV.fasta` recorded the novel allele sequence.

For samples with multiple assembly. Consensus allele result can be get from `allele_consensus.sh` pipeline. Taking `database_novel_TRV.tsv` and `database_novel_TRV.fasta` as input, `allele_consensus.sh` will generate `database_novel_TRV_consensus.tsv` and `database_novel_TRV_consensus.fasta` as output according to `./example/samples/consensus_name_HGSVC.log`.

In `./example/samples/consensus_name_HGSVC.log`, terms are separated by space. The first term is the consensus name while the later terms indicate the samples' different assembly id.


## Example

The `example/material/` directory contains IMGT allele sequences and RSS information.
The `example/samples/` containts two miniature samples. `HG002_part_gAIRR-seq_R1.fasta` and `HG002_part_gAIRR-seq_R2.fasta` are a small part of the pair-end gAIRR-seq reads sequenced from HG002. `HG002-S22-H1-000000F_1900000-2900000.fasta` is a genome assembly sequence extracted from (Garg, S. *et al*, 2021). The genome sequence is the 1900000:2900000 segment from the contig HG002-S22-H1 of HG002's maternal haplotype assembly.

In the example settings. Running 
```./scripts/AIRRCall.sh```
will gAIRR-call the HG002's AIRR alleles based on `HG002_part_gAIRR-seq_R1.fasta` and `HG002_part_gAIRR-seq_R2.fasta`.
Running 
```./scripts/AIRRAnnotate.sh``` 
will gAIRR-annotate part of the HG002's genome assembly `HG002-S22-H1-000000F_1900000-2900000.fasta`. In `./scripts/AIRRAnnotate.sh` , several shell script commands are commented. The commented commands are the settings to gAIRR-annotate two phased assemblies while in the example is to gAIRR-annotate single strend genome assembly.




<!--
## Verification pipeline

`AIRRCall.sh` and `AIRRAnnotate.sh` can be operated independently; however, `AIRRVerify.sh` can only be operated after both `AIRRCall.sh` and `AIRRAnnotate.sh` be performed successfully.

Usage
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

