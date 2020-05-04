_Updated: May. 3, 2020_
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

We then use `scripts/compare_lists.py` to compare BLASTn calls with annotation which uses the information from personalized reference genomes.

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
- set an identity threshold at 100
- writes reads associated with each allele group for second-step high-confidence allele calling
  - cluster labelling is from `TR_all.cluster`
  - cluster-allele-read information is shown as a dictionary, serialized by pickle, as `all_top600_allelelen_tr.cluster.pickle`

## Cluster alleles using Clustal-omega

Clustal-omega pre-built binary files can be downloaded from: http://www.clustal.org/omega/#Documentation

The following command clusters TCR alleles, and limit each cluster to include less than 20 alleles:

```
~/bin/clustal-omega-1.2.3-macosx -i TR_all.fasta -o TR_msa.fasta --clustering-out=TR_all.cluster --threads 4 --cluster-size 20
```


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
