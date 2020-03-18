_Updated: Mar. 18, 2020_
## BLASTn-based pipeline

BLASTn performs local alignment, which effectively compares local regions of a sequence with the database.
Here we use reads as queries and build databases for TCR/BCR alleles of interest.
Because our reads can be longer than an allele, and we haven't had a deep understanding to the genetic diversity of the flanking regions, local-alignment-based methods can be useful in this application.

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
python3.5 ../../immunogenomics/scripts/parse_blastn_output.py -f NA12878_S46_L001_R1_001.blastn.out.txt -n 80 -o trbv_top80.txt
python3.5 ../../immunogenomics/scripts/compare_lists.py -fc trbv_top80.txt -fa trbv_all.txt -g TRBV
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
