## Usage
**First step:**
Find reads carrying either heptamer or nonamer (both forward and reverse)

`exact_filter.sh <input.fastq> <filtered.fastq>`

**Second step:**
Find reads carrying both heptamer and nonmaer, and their distance is close to either 12-bp or 23-bp.
The output of this step is a list of read names that potentially carrying a recombination signal sequences (RSS).

`find_rss.py -f <filtered.fastq> -o <rss.names>`
