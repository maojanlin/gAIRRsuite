## Graphs

Visualization of results from different stages of the pipeline.

With the `compare_with_annotation.py` -o option, we can generate a `trbv_top600.labelled.txt` file.
```
python ../../immunogenomics/scripts/compare_with_annotation.py -fc all_top600_allelelen.txt -fa NA12878_annotated_all.txt -g TR -o trbv_top600.labelled.txt
```
With the `trbv_top600.labelled.txt` annotation and `/scripts/plotting.Rmd` notebook, we can plot the scoring-ranking graph with annotated alleles and false positives in different color. Below graphs mentioned are our analysis on NA12878 TCRV and BCRV alleles.

`parse_blastn_output.py`
In mid-April, our pipeline is using the accumulated bit_score of BLASTn align all reads as query to all alleles as database.
If one read is aligned to multiple alleles with the same bit_score, we average the bit_score for all the aligned alleles.
The results are shown in `accumulate_average_tcrv.png` and `accumulate_average_bcrv.png`.

`coverage_analysis.py`
In mid-May, we analyze the read-depth coverage on each alleles with identity alignment. We notice that the minimum read-depth of each alleles differentiate the annotated alleles and unannotated ones. In TCRV alleles, we can get 100% accuracy and precision with the new pipeline; however, in BCRV alleles, we only get some improvement. Some false positive and false negative are still in our analysis pipeline. The results are shown in `min_read_depth_tcrv.png` and `min_read_depth_bcrv.png`.
