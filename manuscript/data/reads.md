## Read Statistics with NanoStats

#### Uncorrected

```

File: B03_uncorrected.fastq

General summary:

Number of reads:        151097
Total bases:            962143036
Median read length:     3041.0
Mean read length:       6367.7
Read length N50:        13500
Mean read quality:      10.7
Median read quality:    11.2

Top 5 longest reads and their mean basecall quality score:

1:      482719 (9.1)
2:      470466 (6.9)
3:      451096 (8.1)
4:      342367 (6.9)
5:      336161 (6.9)

Top 5 highest mean basecall quality scores and their read lengths:

1:      15.6 (988)
2:      15.1 (826)
3:      15.0 (503)
4:      14.9 (4726)
5:      14.9 (898)

Number and percentage of reads above quality cutoffs:

>Q5:    150912 (99.9%)
>Q10:   102560 (67.9%)
>Q15:   2 (0.0%)
>Q20:   0 (0.0%)
>Q25:   0 (0.0%)

```

### Corrected without Reference Genome

```

File: B03_uncorrected.fastq

Remove adapters with Porechop using --discard-middle for later compatibility with nanopolish index:

porechop -i B03_uncorrected.fastq -o B03_chopped.fastq --discard-middle

start:    53,396 / 151,097 reads      1,425,798 bp removed
end:      25,387 / 151,097 reads        346,282 bp removed
middle:       15 / 151,097 reads

File: B03_chopped.fastq

Filter reads without using reference data:

filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 B03_chopped.fastq > B03_corrected.fastq


General summary:

Number of reads:        24526
Total bases:            500011337
Median read length:     17883.5
Mean read length:       20387.0
Read length N50:        21447
Mean read quality:      11.9
Median read quality:    12.1

Top 5 longest reads and their mean basecall quality score:

1:      482719 (9.1)
2:      470466 (6.9)
3:      451096 (8.1)
4:      342367 (6.9)
5:      336161 (6.9)

Top 5 highest mean basecall quality scores and their read lengths:

1:      14.8 (12328)
2:      14.6 (9247)
3:      14.4 (15245)
4:      14.3 (9829)
5:      14.3 (9807)

Number and percentage of reads above quality cutoffs:

>Q5:    24526 (100.0%)
>Q10:   22521 (91.8%)
>Q15:   0 (0.0%)
>Q20:   0 (0.0%)
>Q25:   0 (0.0%)



```
