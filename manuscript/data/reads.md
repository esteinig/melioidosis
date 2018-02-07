## Read Statistics with NanoStats

#### Uncorrected

```

File: B03_uncorrected.fastq

General Statistics

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





```
