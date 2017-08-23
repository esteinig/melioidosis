# Meliodosis

Nanopore sequencing and assembly of *B. pseudomallei* reference genome B03 from Balimo, Western Province, PNG.

This summary report is the initial test run for quality control and assembly of B03.

## Quality Control

Description:

Workflow:

Programs:

- nanoplot
- porechop
- filtlong
- minimap2
- samtools

Scripts:

- read_length_identity.py
- plot_identity_length.py
- plot_coverage.py

Configuration for QC:

```json
{
  "nanoplot": {
    "cpu": 32,
    "readtype": "1D",
    "color": "blue",
    "drop_outliers": true
  },
  
  "porechop": {
    "cpu": 32,
    "adapter_threshold": 90.0,
    "check_reads": 20000,
    "scoring_scheme": "3,-6,-5,-2"
  },
  
  "filtlong": {
    "other": "",
    "target_bases": 500000000,
    "min_length": 1000,
    "keep_percent": 90.0,
    "trim": true,
    "split": 250,
    "assembly": "ref/B03.fasta"
  },
  
  "minimap2": {
    "cpu": 32,
    "preset": "map10k",
    "reference": "ref/B03.fasta"
  }
}
```

Read summary **before filtering**:

```
Number of reads:	151097
Total bases:	962143036
Median read length:	3041.0
Mean read length:	6367.72
Readlength N50:	8152

Top 5 read lengths and their average basecall quality score:
Length: 482719bp	Q: 9.14
Length: 470466bp	Q: 6.88
Length: 451096bp	Q: 8.08
Length: 342367bp	Q: 6.9
Length: 336161bp	Q: 6.89

Top 5 average basecall quality scores and their read lengths:
Length: 988bp	Q: 15.61
Length: 826bp	Q: 15.08
Length: 503bp	Q: 14.98
Length: 4726bp	Q: 14.94
Length: 898bp	Q: 14.89

Number of reads and fraction above quality cutoffs:
Q5:	150912	99.88%
Q10:	102560	67.88%
Q15:	2	0.0%
Q20:	0	0.0%
```

Read summary **after filtering**:

```
Number of reads:	38105
Total bases:	500017486
Median read length:	10842.0
Mean read length:	13122.1
Readlength N50:	15819

Top 5 read lengths and their average basecall quality score:
Length: 82811bp	Q: 12.51
Length: 81737bp	Q: 11.75
Length: 80480bp	Q: 11.16
Length: 79874bp	Q: 11.56
Length: 79319bp	Q: 12.38

Top 5 average basecall quality scores and their read lengths:
Length: 4712bp	Q: 15.0
Length: 12325bp	Q: 14.78
Length: 2546bp	Q: 14.72
Length: 5104bp	Q: 14.65
Length: 9244bp	Q: 14.65

Number of reads and fraction above quality cutoffs:
Q5:	38105	100.0%
Q10:	38089	99.96%
Q15:	0	0.0%

```

Read quality (Q) vs. length (bp) **before filtering** (outliers removed):

Read quality (Q) vs. length (bp) **after filtering** (outliers removed):

Read identity (%) vs. length (bp) **before filtering**:

Read identity (%) vs. length (bp) **after filtering**:

Mean coverage against reference chromosomes **after filtering**:

Median coverage against reference chromosomes **after filtering**:
