# SeqTools
Tools for NGS and third-generation sequencing (PacBio, Oxford Nanopore) data analysis.

This repository enumerates a set of basic python scripts used for sequencing analysis of NGS and third generation sequencing data like PacBio or ONT, as well as synthetic long reads (e.g. LoopSeq). Before running these scripts make sure you have all the dependencies installed! All of the tools were designed to run under UNIX with Python3 or later!

## Installation
Clone or download the repository.

```sh
git clone https://github.com/moldovannorbert/seqtools.git
cd seqtools
```

## Contents
- [MotifFinder](#MotifFinder)
- [ReadStatistics](#ReadStatistics)
- [TranscriptAnalyzer](#TranscriptAnalyzer)
- [BulkMapper](#BulkMapper)

## <a name="MotifFinder"></a>MotifFinder
MotifFinder is a tool for searching TATA, CAAT and GC-boxes as well as polyadenylation signals on a given reference genome using transcriptional start site (TSS) transcriptional end site (TES) or transcript coordinates. It will search for a specific pattern (perfect match of strings) found in `motif.txt` on the reference, and will report the sequence of the motif and the distance of the motif from the TSS or TES. It will also reproduce the sequence surroundings of TSS and TES positions, which can be used for initiator region (INR) or terminator region analysis.

### Dependencies
Before running this script please install the following dependencies:
-	[pandas]

### Usage
You will need:
- A `.fasta` file containing the reference genome.
- A directory containing `.gff3` file(s). The third column can contain either the text `tss` or `tes`. In these cases, the fourth and fifth columns must be the same, as they represent single nucleotide positions. The third column can also be `mRNA`, for which the following two columns must be different numbers. In this case the script assumes TES and TSS positions based on the fourth and fifth columns. 
The orientation (+/- in the seventh column) is always taken into account.

To see the help menu type `-h`.
```sh
./motiffinder.py –h
```

This line will create a directory containing the analysis.

```sh
./motiffinder.py –g /path/to/genome.fasta –i /path/to/annotation_directory/ -s /path/to/output_directory/ 
```

To change the number of nucleotides reproduced for the TSS use the `-ts` option. The default is: 5.
This line will reproduce a total of 21 letters from upstream and downstream of the TSS:
```sh
./motiffinder.py –g /path/to/genome.fasta –i /path/to/annotation_directory/ -s /path/to/output_directory/ -ts 10
```

To change the number of nucleotides reproduced for the TES use the `-te` option. The default is: 50.
This line will reproduce a total of 21 letters from upstream and downstream of the TES:
```sh
./motiffinder.py –g /path/to/genome.fasta –i /path/to/annotation_directory/ -s /path/to/output_directory/ -te 10
```

If your `.gff3` file contains read counts in integer format in the sixth column, you can filter your analysis by setting a minimum read count. To filter for read count us the `-mrc` option.
This line will output features with 10+ reads:
```sh
./motiffinder.py –g /path/to/genome.fasta –i /path/to/annotation_directory/ -s /path/to/output_directory/ -mrc 10
```

To exclude one contig from your analysis use the `-e` option followed by the contig identifier.
This line will exclude lines with the contig identifier ‘ABCD1’ from the analysis:
```sh
./motiffinder.py –g /path/to/genome.fasta –i /path/to/annotation_directory/ -s /path/to/output_directory/ -e ABCD1
```

## <a name="ReadStatistics"></a>ReadStatistics
ReadStatistics is a tool for calculating and visualising descriptive statistics of mapped reads in bulk.

### Dependencies
Before running this script please install the following dependencies:
-	[pandas]
-	[seaborn]

### Usage
You will need:
- A directory containing `.sam` files.

To see the help menu type ‘-h’.
```sh
./readstats.py –h
```
This line will create directories containing data files, the statistical analysis and figures:
```sh
./readstats.py /path/to/sam_directory/ /path/to/output_directory/
```

To exclude reads of one or more contigs from the analysis use the ‘-e’ option followed by a list of contigs separated by a coma.
This line will exclude reads mapping to contigs ABCD1, ABCD2 and ABCD3 from the analysis:
```sh
./readstats.py /path/to/sam_directory/ /path/to/output_directory/ -e ABCD1,ABCD2,ABCD3
```

If you already have the data files produced by the script and you want to rerun the statistics and plotting use the ‘-st’ option and replace the first path with the path to the folder containing the data files.
This line will use the DATA files in the specified directory to produce the output:
```sh
./readstats.py /path/to/DATA_directory/ /path/to/output_directory/ -st
```

## <a name="TranscriptAnalyzer"></a>TranscriptAnalyzer
TranscriptAnalyzer is a script for calculating the transcript length and exon count of transcripts from a ‘.gff3’ files and to determine the isoforms of transcripts. It can be used to determine differential transcript usage.

### Dependencies
Before running this script please install the following dependencies:
-	[pandas]

### Usage
You will need:
- A directory with one or more`.gff3` files containing only `mRNA` and `exon` lines. The files can rerpesent control data or experimental datasets. The control datasets have the `_control` suffix in their file name.

To see the help menu type ‘-h’.
```sh
./tranalyzer.py –h
```
This line will create directories containing per-transcript-per-sample reand lengths and exon counts and a statistical summary calculated from these:
```sh
./tranalyzer.py -i /path/to/input_directory/ -o /path/to/output_directory/
``` 
To determine differential transcript usage, use the `-c` option. You will need at least one contorl file and one experimental file for this in your input directory!
This line will output additional information about the type of a transcript in the experimental file compared to the control:
```sh
./tranalyzer.py -i /path/to/input_directory/ -o /path/to/output_directory/ -c
``` 
To change the default suffix for the control file use the `-s` option. The default is `_control`.
This line will change the default suffix of the control files for which the script is looking for in the input directory:
```sh
./tranalyzer.py -i /path/to/input_directory/ -o /path/to/output_directory/ -s _suffix
``` 

[pandas]: https://pandas.pydata.org/pandas-docs/stable/install.html
[seaborn]: https://seaborn.pydata.org/installing.html
