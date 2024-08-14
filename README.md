# InDel_Toolkit

## Overview

Useful scripts for dealing with InDel Data. 

## Features

- stickleback: A python script for mapping long-reads sequencing data with engineered insertion sequences.
- smelt: A python script for parsing and counting insertions and deletions of defined size in sam files. 

## Installation

### Prerequisites

- Python 3.0+

- And python packages, via *pip*:
  
```bash
pip install numpy
pip install pandas
pip install multiprocessing
pip install Bio
pip install levenshtein
```

### Installing

1. Clone the repository:
    ```bash
    git clone https://github.com/QVEU/InDel_Toolkit.git
    ```
2. Navigate to the project directory:
    ```bash
    cd InDel_Toolkit
    ```
    
## Command-Line Usage

### *stickleback* Example

*stickleback* was design to map long-read sequencing libraries from nanopore or pacbio containing engineered libraries containing any defined insertion introduced via the [SPINE](https://academic.oup.com/nar/article/48/2/e11/5634037) pipeline.

```bash
python stickleback.py <pathto/input.sam (str)> <query sequence (str)> </path/to/templateFasta (str)> [Min Read Length (int)] [Max Read Length (int)]

# <> = req'd argument, []= optional argument

```

To identify reads with 9 bp deletions from an example SAM:

```bash
> cd InDel_Toolkit/
> python lib/stickleback.py test/stickleback_test.sam agcgggagaccggggtctctgagcg lib/templates/puc19-ev71-twtainan1998_4643-bsmbi-and-bsai-free-deleted-1-annotations-1-7471.fasta 


----------------=============-----------------
--==--==--==--==   ><```º>   ==--==--==--==--=
==--==--==--==-- stickleback --==--==--==--==-
----------------=============-----------------


Query Length: 25
Template Length: 7471
Template-Query Distance: 9
Mapping Reads of Size 25 to 7471 with a cutoff of 7

Loading SAM: test/stickleback_test.sam
Total Candidate Reads: 250

Outfile: test/stickleback_test_stickleback.csv

1. Computing minimum distance hit position for 247 reads.
Query: AGCGGGAGACCGGGGTCTCTGAGCG
2. Mapping minimum distance site on template sequence...
done.

Mapped hits in 36 reads.
Done in 0.03501284917195638 minutes.
Wrote test/stickleback_test_stickleback.csv.

```

### *smelt* Example

*smelt* compliments Stickleback, and is used to tabulate *deletions* across coding sequences. It takes an SAM file (from Illumina Sequencing, or other high-accuracy NGS methods), reads the [cigar strings](https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/), and identifies reads with deletions of a given size relative to the reference. It identifies these deletions on the translated sequence relative to a translated nucelotide reference to resolve ambiguities, and is therefore specifically designed for engineered libraries where specific codons are deleted. 

```bash
python smelt.py <inputfile.sam (str)> <Deletion Size (int)> <outputfile.csv (str)>

# <> = req'd argument

```

To identify reads with 3 bp deletions from an example SAM:

```bash
cd InDel_Toolkit/lib/
python smelt.py ../test/smelt_test.sam 3 ../test/smelt_test.csv
```

### Referencing
Please cite: Bakhache, W., Orr, W., McCormick, L. & Dolan, P., Uncovering Structural Plasticity of Enterovirus A through Deep Insertional and Deletional Scanning. Research Square (Preprint) (2024).

### Contributing

Steps to contribute: 

- Fork the repository.
- Create a new branch (git checkout -b feature-branch).
- Commit your changes (git commit -am 'Add new feature').
- Push to the branch (git push origin feature-branch).
- Open a pull request.

### License

This project is licensed under the MIT License - see the LICENSE.md file for details.

### Contact

Email: Patrick.Dolan@nih.gov

Twitter: @drptdolan

GitHub Issues: Submit an issue

### Acknowledgements

The authors are thankful for minimap2.
