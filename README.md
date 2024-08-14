# InDel_Toolkit

## Overview

Useful scripts for dealing with InDel Data. 

## Features

- stickleback: A python script for mapping long-reads sequencing data with engineered insertion sequences.
- smelt: A python script for parsing and counting insertions and deletions of defined size in sam files. 

## Installation

### Prerequisites

- Python 3.0+

```bash
pip install numpy
pip install pandas
pip install time
pip install re
pip install multiprocessing
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

### Stickleback Example

### DelMapper Example

DelMapper takes an SAM file, reads the [cigar strings](https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/), and identifies reads with deletions of a given size relative to the reference. It identifies these deletions on the translated sequence relative to a translated nucelotide reference to resolve ambiguities, and is therefore specifically designed for engineered libraries where specific codons are deleted. 

```bash
python DelMapper_v0.1.py < inputfile.sam (str)> < delSize (int)> < outputfile.csv (str)> 
```

To identify reads with 9 bp deletions from an example SAM:

```bash
python DelMapper.py output_test.sam 9 /nsPs_9D_PTD.csv

```

### Contributing

Explain how others can contribute to the project. Include guidelines for submitting pull requests, reporting issues, or suggesting features.

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
