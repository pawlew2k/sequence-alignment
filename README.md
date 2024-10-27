# Sequence alignment

## Description

Implementation of 2 sequence alignment algorithms used in bioinformatics to align protein or nucleotide sequences:

- [Needleman–Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) finds the optimal local alignment
- [Smith–Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) finds the optimal global alignment

Program designates N alignments with the highest score based on substitution matrix and gap penalty.

## Installation

### Prerequisites
- [![python](https://upload.wikimedia.org/wikipedia/commons/1/16/Blue_Python_3.10%2B_Shield_Badge.svg)](https://www.python.org)

- Python packages:
[![NumPy](https://img.shields.io/badge/numpy-%23013243.svg?style=for-the-badge&logo=numpy&logoColor=white)](https://numpy.org/)
[![Pandas](https://img.shields.io/badge/pandas-%23150458.svg?style=for-the-badge&logo=pandas&logoColor=white)](https://pandas.pydata.org/)

### Setup

1. **Clone the repository**:
``` bash
git clone https://github.com/yourusername/sequence-alignment.git
cd sequence-alignment
```

1. **Install dependencies**:
``` bash
pip install -r requirements.txt
```

## Usage

### Run

```
python needleman_wunsch.py|smith_waterman.py
  -h, --help            show this help message and exit
  --seq1 SEQ1           first sequence
  --seq2 SEQ2           second sequence
  -n N                  count of maximum optimal alignments
  -s SUBSTITUTION_MATRIX, --substitution_matrix SUBSTITUTION_MATRIX
                        filepath to the substitution matrix in CSV format
  -g GAP_PENALTY, --gap_penalty GAP_PENALTY
                        value of constant penalty for insertion a gap ('-') incurs
  -o OUTPUT, --output OUTPUT
                        output filename, file will be saved in '/output' directory```
```
### Example

**Needleman–Wunsch algorithm**

``` python
python needleman_wunsch.py --seq1 TATA --seq2 ATAT -n 3 -s substitution_matrix.csv -o global.txt
```

Output:

```
Global alignment no. 1:
TATA-
-ATAT
Score: 11

Global alignment no. 2:
-TATA
ATAT-
Score: 11
```

**Smith–Waterman algorithm**

``` python
python smith_waterman.py --seq1 TATA --seq2 ATAT -n 3 -s substitution_matrix.csv -o local.txt
```

Output:

```
Local alignment no. 1:
TATA
TAT-
Score: 13

Local alignment no. 2:
ATA-
ATAT
Score: 13
```

### File Structure

```
sequence-alignment/
├── needleman_wunsch.py             # Main script for running Needleman–Wunsch algorithm via CLI
├── smith_waterman.py               # Main script for running Smith–Waterman algorithm via CLI
│
├── substitution_matrix.csv         # Substitution scoring matrix used by algoritmhs
│
├── output/                         # Directory to save the results from the main scripts
│   ├── global.txt               
│   ├── local.txt
│
├── src/
│   ├── directional_cell.py         # Defines cell with directions used in alignment matrices
│   ├── parser.py                   # Parses input parameters
│   ├── sequence_alignment.py       # Implements Needleman-Wunsch and Smith-Waterman algorithms;
│                                   # Includes alignment logic, score matrix creation, traceback and output handling       
├── tests/
│   ├── unit_tests.py               # Unit tests
│
├── task_content.pdf                # Project description and requirements 
├── requirements.txt                # Dependencies required to run the project
└── README.md                       # Project description and usage instructions
```
