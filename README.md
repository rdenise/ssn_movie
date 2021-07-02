# ssn_movie

Create a image by image visualisation of a SSN based on user annotation

# Requirement

python == 3.8
graphviz == 2.47.2

## Python libraries

networkxgmml==0.1.6
networkx==2.5.1
matplotlib==3.4.2
pandas==1.2.4
numpy==1.20.3

# Help

```
usage: ssn_movie.py [-h] -g <XGMML> [-k <KOFAM>] [-e <EGGNOG>] [-a <annotation>] [-o <OUTPUT>]

See the effect of alignement threshold based on annotation color

optional arguments:
  -h, --help            show this help message and exit

General input dataset options:
  -g <XGMML>, --xgmml <XGMML>
                        XGMML file from the analysis of EFI-EST : https://efi.igb.illinois.edu/efi-est/
  -k <KOFAM>, --kofam <KOFAM>
                        KOFAM file from the analysis of kofamscan
  -e <EGGNOG>, --eggnog <EGGNOG>
                        EGGNOG file from the analysis of eggnog-mapper
  -a <annotation>, --annotation <annotation>
                        Tabulated file with your own annotation, need to have a header with the columns 'Hit_Id' and 'Gene' where hit_id id the name of the sequence in your network and gene is the annotation on the sequence
  -o <OUTPUT>, --output <OUTPUT>
                        Name of the output file (default: [NAME_OF_XGMML])
```

# Input annotation file format

For the annotation file (that do not come from EGGNOG or KOFAMScan) the file need to be a tabulated separated value file as follow

The first column contains the name of the sequence as present in the fasta file that was used to create the SSN
The second column contains the annotation for this sequence

```
#Hit_Id  Gene
zzzz0001  pdxK
ARM001  pdxX
```

