SOUFIR Emma M2BI - 2023/2024
Projet court - Assignation et détection des parties transmembranaires d'une protéine

# TM_detect

TM_detect is a programe which locates the membrane in a transmbrane protein and detects transmebrane segments. 

## Dependencies

This program uses the tool DSSP. You can install it : 


## Usage

```TM_detect [-h] [-n N] [-w WIDTH] [-g GAP] [-m GAP_MEMBRANE] filename
positional arguments:
  filename         A PDB file

optional arguments:
  -h, --help       show this help message and exit
  -n N             Number of points to place on the sphere. (default is 15)
  -w WIDTH         Initial width of the membrane. (default is 14 A)
  -g GAP           Gap of sliding membrane along an axis. (default is 1 A)
  -m GAP_MEMBRANE  Gap of optimising membrane's width. (default is 1 A)
```


## Algorithm


## Structure

This program is implemented using the OOP paradigm. 
