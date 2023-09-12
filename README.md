SOUFIR Emma M2BI - 2023/2024  
Projet court - Assignation et détection des parties transmembranaires d'une protéine

# TM_detect

TM_detect is a program which locates the membrane in a transmbrane protein and detects transmebrane segments. 

## Installation

### Dependencies

This program uses the tool DSSP combine with BioPython. You can install it : 

`sudo apt-get install dssp`

This program also uses PyMol. 

### Virtual environment

You can use the virtual environment located [Link text Here](https://link-url-here.org). 

## Usage

```
TM_detect [-h] [-n N] [-w WIDTH] [-g GAP] [-m GAP_MEMBRANE] filename

positional arguments:
  filename         A PDB file

optional arguments:
  -h, --help       show this help message and exit
  -n N             Number of points to place on the sphere. (default is 15)
  -w WIDTH         Initial width of the membrane. (default is 14 A)
  -g GAP           Gap of sliding membrane along an axis. (default is 1 A)
  -m GAP_MEMBRANE  Gap of optimising membrane's width. (default is 1 A)
```


## Structure
The source code is located in the `./src/` folder. 
This program is implemented using the OOP paradigm. 
The files [Protein.py](https://github.com/esoufir/TM_detect/blob/main/src/Protein.py) and [Geometry.py](https://github.com/esoufir/TM_detect/blob/main/src/Geometry.py) contain the classes. 
The file [TM_detect.py](https://github.com/esoufir/TM_detect/blob/main/src/TM_detect.py) contains the main instructions of the program. 
