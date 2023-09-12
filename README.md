SOUFIR Emma M2BI - 2023/2024  
Projet court - Assignation et détection des parties transmembranaires d'une protéine

# TM_detect

TM_detect is a program which locates the membrane in a transmbrane protein and detects transmebrane segments. 

## Installation

### Virtual environment

You can use the virtual environment [TM_detect.yml](https://github.com/esoufir/TM_detect/blob/main/src/TM_detect.yml). 

### Dependencies

This program uses the tool DSSP combine with BioPython. You can install it : 

`sudo apt-get install dssp`

This program also uses PyMol. Evrything is detailed in the virtual environment [TM_detect.yml](https://github.com/esoufir/TM_detect/blob/main/src/TM_detect.yml). 




## Usage

When in the `./src/` folder, run the following command : 
```
python TM_detect [-h] [-n N] [-w WIDTH] [-g GAP] [-m GAP_MEMBRANE] filename

positional arguments:
  filename         A PDB file

optional arguments:
  -h, --help       show this help message and exit
  -n N             Number of points to place on the sphere. (default is 15)
  -w WIDTH         Initial width of the membrane. (default is 14 A)
  -g GAP           Gap of sliding membrane along an axis. (default is 1 A)
  -m GAP_MEMBRANE  Gap of optimising membrane's width. (default is 1 A)
```
You can use some data examples located in the folder `./data/` : 

For instance, the following command : 
```
python TM_detect.py ../data/1k24.pdb -n 15 -w 14 -g 1 -m 1
```
gives the output : 


## Structure
The source code is located in the `./src/` folder. 
This program is implemented using the OOP paradigm. 
The files [Protein.py](https://github.com/esoufir/TM_detect/blob/main/src/Protein.py) and [Geometry.py](https://github.com/esoufir/TM_detect/blob/main/src/Geometry.py) contain the classes. 
The file [TM_detect.py](https://github.com/esoufir/TM_detect/blob/main/src/TM_detect.py) contains the main instructions of the program. 
