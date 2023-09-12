SOUFIR Emma M2BI - 2023/2024  
Projet court - Assignation et détection des parties transmembranaires d'une protéine

# TM_detect

TM_detect is a program which locates the membrane in a transmbrane protein and detects transmebrane segments. 

## Installation

### Dependencies

This program uses the tool DSSP combine with BioPython. You can install it : 

`sudo apt-get install dssp`

This program also uses PyMol. Evreything is detailed in the virtual environment [TM_detect.yml](https://github.com/esoufir/TM_detect/blob/main/src/TM_detect.yml). 


Download the ZIP repository.   
In a proper folder, unzip it.   
Change your current directory into the root of the project.   
```
conda env create --prefix ./mypymolenv --file src/TM_detect.yml
conda activate mypymolenv
```



### Virtual environment

You can use the virtual environment [TM_detect.yml](https://github.com/esoufir/TM_detect/blob/main/src/TM_detect.yml). 


## Usage

After activating the virtual environment and when in the `./src/` folder, run the following command : 
```
python TM_detect [-h] [-n N] [-w WIDTH] [-g GAP] [-m GAP_MEMBRANE] path/to/filename

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
python TM_detect.py ../data/1prn.pdb -n 15 -w 14.0 -g 1.0 -m 1.0
```
gives the output : 
![Example](https://github.com/esoufir/TM_detect/blob/main/results/example_1prn.png)

An output file is also generated in the folder `./results/`. It gives information on the transmembrane segments. For example : 

``` 
Transmembrane segment from residue 10 to 10
Transmembrane segment from residue 12 to 12
Transmembrane segment from residue 27 to 27
Transmembrane segment from residue 31 to 31
Transmembrane segment from residue 33 to 33
Transmembrane segment from residue 51 to 51
Transmembrane segment from residue 53 to 53
Transmembrane segment from residue 55 to 55
Transmembrane segment from residue 61 to 66
Transmembrane segment from residue 68 to 69
...

```



## Structure
The source code is located in the `./src/` folder. 
This program is implemented using the OOP paradigm. 
The files [Protein.py](https://github.com/esoufir/TM_detect/blob/main/src/Protein.py) and [Geometry.py](https://github.com/esoufir/TM_detect/blob/main/src/Geometry.py) contain the classes. 
The file [TM_detect.py](https://github.com/esoufir/TM_detect/blob/main/src/TM_detect.py) contains the main instructions of the program. 
