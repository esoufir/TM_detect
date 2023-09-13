SOUFIR Emma M2BI - 2023/2024  
Projet court - Assignation et détection des parties transmembranaires d'une protéine

# TM_detect

TM_detect is a program which locates the membrane in a transmbrane protein and detects transmebrane segments. 

## Installation

### Dependencies


Download the ZIP repository.   
In a proper folder, unzip it.  

Change your current directory to the root of the project. The `ls -l` command should give you : 

```
-rw-r--r-- 1 emmas emmas   2800 Sep 13 11:15  README.md
-rw-r--r-- 1 emmas emmas 445283 Sep 13  2023  SOUFIR_Emma_Rapport_projet_court.pdf
drwxr-xr-x 2 emmas emmas   4096 Sep 13 09:17  data/
drwxr-xr-x 2 emmas emmas   4096 Sep  4 17:23  doc/
-rw-r--r-- 1 emmas emmas   2675 Sep 13 10:50  environment.yml
drwxr-xr-x 2 emmas emmas   4096 Sep 13 10:13  results/
drwxr-xr-x 3 emmas emmas   4096 Sep 13 10:47  src/
```

Then, load and activate the virtual environment [environment.yml](https://github.com/esoufir/TM_detect/blob/main/environment.yml) :
```
conda env create --file environment.yml
conda activate mypymolenv
```



### Virtual environment

You can use the virtual environment [environment.yml](https://github.com/esoufir/TM_detect/blob/main/environment.yml). 
This program uses the tool DSSP and PyMol. Everything is detailed in the virtual environment [environment.yml](https://github.com/esoufir/TM_detect/blob/main/environment.yml).
Another environment using pymol-open-source is available [environment-open-source.yml](https://github.com/esoufir/TM_detect/blob/main/environment-open-source.yml). 

```
conda env create --file environment_open_source.yml
conda activate mypymolenv
```

## Usage

After activating the virtual environment, run the following command : 
```
python src/TM_detect.py [-h] [-n N] [-w WIDTH] [-g GAP] [-m GAP_MEMBRANE] path/to/filename

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
python src/TM_detect.py data/1prn.pdb -n 15 -w 14.0 -g 1.0 -m 1.0
```
opens PyMol GUI and gives the output : 
![Example](https://github.com/esoufir/TM_detect/blob/main/results/example_1prn.png)

An output file is also generated at the root of the project. It gives the coordinates of the planes representing the membrane. For example : 

``` 
CA	 -24.43600082397461	 14.550000190734863	 19.768967474151072
CA	 -24.43600082397461	 17.550000190734863	 20.644985246959646
CA	 -24.43600082397461	 20.550000190734863	 21.52100301976822
CA	 -24.43600082397461	 23.550000190734863	 22.397020792576793
CA	 -24.43600082397461	 26.550000190734863	 23.273038565385363
CA	 -24.43600082397461	 29.550000190734863	 24.149056338193937
CA	 -24.43600082397461	 32.55000019073486	 25.02507411100251
...
```

You can find some examples of outputs in the `results/` folder, where the script was used on the `data/` files with default arguments. 

## Structure

The source code is located in the `./src/` folder. 
This program is implemented using the OOP paradigm. 
The files [Protein.py](https://github.com/esoufir/TM_detect/blob/main/src/Protein.py) and [Geometry.py](https://github.com/esoufir/TM_detect/blob/main/src/Geometry.py) contain the classes. 
The file [TM_detect.py](https://github.com/esoufir/TM_detect/blob/main/src/TM_detect.py) contains the main instructions of the program. 
