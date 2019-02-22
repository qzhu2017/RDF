# RDF
A code to compute the radial distribution function
#edit
This branch is initally planned for the undergraduate student's project on machine learning materials properties.

Currently, there are four classes,
- Element
- crystal
- cif
- RDF

One could load the crystal from 
- dictionary
- POSCAR
- CIF file 

To perform RDF calculation, one needs to provide the following info
- crystal structure
- R_max (default 12 \AA)
- step length (defult: 0.08 \AA)
- sigma width for Gaussian smearing (default: 0.10 \AA)


## Requirement
This code need the installation of Python3, scipy, numpy, pymatgen. The easiest way is to install them via anaconda or miniconda.

$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
$ conda install -c matsci pymatgen


## usage
```
$ pythonw RDF.py -h
Usage: RDF.py [options]

Options:
  -h, --help            show this help message and exit
  -c STRUCTURE, --crystal=STRUCTURE
                        crystal from file, cif or poscar, REQUIRED
  -r Rmax, --Rmax=Rmax  Rmax, default: 12 A
  -s SIGMA, --sigma=SIGMA
                        sigma width for Gaussian smearing, default: 0.20
  -d DELTA, --delta=DELTA
                        step length, default: 0.08
  -o MSTYLE, --output=MSTYLE
                        matplotlib style, fivethirtyeight, bmh, grayscale,
                        dark_background, ggplot
  -p PLOT, --plot=PLOT  generate the plot to file, default: None
 ```
 ## execution 
 one just needs to run the followings,
```
$ python RDF.py -c POSCAR-Ar
```
It will generate a png file with RDF plot as follows
![FCC-Ar](https://github.com/qzhu2017/RDF/blob/master/images/Ar.png)
