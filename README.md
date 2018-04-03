# RDF
A code to compute the radial distribution function

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


More detailed usage could be found in the [jupyter notebook](to add in future).

## usage
```
$ python RDF.py -h
Usage: RDF.py [options]

Options:
  -h, --help            show this help message and exit
  -c crystal, --crystal=crystal
                        crystal from file, cif or poscar, REQUIRED
  -r Rmax, --Rmax=Rmax  Rmax, default: 12 A
  -s sigma, --sigma=sigma
                        sigma width for Gaussian smearing, default: 0.20
  -d R_bin, --delta=R_bin
                        step length, default: 0.08
  -o mstyle, --output=mstyle
                        matplotlib style, fivethirtyeight, bmh, grayscale,
 ```
 ## execute 
 one just needs to run the followings,
```
$ python RDF.py -c POSCAR-Ar
[[  8.00000000e-02   0.00000000e+00]
 [  1.60000000e-01   0.00000000e+00]
 [  2.40000000e-01   0.00000000e+00]
 [  3.20000000e-01   0.00000000e+00]
 [  4.00000000e-01   0.00000000e+00]
 [  4.80000000e-01   0.00000000e+00]
 [  5.60000000e-01   0.00000000e+00]
 [  6.40000000e-01   0.00000000e+00]
 [  7.20000000e-01   0.00000000e+00]
 [  8.00000000e-01   0.00000000e+00]
 [  8.80000000e-01   0.00000000e+00]
 [  9.60000000e-01   0.00000000e+00]
 [  1.04000000e+00   0.00000000e+00]
 [  1.12000000e+00   0.00000000e+00]
 [  1.20000000e+00   0.00000000e+00]
 [  1.28000000e+00   0.00000000e+00]
 [  1.36000000e+00   0.00000000e+00]
 [  1.44000000e+00   0.00000000e+00]
 [  1.52000000e+00   0.00000000e+00]
 [  1.60000000e+00   0.00000000e+00]
 [  1.68000000e+00   0.00000000e+00]
 [  1.76000000e+00   0.00000000e+00]
 [  1.84000000e+00   0.00000000e+00]
 [  1.92000000e+00   0.00000000e+00]
 [  2.00000000e+00   0.00000000e+00]
 [  2.08000000e+00   0.00000000e+00]
 [  2.16000000e+00   0.00000000e+00]
 [  2.24000000e+00   0.00000000e+00]
 [  2.32000000e+00   0.00000000e+00]
 [  2.40000000e+00   0.00000000e+00]
 [  2.48000000e+00   0.00000000e+00]
 [  2.56000000e+00   0.00000000e+00]
 [  2.64000000e+00   0.00000000e+00]
 [  2.72000000e+00   0.00000000e+00]
 [  2.80000000e+00   0.00000000e+00]
 [  2.88000000e+00   0.00000000e+00]
 [  2.96000000e+00   0.00000000e+00]
 [  3.04000000e+00   0.00000000e+00]
 [  3.12000000e+00   0.00000000e+00]
 [  3.20000000e+00   0.00000000e+00]
 [  3.28000000e+00   0.00000000e+00]
 [  3.36000000e+00   0.00000000e+00]
 [  3.44000000e+00   0.00000000e+00]
 [  3.52000000e+00   0.00000000e+00]
 [  3.60000000e+00   0.00000000e+00]
...
...
```
It will also generate a png file with RDF plot as follows
![FCC-Ar](https://github.com/qzhu2017/RDF/blob/master/images/Ar.png)
