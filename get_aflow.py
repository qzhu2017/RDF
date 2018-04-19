from RDF import *
from aflow import *
from math import pi
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

results = search(batch_size=20
        ).filter(K.Egap > 1
        ).filter(K.nspecies == 2
        ).filter(K.natoms < 5)
print(len(results))

X = []
Y = []
for result in results:
    crystal = Structure.from_str(result.files['CONTCAR.relax.vasp'](), fmt='poscar')
    X.append(RDF(crystal).RDF[:,1])
    Y.append(result.Egap)

#Store this file
#Meachine Learning
