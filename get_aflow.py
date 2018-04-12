from RDF import *
from aflow import *
from math import pi
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

results = search(batch_size=20
        ).filter(K.Egap > 1
        ).filter(K.nspecies == 2
        ).filter(K.natoms == 2)

for result in results:
    para = result.geometry
    para[3:] = para[3:]*pi/180
    lattice = crystal.para2matrix(para)
    atom_type = result.species
    composition = result.composition
    coordinate = result.positions_fractional
    site = []
    
    for numIon, ele in zip(composition, atom_type):
        ele = ele.replace('\n','')
        for x in range(numIon):
            site.append(ele)
   
    struc = Structure(lattice, site, coordinate)
    finder = SpacegroupAnalyzer(struc, symprec=0.06, angle_tolerance=5)
    struc = finder.get_conventional_standard_structure()
    composition=[]
    dict = struc.composition.as_dict()
    for ele in struc.composition.elements:
        composition.append(int(dict[ele.name]))

    crystal0 = crystal(fileformat='poscar', 
        lattice= struc.lattice.matrix,
        atom_type = atom_type,
        composition = composition,
        coordinate = struc.frac_coords)
    rdf = RDF(crystal0).RDF
    print(rdf)  # this is X
    print(result.Egap) #this is Y
