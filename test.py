from RDF import *
from aflow import *
from math import pi
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

def convert(struc):
    """
    convert pymatgen struct class to our own struc class
    """

    composition=[]
    atom_type = []
    dict = struc.composition.as_dict()

    for ele in struc.composition.elements:
        composition.append(int(dict[ele.name]))
        atom_type.append(ele.name)

    crystal0 = crystal(fileformat='poscar',
        lattice= struc.lattice.matrix,
        atom_type = atom_type,
        composition = composition,
        coordinate = struc.frac_coords)
    return crystal0


struc = Structure.from_file("POSCAR-NaCl")
rdf = RDF(convert(struc)).RDF
print(rdf)

struc.make_supercell([2, 2, 2])
rdf = RDF(convert(struc)).RDF
print(rdf)
