from RDF import *
from pymatgen.core.structure import Structure
import time

time0=time.time()
struc = Structure.from_file("POSCAR-NaCl")
rdf1 = RDF(struc).RDF
print('Time Elasped (seconds):  ', time.time()-time0)

struc.make_supercell([2, 2, 2])
rdf2 = RDF(struc).RDF
print('Time Elasped (seconds):  ', time.time()-time0)
print('difference:           :  ', np.max(rdf1-rdf2))
