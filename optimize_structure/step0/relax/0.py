import csv
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.core.bonds import CovalentBond
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

ratio=1.0

with open('../dataset/MoS2_96/structures/structures_1', 'r') as f:
    cont = f.read()
struc1 = Structure.from_str(cont, fmt='poscar')
fx,fy,fz = [],[],[]
with open('forces', 'r') as f:
    for line in f:
        fx.append(float(line.split()[-3]))
        fy.append(float(line.split()[-2]))
        fz.append(float(line.split()[-1]))
fx = np.array(fx)
fy = np.array(fy)
fz = np.array(fz)

m=struc1.lattice.matrix
lattice=Lattice(struc1.lattice.matrix)
num=struc1.num_sites
atom1, atom2=[],[]
for i in range(num):
    atom1.append(struc1.frac_coords[i])
    atom2.append(struc1.frac_coords[i])

f1 = open('new', 'w');
f2 = open('maxf', 'w');

f1.write("Mo67  S134 \n")
f1.write("1.0 \n")
f1.write("19.125407628  0.0000000000  0.0000000000 \n")
f1.write("-9.562703814  16.5630888636  0.0000000000 \n")
f1.write("0.0000000000  0.0000000000  39.3884239200 \n")
f1.write("Mo  S \n")
f1.write("67  134 \n")
f1.write("Direct\n")

mf=-10
for i in range(num):
    a1=lattice.get_cartesian_coords(atom1[i])
    atom2[i][0]=a1[0]+fx[i]*ratio
    atom2[i][1]=a1[1]+fy[i]*ratio
    atom2[i][2]=a1[2]+fz[i]*ratio
    temp=np.sqrt(fx[i]*fx[i])
    if(temp>mf):
        mf=temp
    temp=np.sqrt(fy[i]*fy[i])
    if(temp>mf):
        mf=temp
    temp=np.sqrt(fz[i]*fz[i])
    if(temp>mf):
        mf=temp
    a2=lattice.get_fractional_coords(atom2[i])
    f1.write("  %14.10f   %14.10f   %14.10f\n" % (a2[0],a2[1],a2[2]))
f2.write("  %14.10f \n" % (mf))
print(mf)
