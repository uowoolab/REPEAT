#!/usr/bin/env python

"""
Generate a REPEAT compatible gaussian .cube file from a VASP
electrostatic potential (LOCPOT and, for older versions, the OUTCAR
too). Just run the script with no arguments in the working directory.

"""

__version__ = "2.0"

import bz2
import gzip
import re
from glob import glob


def compressed_open(filename):
    """Return file objects for either compressed and uncompressed files"""
    filename = glob(filename) + glob(filename+".gz") + glob(filename+".bz2")
    try:
        filename = filename[0]
    except IndexError:
        print("File %s not found" % filename)
        return None
    if filename[-4:] == ".bz2":
        return bz2.BZ2File(filename)
    elif filename[-3:] == ".gz":
        return gzip.open(filename)
    else:
        return open(filename, "r")


def outcar_types(ntypes=None, file_name='OUTCAR'):
    """Find the species present in the simulation from the OUTCAR."""
    otypes = []
    outcar = compressed_open(file_name)
    if ntypes is None:
        for oline in outcar:
            if 'ions per type' in oline:
                ntypes = len(oline.split()) - 4
                outcar.seek(0)
                break
    for oline in outcar:
        if 'INCAR:' in oline:
            for _type_idx in range(ntypes):
                otype = outcar.next().split()[2].split('_')[0]
                otypes.append(otype)
    return otypes

# Find the data files
in_file_name = 'LOCPOT'
locpot = compressed_open(in_file_name)

# Top line identifies system -- compress spaces
name = re.sub('\\s+', ' ', locpot.readline().strip())
out_name = "%s.cube" % name
cube_out = open(out_name, 'w')

# Cube is in Bohr, so we can alter the scale here to account for it
scale = float(locpot.readline())
ANG_TO_BOHR = 1.889725985
scale = scale*ANG_TO_BOHR

lattice = [[scale*float(i) for i in locpot.readline().split()]
           for _ in ['a', 'b', 'c']]

types = locpot.readline().split()
if types[0].isalpha():
    atom_counts = [int(x) for x in locpot.readline().split()]
else:
    # Old vasp does not have atom types in CONTCAR-like output
    atom_counts = [int(x) for x in types]
    types = outcar_types(ntypes=len(atom_counts))

if locpot.readline().strip()[0].lower() in 'ck':  # Cartesian coords
    mcell = [[scale, 0.0, 0.0], [0.0, scale, 0.0], [0.0, 0.0, scale]]
else:  # Direct (fractional)
    mcell = lattice

atoms = []
for _atom_idx in range(sum(atom_counts)):
    coord = [float(x) for x in locpot.readline().split()]
    coord = [mcell[0][0]*coord[0] + mcell[1][0]*coord[1] + mcell[2][0]*coord[2],
             mcell[0][1]*coord[0] + mcell[1][1]*coord[1] + mcell[2][1]*coord[2],
             mcell[0][2]*coord[0] + mcell[1][2]*coord[1] + mcell[2][2]*coord[2]]
    atoms.append(coord)

locpot.readline()  # Blank line
ngrid = [int(x) for x in locpot.readline().split()]

# Just read everything in to one long list
factor_e = -1/27.212  # LOCPOT in eV
v = [factor_e*float(x) for line in locpot for x in line.split()]
locpot.close()

# Extra data for cube file

origin = (0, 0, 0)
cube_lattice = [[(x/ngrid[0]) for x in lattice[0]],
                [(y/ngrid[1]) for y in lattice[1]],
                [(z/ngrid[2]) for z in lattice[2]]]

ATOMIC_NUMBER = [
    "ZERO", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Uut", "Uuq", "Uup", "Uuh", "Uuo"]

# Start writing here

header = ['--------------------REPEAT charges--------------------\n',
          '---cube file created from VASP (LOCPOT and OUTCAR)----\n',
          '% 5s' % sum(atom_counts),
          '% 12.6f% 12.6f% 12.6f\n' % origin,
          '% 5s' % ngrid[0],
          '% 12.6f% 12.6f% 12.6f\n' % tuple(cube_lattice[0]),
          '% 5s' % ngrid[1],
          '% 12.6f% 12.6f% 12.6f\n' % tuple(cube_lattice[1]),
          '% 5s' % ngrid[2],
          '% 12.6f% 12.6f% 12.6f\n' % tuple(cube_lattice[2])]

cube_out.writelines(header)

for atype, count in zip(types, atom_counts):
    for _idx in range(count):
        atom = atoms.pop(0)
        cube_out.write('%5.0f% 12.6f% 12.6f% 12.6f% 12.6f\n' %
                       (ATOMIC_NUMBER.index(atype), ATOMIC_NUMBER.index(atype),
                        atom[0], atom[1], atom[2]))

# ESP needs newlines every 6 or each z loop
count = 0
for idx in range(ngrid[0]):
    for idy in range(ngrid[1]):
        for idz in range(ngrid[2]):
            cube_out.write('%13.5E' % v[idx+idy*ngrid[0]+idz*ngrid[0]*ngrid[1]])
            count += 1
            if not count % ngrid[2]:
                cube_out.write('\n')
                count = 0
            elif not count % 6:
                cube_out.write('\n')

cube_out.close()

print("Returned file %s" % out_name)
