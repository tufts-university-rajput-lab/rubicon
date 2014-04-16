from unittest import TestCase
from pymatgen import Molecule
from rubicon.gff.antechamberio import AntechamberRunner
from rubicon.gff.gff import Gff
from rubicon.gff.topology import TopMol

__author__ = 'navnidhirajput'

coords = [[0.000000, 0.000000, 0.000000],
          [0.000000, 0.000000, 1.089000],
          [1.026719, 0.000000, -0.363000],
          [-0.513360, -0.889165, -0.363000],
          [-0.513360, 0.889165, -0.363000]]

mol = Molecule(["C", "H", "H", "H", "H"], coords)


class TestAntechamber(TestCase):



    def test_get_ff_bonds(self):
        my_ant = AntechamberRunner(mol)
        my_gff = Gff()
        top = TopMol.from_file('mol.rtf')

        my_ant._get_ff_bonds(my_gff.bonds, top.bonds, my_ant.atom_gaff)
        ans_bond={'C-H1': ('c3-hc', (337.3, 1.092)), 'C-H2': ('c3-hc', (337.3, 1.092)), 'C-H3': ('c3-hc', (337.3, 1.092)), 'C-H': ('c3-hc', (337.3, 1.092))}
        self.assertEquals(ans_bond,my_ant.topbondff)

    def test_get_ff_angles(self):
        my_ant = AntechamberRunner(mol)
        my_gff = Gff()
        top = TopMol.from_file('mol.rtf')

        my_ant._get_ff_angles(my_gff.angles, top.angles, my_ant.atom_gaff)
        ans_angle={'H1-C-H3': ('hc-c3-hc', (39.43, 108.35)), 'H1-C-H2': ('hc-c3-hc', (39.43, 108.35)), 'H2-C-H3': ('hc-c3-hc', (39.43, 108.35)), 'H-C-H3': ('hc-c3-hc', (39.43, 108.35)), 'H-C-H2': ('hc-c3-hc', (39.43, 108.35)), 'H-C-H1': ('hc-c3-hc', (39.43, 108.35))}
        self.assertEquals(ans_angle,my_ant.topangleff)

    def test_get_ff_dihedrals(self):
        my_ant = AntechamberRunner(mol)
        my_gff = Gff()
        top = TopMol.from_file('mol.rtf')

        my_ant._get_ff_dihedrals(my_gff.dihedrals, top.dihedrals, my_ant.atom_gaff)
        ans_dihedral={}
        self.assertEquals(ans_dihedral,my_ant.topdihedralff)

    def test_get_ff_imdihedrals(self):
        my_ant = AntechamberRunner(mol)
        my_gff = Gff()
        top = TopMol.from_file('mol.rtf')

        my_ant._get_ff_imdihedrals(my_gff.imdihedrals, top.imdihedrals, my_ant.atom_gaff)
        ans_imdihedral={}
        self.assertEquals(ans_imdihedral,my_ant.topimdihedralff)

    def test_run_antechamber(self):
        ant=AntechamberRunner(mol)
        return_cmd,my_gff,ant,my_lammps_list,gff_list=ant._run_antechamber('mol.pdb',mol)
        output= '\n'.join(my_lammps_list.lines)
        ans='''LAMMPS Data File


2 atom type
1 bond type
1 angle type
0 dihedral type
0 improper type

15 atoms
12 bonds
18 angles
0 dihedrals
0 impropers

0.0 40.0 xlo  xhi
0.0 40.0 ylo  yhi
0.0 40.0 zlo  zhi

Masses

1 12.01 # c3
2 1.008 # hc


Pair Coeffs

1 3.3996744 0.1094 # c3
2 2.6495366 0.0157 # hc


Bond Coeffs

1 337.3 1.092 # c3 hc


Angle Coeffs

1 39.43 108.35 # c3 hc hc


Dihedral Coeffs



Imp Dihedral Coeffs




Atoms

1  1  1  1.105  39.124  38.656 # 1  c3 C
2  1  2  1.113  39.506  37.636 # 1  hc H
3  1  2  0.188  38.561  38.825 # 1  hc H1
4  1  2  1.964  38.472  38.805 # 1  hc H2
5  1  2  1.155  39.956  39.356 # 1  hc H3
6  2  1  0.427  3.75  3.81 # 2  c3 C
7  2  2  0.004  4.753  3.806 # 2  hc H
8  2  2  1.508  3.811  3.686 # 2  hc H1
9  2  2  0.198  3.263  4.757 # 2  hc H2
10  2  2  0.0  3.172  2.993 # 2  hc H3
11  3  1  1.775  0.712  3.992 # 3  c3 C
12  3  2  1.329  1.376  4.732 # 3  hc H
13  3  2  2.498  0.059  4.479 # 3  hc H1
14  3  2  0.995  0.11  3.53 # 3  hc H2
15  3  2  2.277  1.305  3.229 # 3  hc H3

Bonds

1  1  1  2  #  1  C  H
2  1  1  3  #  1  C  H1
3  1  1  4  #  1  C  H2
4  1  1  5  #  1  C  H3
5  1  1  2  #  2  C  H
6  1  1  3  #  2  C  H1
7  1  1  4  #  2  C  H2
8  1  1  5  #  2  C  H3
9  1  1  2  #  3  C  H
10  1  1  3  #  3  C  H1
11  1  1  4  #  3  C  H2
12  1  1  5  #  3  C  H3

Angles

1  1  2  1  3  #  1  H  C  H1
2  1  2  1  4  #  1  H  C  H2
3  1  2  1  5  #  1  H  C  H3
4  1  3  1  4  #  1  H1  C  H2
5  1  3  1  5  #  1  H1  C  H3
6  1  4  1  5  #  1  H2  C  H3
7  1  2  1  3  #  2  H  C  H1
8  1  2  1  4  #  2  H  C  H2
9  1  2  1  5  #  2  H  C  H3
10  1  3  1  4  #  2  H1  C  H2
11  1  3  1  5  #  2  H1  C  H3
12  1  4  1  5  #  2  H2  C  H3
13  1  2  1  3  #  3  H  C  H1
14  1  2  1  4  #  3  H  C  H2
15  1  2  1  5  #  3  H  C  H3
16  1  3  1  4  #  3  H1  C  H2
17  1  3  1  5  #  3  H1  C  H3
18  1  4  1  5  #  3  H2  C  H3

Dihedrals


Impropers
'''

        self.assertEquals(ans,output)
