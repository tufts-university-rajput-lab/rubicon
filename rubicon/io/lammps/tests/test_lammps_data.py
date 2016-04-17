# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import unittest
import os

import numpy as np

from pymatgen.core.structure import Molecule

from rubicon.io.lammps.data import LammpsData


__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestLammpsData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        h2o_coords = [[9.626, 6.787, 12.673],
                      [9.626, 8.420, 12.673],
                      [10.203, 7.604, 12.673]]
        molecule = Molecule(["H", "H", "O"], h2o_coords)
        box_size = [[0.0, 10.0],[0.0,10.0],[0.0, 10.0]]
        cls.lammps_data = LammpsData.from_structure(molecule, box_size)

    def test_system_info(self):
        atomic_masses = [[1, 1.00794], [2, 15.9994]]
        atoms_data = [[1, 1, 1, 0.0, 4.4875653445297559, 4.1830559491720365, 5.0000000000000018],
                      [2, 1, 1, 0.0, 4.4875653445297559, 5.8160559491720365, 5.0000000000000018],
                      [3, 1, 2, 0.0, 5.0645653445297558, 5.0000559491720367, 5.0000000000000018]]
        natom_types = 2
        natoms = 3
        np.testing.assert_almost_equal(self.lammps_data.atomic_masses, atomic_masses, decimal=10)
        np.testing.assert_almost_equal(self.lammps_data.atoms_data, atoms_data, decimal=10)
        self.assertEqual(self.lammps_data.natom_types, natom_types)
        self.assertEqual(self.lammps_data.natoms, natoms)

    def test_string_representation(self):
        string_rep = 'Data file generated by rubicon\n\n' \
                     '3 atoms\n\n' \
                     '2 atom types\n\n' \
                     '0.0 10.0 xlo xhi\n' \
                     '0.0 10.0 ylo yhi\n' \
                     '0.0 10.0 zlo zhi\n\n' \
                     'Masses \n\n' \
                     '1 1.00794\n' \
                     '2 15.9994\n\n' \
                     'Atoms \n\n' \
                     '1 1 1 0.0 4.48756534453 4.18305594917 5.0\n' \
                     '2 1 1 0.0 4.48756534453 5.81605594917 5.0\n' \
                     '3 1 2 0.0 5.06456534453 5.00005594917 5.0'
        self.assertEqual(str(self.lammps_data), string_rep)

    def test_from_file(self):
        self.lammps_data.write_data_file(os.path.join(module_dir, "lammps_data.dat"))
        lammps_data = LammpsData.from_file(os.path.join(module_dir, "lammps_data.dat"))
        self.assertEqual(str(lammps_data), str(self.lammps_data))

    def tearDown(self):
        for x in ["lammps_data.dat"]:
            if os.path.exists(os.path.join(module_dir, x)):
                os.remove(os.path.join(module_dir, x))


if __name__ == "__main__":
    unittest.main()
