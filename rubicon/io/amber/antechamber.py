# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
A wrapper for AntechamberRunner which generates force field files
or a specified molecule using gaussian output file as input
"""

import shlex
import subprocess
import tempfile
from collections import namedtuple

from monty.dev import requires
from monty.os.path import which
from monty.tempfile import ScratchDir

from pymatgen import Molecule
from pymatgen.io.lammps.data import Topology
from rubicon.io.lammps.topology import correct_corrupted_top_files, \
    TopCorruptionException
from pymatgen.io.lammps.data import ForceField
from rubicon.io.lammps.force_field import correct_corrupted_frcmod_files, \
    FFCorruptionException


__author__ = 'Navnidhi Rajput, Kiran Mathew'


class AntechamberRunner(object):
    """
    A wrapper for AntechamberRunner software
    """

    @requires((which('parmchk') or which('parmchk2')), "Requires the binary parmchk."
                                "Install AmberTools from http://ambermd.org/#AmberTools")
    @requires(which('antechamber'), "Requires the binary antechamber."
                                    "Install AmberTools from http://ambermd.org/#AmberTools")

    def __init__(self, mols):
        """
        Args:
            mols: List of molecules
        """
        self.mols = mols
        if which('parmchk'):
            self.parmchk_version = 'parmchk'
        else:
            self.parmchk_version = 'parmchk2'

    def _run_parmchk(self, filename="mol.mol2", format="mol2", outfile_name="mol.frcmod",
                     print_improper_dihedrals="Y"):
        """
        run parmchk
        """
        command = parmchk_version + " -i {} -f {} -o {} -w {}".format(filename, format, outfile_name,
                                                           print_improper_dihedrals)
        exit_code = subprocess.call(shlex.split(command))
        return exit_code

    def _run_antechamber(self, filename, infile_format="gout", outfile_name="mol",
                         outfile_format="mol2", charge_method="resp", status_info=2):
        """
        run antechamber using the provided gaussian output file
        """
        command = "antechamber -i {} -fi {} -o {}.{} -fo {} -c {} -s {}".format(filename,
                                                                             infile_format,
                                                                             outfile_name,
                                                                             outfile_format,
                                                                             outfile_format,
                                                                             charge_method,
                                                                             status_info)
        # dont think 'charmm' is even an option for -fo
        # GeneralizedForceFiled tries to read in *.ac(antechamber format) file and Toplogy
        # is trying to readin *.rtf(charmm format topology) file !!! WHY?!!
        # command = 'antechamber -i ' + filename + " -fi gout -o mol -fo charmm -c resp -s 2"
        exit_code = subprocess.call(shlex.split(command))
        return exit_code

    def _run_tleap(self, mol_name='mol', amb_ff_version='ff14SB'):
        '''
        run tleap
        '''
        lines = []
        lines.append('source leaprc.{}'.format(amb_ff_version))
        lines.append('source leaprc.gaff')
        lines.append('{} = loadmol2 {}.mol2'.format(mol_name,mol_name))
        lines.append('check {}'.format(mol_name))
        lines.append('loadamberparams {}.frcmod'.format(mol_name))
        lines.append('saveoff {} {}.lib'.format(mol_name,mol_name))
        lines.append('saveamberparm {} {}.prmtop {}.inpcrd'.format(mol_name,mol_name,mol_name))
        lines.append('quit')

        text = '\n'.join(lines)

        file = open('tleap.in', 'w')
        file.write(text)
        file.close()

        command = 'tleap -f tleap.in'

        exit_code = subprocess.call(shlex.split(command))
        return exit_code

    def _get_gaussian_ff_top_single(self, filename=None):
        """
        run antechamber using gaussian output file, then run parmchk
        to generate missing force field parameters. Store and return
        the force field and topology information in ff_mol.

        Args:
            filename: gaussian output file of the molecule

        Returns:
            Amberff namedtuple object that contains information on force field and
            topology
        """
        scratch = tempfile.gettempdir()
        Amberff = namedtuple("Amberff", ["force_field", "topology"])
        with ScratchDir(scratch, copy_from_current_on_enter=True,
                        copy_to_current_on_exit=True) as d:
            # self._convert_to_pdb(mol, 'mol.pdb')
            # self.molname = filename.split('.')[0]
            self._run_antechamber(filename)
            self._run_parmchk()
            # if antechamber can't find parameters go to gaff_example.dat
            try:
                mol = Molecule.from_file('mol.rtf')
                print('mol.rtf file exists')
            except TopCorruptionException:
                correct_corrupted_top_files('mol.rtf', 'gaff_data.txt')
                top = Topology.from_file('mol.rtf')
                print('mol.rtf file does not exist')
            try:
                gff = ForceField.from_file('mol.frcmod')
            except FFCorruptionException:
                correct_corrupted_frcmod_files('ANTECHAMBER.FRCMOD', 'gaff_data.txt')
                gff = ForceField.from_file('ANTECHAMBER.FRCMOD')
            # gff.set_atom_mappings('ANTECHAMBER_AC.AC')
            # gff.read_charges()
            # decorate the molecule with the sire property "atomname"
            #mol.add_site_property("atomname", (list(gff.atom_index.values())))
        return Amberff(gff, top)

    def get_gaussian_ff_top(self, filenames):
        """
        return a list of amber force field and topology for the list of
        gaussian output filenames corresponding to each molecule in mols list.

        Args:
            filenames (list): list of gaussian output files for each type of molecule

        Returns:
            list of Amberff namedtuples
        """
        amber_ffs = []
        for fname in filenames:
            amber_ffs.append(self._get_gaussian_ff_top_single(filename=fname))
        return amber_ffs
