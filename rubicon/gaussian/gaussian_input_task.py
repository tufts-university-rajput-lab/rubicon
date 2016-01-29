import shlex
import subprocess


__author__ = 'navnidhirajput'

from fireworks import FireTaskBase, explicit_serialize, FWAction
from pymatgen import write_mol, Molecule
import glob
from pymatgen.io import gaussian
import os


@explicit_serialize
class WritegaussianGeoTask(FireTaskBase):
    """
    Writes Gaussian Input files for geometry optimization.

    """

    _fw_name = "Gaussian Geometry Writer"

    def run_task(self, fw_spec):
        print "dir",os.getcwd()
        mol = fw_spec["molecule"]
        mol_name = fw_spec["mol_name"]
        charge = fw_spec["charge"]
        spin_multiplicity = fw_spec["spin_multiplicity"]

        gaus_lines = gaussian.GaussianInput(mol, charge = charge,
                                            spin_multiplicity = spin_multiplicity,
                                            title='created by gaussian_geo_task from' + ' ' + mol_name,
                                            functional="b3lyp",
                                            basis_set="aug-cc-pvdz",
                                            route_parameters={
                                                'opt': "(calcfc,tight)",
                                                'int': "ultrafine",
                                                "\n# SCF": "tight  nosymm test"},
                                            input_parameters=None,
                                            link0_parameters={
                                                "%mem": "256MW",
                                                "%NProcShared": 4,
                                                "%LindaWorker": "localhost",
                                                "%chk": mol_name + ".chk"},
                                            dieze_tag="#",
                                            gen_basis=None)


        gaus_lines.write_file('mol_geo.gau', cart_coords=True)


        with open('mol_geo.gau') as f, open("mol_geo.out", 'w') as fo :
           subprocess.call(shlex.split("g09launch"), stdin=f, stdout = fo)

        prev_gaussian_geo = shlex.os.path.join(shlex.os.getcwd(), mol_name+'.out')
        update_spec = {'prev_gaussian_geo': prev_gaussian_geo}

        return FWAction(update_spec=update_spec)



@explicit_serialize
class WritegaussianFreqTask(FireTaskBase):
    """
    Writes Gaussian Input files for frequency and charge calculation.

    """

    _fw_name = "Gaussian frequency and charge Writer"

    def run_task(self, mol, charge=None, spin_multiplicity=None):
        moleculelist = glob.glob(
            "/Users/navnidhirajput/Dropbox/solvent_molecules/*")
        for filename in moleculelist:
            mol = Molecule.from_file(filename)

            gaus_freq_charge = gaussian.GaussianInput(mol, charge = charge,
                                                      spin_multiplicity=spin_multiplicity,
                                                      title='created by gaussian_frq_task from' + ' ' + filename[
                                                                                                          48:],
                                                      functional="b3lyp",
                                                      basis_set="aug-cc-pvdz  freq",
                                                      route_parameters={
                                                          '\n# geom': "allcheck",
                                                          'guess': "read",
                                                          "SCF": "tight  nosymm test",
                                                      "pop":"MK iop(6/33=2,6/41=10,6/42=10,7/33=1)"},
                                                      input_parameters=None,
                                                      link0_parameters={
                                                          "%mem": "256MW",
                                                          "%NProcShared": 4,
                                                          "%LindaWorker": "localhost",
                                                          "%chk": filename[
                                                                  48:-4] + ".chk"},
                                                      dieze_tag="#",
                                                      gen_basis=None)
        with open('mol_freq.gau', 'w') as f:
            gaus_freq_charge.write_file('mol_freq.gau', cart_coords=True)


        subprocess.check_call(shlex.split("g09 "), stdin=f)

        return FWAction()

if __name__ == '__main__':
    task_geo = WritegaussianGeoTask()
    moleculelist = glob.glob("/Users/navnidhirajput/Dropbox/solvent_molecules/*")
    for filename in moleculelist:
        mol = Molecule.from_file(filename)
        file_name = os.path.basename(filename)
        task_geo.run_task({"molecule":mol, "mol_name": os.path.splitext(file_name)[0], "charge": 1,"spin_multiplicity":-1})