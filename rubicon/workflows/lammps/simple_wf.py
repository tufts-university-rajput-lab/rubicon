# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from fireworks import Workflow, LaunchPad
from fireworks.core.firework import Firework

from pymatgen.io.lammps.input import NVTLammpsInput

from rubicon.firetasks.lammps.input_tasks import WritelammpsInputFromDictInput
from rubicon.firetasks.lammps.run_tasks import RunLammpsDirect


__author__ = 'Kiran Mathew'


if __name__ == '__main__':
    data_file = "nvt.data"
    input_file = "nvt.inp"
    lammps_dict_input = NVTLammpsInput(data_file=data_file, is_forcefield=True)
    task1 = WritelammpsInputFromDictInput(lammps_dict_input=lammps_dict_input, input_file=input_file)
    task2 = RunLammpsDirect(lammps_cmd="pwd "+input_file)
    fw1 = Firework([task1], name='write input', fw_id=1)
    fw2 = Firework([task2], name='run lammps', fw_id=2)
    wf = Workflow([fw1, fw2], name="LAMMPS")
    lp = LaunchPad.auto_load()
    lp.add_wf(wf)
