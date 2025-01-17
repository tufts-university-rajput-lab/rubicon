# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import copy
import datetime
import json
import logging
import os
import re
import shlex
import shutil
import socket
import sys
import random
import string

from custodian.custodian import Custodian
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QchemJob
from rubicon.workflows.qchem.wf_settings import MOVE_TO_EG_GARDEN

from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.fw_config import FWData
from fireworks.utilities.fw_serializers import FWSerializable
from pymatgen.core.structure import Molecule
from pymatgen.io.qchem import QcInput, QcTask
from rubicon.utils.eg_wf_utils import move_to_eg_garden

__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Jun 07, 2013'

DATETIME_HANDLER = lambda obj: obj.isoformat() \
    if isinstance(obj, datetime.datetime) else None


class QChemTask(FireTaskBase, FWSerializable):
    """
    Write QChem task and run QChem
    """

    _fw_name = "QChem Task"

    def run_task(self, fw_spec):
        qcinp = self._get_qcinp_from_fw_spec(fw_spec)
        mixed_basis = fw_spec.get("mixed_basis", None)
        mixed_aux_basis = fw_spec.get("mixed_aux_basis", None)
        implicit_solvent = fw_spec.get("implicit_solvent", None)

        custodian_out, prev_qchem_dir = self.run_qchem(
            qcinp, implicit_solvent, mixed_aux_basis, mixed_basis)
        all_errors = set()
        for run in custodian_out:
            for correction in run['corrections']:
                all_errors.update(correction['errors'])

        if MOVE_TO_EG_GARDEN:
            prev_qchem_dir = move_to_eg_garden(prev_qchem_dir)

        stored_data = {'error_list': list(all_errors)}
        update_spec = {'prev_qchem_dir': prev_qchem_dir,
                       'prev_task_type': fw_spec['task_type']}
        propagate_keys = ['egsnl', 'snlgroup_id', 'inchi_root',
                          'mixed_basis', 'mixed_aux_basis', 'mol']

        for k in propagate_keys:
            if k in fw_spec:
                update_spec[k] = fw_spec[k]

        return FWAction(stored_data=stored_data, update_spec=update_spec)

    @classmethod
    def set_nodelist_env(cls, ):
        fw_data = FWData()
        if fw_data.MULTIPROCESSING and fw_data.NODE_LIST is not None:
            nodelist = ",".join(fw_data.NODE_LIST)
            os.environ["QCNODE"] = nodelist

    @classmethod
    def get_hostname(cls):
        if "NERSC_HOST" in os.environ:
            hostname = os.environ["NERSC_HOST"]
        elif 'vesta' in socket.gethostname():
            hostname = "vesta"
        else:
            hostname = "unknow"
        return hostname

    @classmethod
    def get_qchem_cmd(cls, qcinp, mol):
        physical_nproc_map = {"edison": 24,
                              "cori": 32,
                              "matgen": 16,
                              "macqu": 2}
        numa_num_map = {"edison": 2,
                        "cori": 2,
                        "matgen": 2,
                        "macqu": 1}
        hostname = cls.get_hostname()
        if hostname in physical_nproc_map:
            physical_nproc = physical_nproc_map[hostname]
        elif hostname == 'vesta':
            return ALCF_Utils.get_alcf_qchem_cmd(qcinp)
        else:
            return ["qchem"] * 3
        natoms = len(mol)
        fw_data = FWData()
        if fw_data.MULTIPROCESSING and fw_data.SUB_NPROCS is not None:
            physical_nproc = int(fw_data.SUB_NPROCS)
        nproc = min(physical_nproc, natoms)
        numa_num = numa_num_map[hostname]
        nproc = max((nproc // numa_num) * numa_num, 1)
        half_nproc = min(physical_nproc // 2, natoms)
        half_nproc = max((half_nproc // numa_num) * numa_num, 1)
        qc_exe = shlex.split("qchem -np {}".format(nproc))
        half_cpus_cmd = shlex.split("qchem -np {}".format(half_nproc))
        openmp_cmd = shlex.split("qchem -nt {}".format(nproc))
        return qc_exe, half_cpus_cmd, openmp_cmd

    @classmethod
    def get_physical_memory(cls):
        physical_nproc_map = {"edison": 64,
                              "cori": 128,
                              "matgen": 64,
                              "macqu": 8,
                              "vesta": 18.5}
        if "NERSC_HOST" in os.environ:
            hostname = os.environ["NERSC_HOST"]
        elif 'vesta' in socket.gethostname():
            hostname = "vesta"
        else:
            hostname = "unknow"
        return physical_nproc_map.get(hostname, 104) - 4

    @classmethod
    def clean_up(cls, qcinp):
        ALCF_Utils.alcf_clean_up(qcinp)

    @classmethod
    def run_qchem(cls, qcinp, implicit_solvent, mixed_aux_basis, mixed_basis,
                  input_file="mol.qcinp", output_file="mol.qcout", gzipped=True):
        mol = qcinp.jobs[0].mol
        num_atoms = len(mol)
        for qj in qcinp.jobs:
            if qj.params["rem"]["jobtype"] != "sp":
                if mixed_basis is not None:
                    qj.set_basis_set(mixed_basis)
                if mixed_aux_basis is not None:
                    qj.set_aux_basis_set(mixed_aux_basis)
        prev_qchem_dir = os.getcwd()
        qc_exe, half_cpus_cmd, openmp_cmd = cls.get_qchem_cmd(qcinp, mol)
        logging.basicConfig(level=logging.INFO)
        qchem_logger = logging.getLogger('QChemDrone')
        qchem_logger.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setLevel(getattr(logging, 'INFO'))
        qchem_logger.addHandler(sh)
        scf_max_cycles = 200
        geom_max_cycles = 200
        alt_cmd = {"half_cpus": half_cpus_cmd,
                   "openmp": openmp_cmd}
        if cls._is_openmp_only_job(qcinp):
            qc_exe = openmp_cmd
            alt_cmd["half_cpus"] = shlex.split(" ".join(half_cpus_cmd).replace("-np", "-nt"))
            alt_cmd.pop("openmp")
        if num_atoms > 50:
            scf_max_cycles = 300
            geom_max_cycles = 500
        qcinp.write_file(input_file)
        if implicit_solvent is not None:
            solvent_data = implicit_solvent.get('solvent_data', None)
            if solvent_data is not None:
                values = ['{:.4f}'.format(solvent_data[t]) for t in
                          ['Dielec', 'SolN', 'SolA', 'SolB', 'SolG', 'SolC',
                           'SolH']]
                solvent_text = ' '.join(values)
                with open('solvent_data', 'w') as f:
                    f.write(solvent_text)
        qclog_file = os.path.splitext(output_file)[0] + ".qclog"
        total_physical_memory = cls.get_physical_memory()
        job = QchemJob(qc_exe, input_file=input_file, output_file=output_file,
                       qclog_file=qclog_file, alt_cmd=alt_cmd, gzipped=gzipped,
                       total_physical_memory=total_physical_memory)
        handler = QChemErrorHandler(qchem_job=job, input_file=input_file,
                                    output_file=output_file,
                                    scf_max_cycles=scf_max_cycles,
                                    geom_max_cycles=geom_max_cycles)
        c = Custodian(handlers=[handler], jobs=[job], max_errors=50)
        custodian_out = c.run()
        cls.clean_up(qcinp)
        return custodian_out, prev_qchem_dir

    @staticmethod
    def _is_openmp_only_job(qcinp):
        for qctask in qcinp.jobs:
            rems = qctask.params["rem"]
            theor_method = rems.get("method", rems.get("exchange", "hf"))
            for mp2task in ["mp2", "xyg", "ccsd"]:
                if mp2task in theor_method:
                    return True
        return False

    @staticmethod
    def _get_qcinp_from_fw_spec(fw_spec):
        if isinstance(fw_spec["qcinp"], dict):
            qcinp = QcInput.from_dict(fw_spec["qcinp"])
        else:
            qcinp = fw_spec["qcinp"]
        if 'mol' in fw_spec:
            if isinstance(fw_spec["mol"], dict):
                mol = Molecule.from_dict(fw_spec["mol"])
            else:
                mol = fw_spec["mol"]
            for qj in qcinp.jobs:
                if isinstance(qj.mol, Molecule):
                    qj.mol = copy.deepcopy(mol)
        return qcinp



class ALCF_Utils(object):

    @staticmethod
    def _calibrate_alcf_cmd(input_file="mol.qcinp", max_minutes=60,
                            num_nodes=8, ranks_per_node=1,
                            num_threads=64, scr_size_GB=2.0, use_ramdisk=True,
                            fs_scr_dir="", use_runjob=False):
        qc_package_path = "/projects/JCESR/pkcoff/qchem43.mod1"
        qcaux_path = "/projects/JCESR/qcaux"
        qc_exe_path = "/projects/JCESR/pkcoff/public/qcprog-43-optv5.exe"
        qc_scr_dir = "/dev/local/qchem"
        # on BGQ if the ramdisk is not used then use the file system for scrach space
        if use_ramdisk == False:
            qc_scr_dir = fs_scr_dir
        scr_size_bytes = int(scr_size_GB * (2 ** 30))

        # if not using the ramdisk then do not set QCLOCALFSSIZE
        if use_ramdisk == True:
            qc_envs = {
                "HOME": os.environ["HOME"],
                "QC": qc_package_path,
                "QCAUX": qcaux_path,
                "QCSCRATCH": qc_scr_dir,
                "QCFILEPREF": os.path.join(qc_scr_dir, "qc_job"),
                "QCTMPDIR": qc_scr_dir,
                "QCTHREADS": num_threads,
                "OMP_NUM_THREADS": num_threads,
                "QCLOCALFSSIZE": scr_size_bytes,
                "INT_OMP_MIN_LENSTEP_PATH1": 1,
                "INT_OMP_MAX_LENSTEP": 1,
                "INT_OMP_MAX_LENSTEP_PATH1": 1,
                "INT_OMP_MIN_LENSTEP": 1,
                "BLAS3_DFT_OMP": 1,
                "PAMI_CLIENTS": "MPI,SharedCounter"
            }
        else:
            qc_envs = {
                "HOME": os.environ["HOME"],
                "QC": qc_package_path,
                "QCAUX": qcaux_path,
                "QCSCRATCH": qc_scr_dir,
                "QCFILEPREF": os.path.join(qc_scr_dir, "qc_job"),
                "QCTMPDIR": qc_scr_dir,
                "QCTHREADS": num_threads,
                "OMP_NUM_THREADS": num_threads,
                "INT_OMP_MIN_LENSTEP_PATH1": 1,
                "INT_OMP_MAX_LENSTEP": 1,
                "INT_OMP_MAX_LENSTEP_PATH1": 1,
                "INT_OMP_MIN_LENSTEP": 1,
                "BLAS3_DFT_OMP": 1,
                "PAMI_CLIENTS": "MPI,SharedCounter"
            }

        with open(os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               "data/qc_lenstep_config_set"))) \
                as f:
            lenstep_config_set = json.load(f)
        lenstep_setting = lenstep_config_set["small size"]
        qc_envs.update(lenstep_setting)
        qsub_env_text = ":".join(
            ["{}={}".format(k, v) for k, v in qc_envs.items()])
        runjob_env_text = " ".join(
            ["--envs {}={}".format(k, v) for k, v in qc_envs.items()])
        block_name = os.environ["COBALT_PARTNAME"]
        cur_dir = os.getcwd()
        qc_runjob_cmd = "runjob --block {block_name} -n {num_processes} --ranks-per-node " \
                        "{ranks_per_node} --verbose 2 {envs} --cwd={cur_dir} : {qc_exe_path} {qc_input_file} " \
                        "{scr_dir}".format(
            block_name=block_name, num_processes=num_nodes * ranks_per_node,
            ranks_per_node=ranks_per_node, envs=runjob_env_text,
            cur_dir=cur_dir,
            qc_exe_path=qc_exe_path, qc_input_file=input_file,
            scr_dir=qc_scr_dir)
        qc_qsub_cmd = "qsub -A JCESR -t {max_minutes} -n {num_nodes} --mode c{ranks_per_node} --env {envs} " \
                      "--cwd={cur_dir} {qc_exe_path} {qc_input_file} {scr_dir}".format(
            max_minutes=max_minutes, num_nodes=num_nodes,
            ranks_per_node=ranks_per_node,
            envs=qsub_env_text, cur_dir=cur_dir, qc_exe_path=qc_exe_path,
            qc_input_file=input_file,
            scr_dir=qc_scr_dir)
        if use_runjob:
            qc_cmd = qc_runjob_cmd
        else:
            qc_cmd = qc_qsub_cmd
        return qc_cmd

    @staticmethod
    def _customize_alcf_qcinp(qcinp, num_nodes=8):
        for qj in qcinp.jobs:
            qj.params["rem"]["parallel_tasks"] = num_nodes
            if qj.params["rem"]["jobtype"] == "freq":
                qj.params["rem"]["ideriv"] = 1
            qj.params["rem"]["PDIAG_ON"] = 1
            solvent_params = qj.params.pop("pcm_solvent", None)
            if solvent_params is not None:
                qj.params["solvent"] = solvent_params
        # use Paul Coffman's version of pathtable
        pathtable_src = "/projects/JCESR/pkcoff/public/pathtable"
        shutil.copy(pathtable_src, "pathtable")

    @classmethod
    def get_alcf_qchem_cmd(cls, qcinp):
        prev_qchem_dir = os.getcwd()
        parent_qchem_dir = os.path.abspath(os.path.join(prev_qchem_dir, os.pardir))
        # These are variables used to support the usage of a file system scratch directory
        # on BGQ for pcm solvation jobs which use too much disk space for the ramdisk
        # Put the scratch directory in ../scratch and keep the names as short as possible
        # since the qchem binary has 114 char max for the full path to the scratch save dir
        use_fs_scr = False
        filesys_scr_dir = parent_qchem_dir + "/scratch/s" + ''.join(
            random.choice(string.ascii_lowercase) for i in range(10))

        # ALCF, Blue Gene
        num_nodes = 4
        num_threads = 32
        scr_size_GB = 0.9313
        max_minutes = 8 * 60

        # For BGQ for sp use a file system scratch directory instead of the ramdisk
        # because the amount of storage required can easily exceed the ramdisk for the pcm step

        for qj in qcinp.jobs:
            if qj.params["rem"]["jobtype"] == "sp":
                use_fs_scr = True
        if use_fs_scr == True:
            # for BGQ hang a short savedir name off the scratch dir
            save_dir = filesys_scr_dir + "/sd"
            os.makedirs(save_dir)
            qc_exe = shlex.split(cls._calibrate_alcf_cmd(
                num_nodes=num_nodes, max_minutes=max_minutes,
                num_threads=num_threads,
                scr_size_GB=scr_size_GB, use_ramdisk=False, fs_scr_dir=save_dir))
            half_cpus_cmd = shlex.split(cls._calibrate_alcf_cmd(
                num_nodes=num_nodes, max_minutes=max_minutes,
                num_threads=num_threads / 2,
                scr_size_GB=scr_size_GB, use_ramdisk=False, fs_scr_dir=save_dir))

        else:
            qc_exe = shlex.split(cls._calibrate_alcf_cmd(
                num_nodes=num_nodes, max_minutes=max_minutes,
                num_threads=num_threads,
                scr_size_GB=scr_size_GB))
            half_cpus_cmd = shlex.split(cls._calibrate_alcf_cmd(
                num_nodes=num_nodes, max_minutes=max_minutes,
                num_threads=num_threads / 2,
                scr_size_GB=scr_size_GB))
        cls._customize_alcf_qcinp(qcinp, num_nodes=num_nodes)
        openmp_cmd = None
        return qc_exe, half_cpus_cmd, openmp_cmd

    @classmethod
    def alcf_clean_up(cls, qcinp):
        prev_qchem_dir = os.getcwd()
        parent_qchem_dir = os.path.abspath(os.path.join(prev_qchem_dir, os.pardir))
        # These are variables used to support the usage of a file system scratch directory
        # on BGQ for pcm solvation jobs which use too much disk space for the ramdisk
        # Put the scratch directory in ../scratch and keep the names as short as possible
        # since the qchem binary has 114 char max for the full path to the scratch save dir
        use_fs_scr = False
        filesys_scr_dir = parent_qchem_dir + "/scratch/s" + ''.join(
            random.choice(string.ascii_lowercase) for i in range(10))
        for qj in qcinp.jobs:
            if qj.params["rem"]["jobtype"] == "sp":
                use_fs_scr = True
                # on BGQ if we used the filesystem for the scratch directory clean it up
        if 'vesta' in socket.gethostname():
            if use_fs_scr == True:
                shutil.rmtree(filesys_scr_dir)
