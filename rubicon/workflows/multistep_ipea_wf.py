import copy
import itertools
import json
import os
from fireworks.core.firework import Workflow, FireWork, Tracker
from pymatgen.io.babelio import BabelMolAdaptor
from pymatgen.io.qchemio import QcInput, QcTask
from rubicon.dupefinders.dupefinder_eg import DupeFinderEG
from rubicon.firetasks.qchem_task import QChemTask

__author__ = 'xiaohuiqu'


class QChemFireWorkCreator():
    def __init__(self, mol, molname, mission, additional_user_tags=None,
                 dupefinder=None, priority=1, update_spec=None, large=False):
        self.molname = molname
        self.mol = mol
        self.large = large
        initial_inchi = self.get_inchi(mol)
        user_tags = {'mission': mission,
                     "molname": molname}
        if additional_user_tags:
            user_tags.update(additional_user_tags)
        spec = dict()
        spec['user_tags'] = user_tags
        spec['_priority'] = priority
        spec['_dupefinder'] = dupefinder.to_dict() if dupefinder \
            else DupeFinderEG().to_dict()
        tracker_out = Tracker("mol.qcout", nlines=20)
        tracker_std = Tracker("mol.qclog", nlines=10)
        tracker_joberr = Tracker("FW_job.error", nlines=20)
        tracker_jobout = Tracker("FW_job.out", nlines=20)
        spec["_trackers"] = [tracker_out, tracker_std, tracker_joberr,
                             tracker_jobout]
        spec['run_tags'] = dict()
        spec['implicit_solvent'] = {}
        spec['inchi'] = initial_inchi
        spec['num_atoms'] = len(mol)
        if update_spec:
            spec.update(update_spec)
        self.base_spec = lambda: copy.deepcopy(spec)

    @staticmethod
    def get_inchi(mol):
        bb = BabelMolAdaptor(mol)
        pbmol = bb.pybel_mol
        return pbmol.write("inchi").strip()

    @staticmethod
    def get_state_name(charge, spin_multiplicity):
        charge_state = {-2: "anion_2", -1: "anion", 0: "neutral", 1: "cation", 2: "cation_2"}
        spin_state = {1: "singlet", 2: "doublet", 3: "triplet"}
        return spin_state[spin_multiplicity] + " " + charge_state[charge]

    @staticmethod
    def get_exchange_correlation_basis_auxbasis_remparams(method):
        if method is None:
            method = "B3LYP/6-31+G*"
        aux_basis = None
        correlation = None
        rem_params = None
        theoretical_level, basis_set = method.split('/')
        basis_set = basis_set.lower()
        if theoretical_level.lower() == "b3lyp-xdm":
            exchange = 'b3lyp'
            rem_params = {"dftvdw_jobnumber": 1,
                          "dftvdw_method": 1,
                          "dftvdw_print": 1,
                          "dftvdw_kai": 800,
                          "dftvdw_use_ele_drv": 1}
        elif theoretical_level.lower() == "xygjos":
            exchange = "xygjos"
            if basis_set == "6-31+g*":
                aux_basis = "rimp2-aug-cc-pvdz"
            else:
                aux_basis = "rimp2-aug-cc-pvtz"
        elif theoretical_level.lower() == "pbe-d3":
            exchange = 'pbe'
            correlation = 'pbe'
            rem_params = {"dft_d": "empirical_grimme3",
                          "dft_d3_s6": 1000,
                          "dft_d3_rs6": 1217,
                          "dft_d3_s8": 722,
                          "dft_d3_3body": False}
        else:
            exchange = theoretical_level.lower()
        method_token = [t for t in [basis_set, exchange, aux_basis, correlation, rem_params]
                        if t]
        return exchange, correlation, basis_set,  aux_basis, rem_params, method_token


    def geom_fw(self, charge, spin_multiplicity, fw_id_cal, fw_id_db,
                priority=None, method=None):
        task_type = "geometry optimization"
        state_name = self.get_state_name(charge, spin_multiplicity)
        if not method:
            if self.large:
                method = "PBE-D3/6-31+G*"
            else:
                method = "B3lYP/6-31+G*"
        title = self.molname + " " + state_name + " " + method + " " + task_type
        exchange, correlation, basis_set,  aux_basis, rem_params, method_token = self.\
            get_exchange_correlation_basis_auxbasis_remparams(method)
        qctask = QcTask(self.mol, charge=charge, spin_multiplicity=spin_multiplicity,
                        jobtype="opt", title=title, exchange=exchange, correlation=correlation,
                        basis_set=basis_set, aux_basis_set=aux_basis, rem_params=rem_params)
        if self.large:
            qctask.set_geom_max_iterations(200)
            qctask.set_scf_algorithm_and_iterations(iterations=100)
            qctask.scale_geom_opt_threshold(gradient=1.0,
                                            displacement=10.0,
                                            energy=10.0)
        qcinp = QcInput([qctask])
        spec = self.base_spec()
        spec["qcinp"] = qcinp.to_dict
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        spec['run_tags']['methods'] = method_token
        if priority:
            spec['_priority'] = priority
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemGeomOptDBInsertionTask
        fw_geom_cal = FireWork([QChemTask()],
                               spec=spec, name=task_name, fw_id=fw_id_cal)
        spec_db = copy.deepcopy(spec)
        del spec_db['_dupefinder']
        spec_db['_allow_fizzled_parents'] = True
        spec_db['task_type'] = task_type + ' DB Insertion'
        del spec_db["_trackers"][:2]
        task_name_db = task_name + " DB Insertion"
        fw_geom_db = FireWork([QChemGeomOptDBInsertionTask()],
                              spec=spec_db, name=task_name_db, fw_id=fw_id_db)

        return fw_geom_cal, fw_geom_db

    def freq_fw(self, charge, spin_multiplicity, fw_id_cal, fw_id_db,
                priority=None, method=None):
        if not method:
            if self.large:
                method = "PBE-D3/6-31+G*"
            else:
                method = "B3lYP/6-31+G*"
        task_type = "vibrational frequency"
        state_name = self.get_state_name(charge, spin_multiplicity)
        title = self.molname + " " + state_name + " " + method + " " + task_type
        exchange, correlation, basis_set,  aux_basis, rem_params, method_token = self. \
            get_exchange_correlation_basis_auxbasis_remparams(method)
        if exchange.lower() in ["xygjos"]:
            rem_params["IDERIV"] = 1
        qctask = QcTask(self.mol, charge=charge, spin_multiplicity=spin_multiplicity,
                        jobtype="freq", title=title, exchange=exchange, correlation=correlation,
                        basis_set=basis_set, aux_basis_set=aux_basis, rem_params=rem_params)
        if self.large:
            qctask.set_scf_algorithm_and_iterations(iterations=100)
        qcinp = QcInput([qctask])
        spec = self.base_spec()
        spec["qcinp"] = qcinp.to_dict
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        spec['run_tags']['methods'] = method_token
        if priority:
            spec['_priority'] = priority
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemFrequencyDBInsertionTask
        fw_freq_cal = FireWork([QChemTask()],
                               spec=spec, name=task_name, fw_id=fw_id_cal)
        spec_db = copy.deepcopy(spec)
        del spec_db['_dupefinder']
        spec_db['_allow_fizzled_parents'] = True
        spec_db['task_type'] = task_type + ' DB Insertion'
        del spec_db["_trackers"][:2]
        task_name_db = task_name + " DB Insertion"
        fw_freq_db = FireWork([QChemFrequencyDBInsertionTask()],
                              spec=spec_db, name=task_name_db, fw_id=fw_id_db)
        return fw_freq_cal, fw_freq_db

    @staticmethod
    def get_dielectric_constant(solvent):
        dielec_const_file = os.path.join(os.path.dirname(__file__),
                                         "../utils/data", "dielectric_constants.json")
        with open(dielec_const_file) as f:
            dielec_consts = json.load(f)
        if solvent not in dielec_consts:
            raise Exception("Don't know the dielectric constants for "
                            "solvent '{}'".format(solvent))
        dielectric_constant = dielec_consts[solvent]
        return dielectric_constant

    @staticmethod
    def get_smx_solvent(implicit_solvent, solvent, solvent_method):
        smx_data_file = os.path.join(os.path.dirname(__file__),
                                     "../utils/data", "smx_data.json")
        with open(smx_data_file) as f:
            smx_data = json.load(f)
        if solvent not in smx_data["builtin_solvent"]:
            smx_solvent = "other"
            if solvent not in smx_data["custom_solvent"]:
                raise Exception("Don't know the SMx parameters for "
                                "solvent '{}'".format(solvent))
            implicit_solvent['solvent_data'] = \
                smx_data["custom_solvent"][solvent]
        else:
            smx_solvent = solvent
        rem_options = dict(solvent_method=solvent_method.lower(),
                           smx_solvent=smx_solvent.lower())
        return rem_options, smx_solvent

    def set_solvent_method(self, qctask_sol, solvent, solvent_method):
        implicit_solvent = dict()
        implicit_solvent['solvent_name'] = solvent
        if solvent_method.lower() in ['cpcm', 'ief-pcm']:
            if solvent_method.lower() == 'ief-pcm':
                solvent_theory = 'ssvpe'
            else:
                solvent_theory = 'cpcm'
            dielectric_constant = self.get_dielectric_constant(solvent)
            qctask_sol.use_pcm(solvent_params={"Dielectric": dielectric_constant},
                               pcm_params={'Theory': solvent_theory})
            implicit_solvent['model'] = solvent_method.lower()
            implicit_solvent['dielectric_constant'] = dielectric_constant
            implicit_solvent['radii'] = 'uff'
            implicit_solvent['vdwscale'] = 1.1
        elif solvent_method.lower() == 'cosmo':
            dielectric_constant = self.get_dielectric_constant(solvent)
            qctask_sol.use_cosmo(dielectric_constant)
            implicit_solvent['model'] = 'cosmo'
            implicit_solvent['dielectric_constant'] = dielectric_constant
        elif solvent_method.lower() in ['sm12mk', 'sm12chelpg', 'sm12',
                                        'sm8']:
            implicit_solvent['model'] = solvent_method.lower()
            rem_options, smx_solvent = self.get_smx_solvent(implicit_solvent, solvent, solvent_method)
            qctask_sol.params['rem'].update(rem_options)
            implicit_solvent['smx_solvent'] = smx_solvent
        else:
            raise Exception("Don't know how to setup solvent model '{}'".
                            format(solvent_method))
        if not self.large:
            qctask_sol.set_dft_grid(128, 302)
            qctask_sol.set_integral_threshold(12)
            qctask_sol.set_scf_convergence_threshold(8)
        else:
            qctask_sol.set_scf_algorithm_and_iterations(iterations=100)
            if solvent_method.lower() in ['cpcm', 'ief-pcm']:
                qctask_sol.params['pcm'].update({"hpoints": 194,
                                                 "heavypoints": 194})
        return implicit_solvent

    def sp_fw(self, charge, spin_multiplicity, fw_id_cal, fw_id_db,
              solvent_method="ief-pcm", solvent="water", priority=None, qm_method=None,
              population_method=None, task_type_name=None):
        if not qm_method:
            qm_method = "B3LYP/6-31+G*"
        spec = self.base_spec()
        if priority:
            spec['_priority'] = priority
        task_type = task_type_name if task_type_name else "single point energy"
        state_name = self.get_state_name(charge, spin_multiplicity)
        title = self.molname + " " + state_name + " " + qm_method + " " + task_type
        title += "\n Gas Phase"
        exchange, correlation, basis_set,  aux_basis, rem_params, method_token = self. \
            get_exchange_correlation_basis_auxbasis_remparams(qm_method)
        if population_method:
            if not rem_params:
                rem_params = dict()
            if population_method.lower() == "nbo":
                rem_params["nbo"] = 1
            elif population_method.lower() == "chelpg":
                rem_params["chelpg"] = True
        qctask_vac = QcTask(self.mol, charge=charge, spin_multiplicity=spin_multiplicity,
                            jobtype="sp", title=title, exchange=exchange, correlation=correlation,
                            basis_set=basis_set, aux_basis_set=aux_basis, rem_params=rem_params)
        if not self.large:
            qctask_vac.set_dft_grid(128, 302)
            qctask_vac.set_integral_threshold(12)
            qctask_vac.set_scf_convergence_threshold(8)
        else:
            qctask_vac.set_scf_algorithm_and_iterations(iterations=100)

        title = " Solution Phase, {}".format(solvent)
        qctask_sol = QcTask(self.mol, charge=charge, spin_multiplicity=spin_multiplicity,
                            jobtype="sp", title=title, exchange=exchange, correlation=correlation,
                            basis_set=basis_set, aux_basis_set=aux_basis, rem_params=rem_params)
        qctask_sol.set_scf_initial_guess(guess="read")
        implicit_solvent = self.set_solvent_method(qctask_sol, solvent, solvent_method)

        qcinp = QcInput([qctask_vac, qctask_sol])
        spec["qcinp"] = qcinp.to_dict
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        spec['run_tags']['methods'] = method_token
        spec['implicit_solvent'] = implicit_solvent
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemSinglePointEnergyDBInsertionTask
        fw_sp_cal = FireWork([QChemTask()],
                             spec=spec, name=task_name, fw_id=fw_id_cal)
        spec_db = copy.deepcopy(spec)
        del spec_db['_dupefinder']
        spec_db['_allow_fizzled_parents'] = True
        spec_db['task_type'] = task_type + ' DB Insertion'
        del spec_db["_trackers"][:2]
        task_name_db = task_name + " DB Insertion"
        fw_sp_db = FireWork([QChemSinglePointEnergyDBInsertionTask()],
                            spec=spec_db, name=task_name_db, fw_id=fw_id_db)
        return fw_sp_cal, fw_sp_db

    def vacuum_only_sp_fw(self, charge, spin_multiplicity, fw_id_cal, fw_id_db,
                          priority=None, qm_method=None, population_method=None):
        if not qm_method:
            qm_method = "B3LYP/6-31+G*"
        spec = self.base_spec()
        if priority:
            spec['_priority'] = priority
        task_type = "vacuum only single point energy"
        state_name = self.get_state_name(charge, spin_multiplicity)
        title = self.molname + " " + state_name + " " + qm_method + " " + task_type
        title += "\n Gas Phase"
        exchange, correlation, basis_set,  aux_basis, rem_params, method_token = self. \
            get_exchange_correlation_basis_auxbasis_remparams(qm_method)
        if population_method:
            if not rem_params:
                rem_params = dict()
            if population_method.lower() == "nbo":
                rem_params["nbo"] = 1
            elif population_method.lower() == "chelpg":
                rem_params["chelpg"] = True
        qctask_vac = QcTask(self.mol, charge=charge, spin_multiplicity=spin_multiplicity,
                            jobtype="sp", title=title, exchange=exchange, correlation=correlation,
                            basis_set=basis_set, aux_basis_set=aux_basis, rem_params=rem_params)
        if not self.large:
            qctask_vac.set_dft_grid(128, 302)
            qctask_vac.set_integral_threshold(12)
            qctask_vac.set_scf_convergence_threshold(8)
        else:
            qctask_vac.set_scf_algorithm_and_iterations(iterations=100)

        qcinp = QcInput([qctask_vac])
        spec["qcinp"] = qcinp.to_dict
        spec['task_type'] = task_type
        spec['charge'] = charge
        spec['spin_multiplicity'] = spin_multiplicity
        spec['run_tags']['methods'] = method_token
        task_name = self.molname + ' ' + state_name + ' ' + task_type
        from rubicon.firetasks.multistep_qchem_task \
            import QChemSinglePointEnergyDBInsertionTask
        fw_sp_cal = FireWork([QChemTask()],
                             spec=spec, name=task_name, fw_id=fw_id_cal)
        spec_db = copy.deepcopy(spec)
        del spec_db['_dupefinder']
        spec_db['_allow_fizzled_parents'] = True
        spec_db['task_type'] = task_type + ' DB Insertion'
        del spec_db["_trackers"][:2]
        task_name_db = task_name + " DB Insertion"
        fw_sp_db = FireWork([QChemSinglePointEnergyDBInsertionTask()],
                            spec=spec_db, name=task_name_db, fw_id=fw_id_db)
        return fw_sp_cal, fw_sp_db


def multistep_ipea_fws(mol, name, mission, solvent, solvent_method, ref_charge, spin_multiplicities=(2, 1, 2), dupefinder=None, priority=1,
                       parent_fwid=None, additional_user_tags=None, qm_method=None):
    large = False
    if len(mol) > 50:
        large = True
    energy_method, geom_method = qm_method.split("//") if qm_method else (None, None)
    fw_creator = QChemFireWorkCreator(
        mol=mol, molname=name, mission=mission, dupefinder=dupefinder, priority=priority, large=large,
        additional_user_tags=additional_user_tags)
    fwid_base = 1
    if parent_fwid:
        if not (isinstance(parent_fwid, int) or isinstance(parent_fwid, list)):
            raise ValueError("Parent FireWork ID must be integer or list")
        parent_fwid = parent_fwid if isinstance(parent_fwid, list) \
            else [parent_fwid]
        fwid_base = max(parent_fwid) + 1
    fireworks = []
    links_dict = dict()
    # the task in the order of anion, neutral, cation
    cgi_cal, ngi_cal, agi_cal = (None, None, None)
    cfi_db, nfi_db, afi_db = (None, None, None)
    cgi_db, ngi_db, agi_db = (None, None, None)
    charges = [ref_charge + i for i in (-1, 0, 1)]
    if len(mol) > 1:
        fw_ids = zip(* [iter(range(fwid_base + 0, fwid_base + 6))] * 2)
        fws = (fw_creator.geom_fw(ch, spin, fwid_cal, fwid_db, method=geom_method)
               for ch, spin, (fwid_cal, fwid_db)
               in zip(charges, spin_multiplicities, fw_ids))
        (cgi_cal, cgi_db), (ngi_cal, ngi_db), (agi_cal, agi_db) = fw_ids
        fireworks.extend(itertools.chain.from_iterable(fws))
        links_dict.update(dict(fw_ids))

        if not large:
            fw_ids = zip(* [iter(range(fwid_base + 6, fwid_base + 6 + 6))] * 2)
            fws = (fw_creator.freq_fw(ch, spin, fwid_cal, fwid_db, method=geom_method)
                   for ch, spin, (fwid_cal, fwid_db)
                   in zip(charges, spin_multiplicities, fw_ids))
            (cfi_cal, cfi_db), (nfi_cal, nfi_db), (afi_cal, afi_db) = fw_ids
            fireworks.extend(itertools.chain.from_iterable(fws))
            links_dict.update(dict(fw_ids))
            links_dict.update({cgi_db: cfi_cal,
                               ngi_db: nfi_cal,
                               agi_db: afi_cal})

    fw_ids = zip(* [iter(range(fwid_base + 12, fwid_base + 12 + 6))] * 2)
    fws = (fw_creator.sp_fw(ch, spin, fwid_cal, fwid_db, solvent=solvent, solvent_method=solvent_method,
                            qm_method=energy_method)
           for ch, spin, (fwid_cal, fwid_db)
           in zip(charges, spin_multiplicities, fw_ids))
    (cspi_cal, cspi_db), (nspi_cal, nspi_db), (aspi_cal, aspi_db) = fw_ids
    links_dict.update(dict(fw_ids))
    fireworks.extend(itertools.chain.from_iterable(fws))
    if len(mol) > 1:
        if large:
            links_dict.update({cgi_db: cspi_cal, ngi_db: nspi_cal,
                               agi_db: aspi_cal})
        else:
            links_dict.update({cfi_db: cspi_cal, nfi_db: nspi_cal,
                               afi_db: aspi_cal})
        links_dict.update({nspi_db: [cgi_cal, agi_cal]})
        if parent_fwid:
            for pfw_id in parent_fwid:
                links_dict[pfw_id] = ngi_cal
    else:
        links_dict.update({nspi_db: [cspi_cal, aspi_cal]})
        if parent_fwid:
            for pfw_id in parent_fwid:
                links_dict[pfw_id] = nspi_cal
    return fireworks, links_dict


def mol_to_ipea_wf(mol, name, **kwargs):
    fireworks, links_dict = multistep_ipea_fws(mol, name, **kwargs)
    return Workflow(fireworks, links_dict, name)