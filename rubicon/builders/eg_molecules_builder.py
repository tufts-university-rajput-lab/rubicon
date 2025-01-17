# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
Build molecules collection
Adapted from Dan Gunter and Wei Chen's vasp materials builder
"""

import copy
import datetime
import logging
import re
import sys

import pymongo
from pymongo import ASCENDING

from rubicon.builders import eg_shared
from rubicon.utils.qchem_firework_creator import QChemFireWorkCreator
from six.moves import map
from six.moves import zip

__author__ = "Xiaohui Qu"
__copyright__ = "Copyright 2012-2013, The Electrolyte Genome Project"
__version__ = "1.0"
__maintainer__ = "Xiaohui Qu"
__email__ = "xqu@lbl.gov"
__status__ = "Development"
__date__ = "1/1/14"

_log = logging.getLogger('eg.' + __name__)


class TaskKeys:
    """Keys we need to project from task collection to do
       the work of building the materials collection.
    """

    def __init__(self):
        pass

    fields = (
        'task_id', 'snlgroup_id_final', 'inchi_final', 'task_type', 'elements',
        'can', 'smiles', 'charge', 'spin_multiplicity', 'implicit_solvent',
        'user_tags', 'run_tags', 'snl_final', 'task_id', "molecule_final",
        'nelements', 'reduced_cell_formula_abc', 'pretty_formula',
        'pointgroup', 'inchi_root',
        'calculations.scf.energies', 'calculations.scf_pcm.energies',
        'calculations.scf_sm12mk.energies',
        'formula', 'task_id_deprecated', 'svg', 'xyz', "last_updated")

    base_molecules = ("quinoxaline", "anthrachinon", "thiane", "viologen")
    lei_1_group_pattern = re.compile('(?P<base_mol>\w+)_wfs_(?P<position>\d+)_'
                                     '(?P<group_name>\w+)')
    lei_2_group_pattern = re.compile('(?P<base_mol>\w+)_(?P<pos1>\d+)_'
                                     '(?P<group1>\w+)_(?P<pos2>\d+)(?P<group2>\w+)')
    literal_to_formula_group_name = {"nitro": "-NO2", "cyano": "-CN",
                                     "trichloromethyl": "-CCl3",
                                     "carboxyl": "-COOH",
                                     "fluoro": "-F", "ethynyl": "-CCH",
                                     "methyl": "-CH3", "ethyl": "-CH2CH3",
                                     "hydroxyl": "-OH", "vinyl": "-C=CH2",
                                     "methoxyl": "-OCH3",
                                     "ethanamide": "-NHC(O)CH3",
                                     "benzene": "-C5H6", "amine": "-NH2",
                                     "methylamine": "-NHCH3",
                                     "dimethylamine": "-N(CH3)2",
                                     "Aceto": "-C(O)CH3", "Amine": "-NH2",
                                     "Benzene": "-C5H6", "Butyl": "-C4H9",
                                     "Carboxyl": "-COOH", "Chloro": "-Cl",
                                     "Cyano": "-CN",
                                     "Dimethylamine": "-N(CH3)2",
                                     "Ethanamide": "-NHC(O)CH3",
                                     "Ethyl": "-CH2CH3", "Ethynyl": "-CCH",
                                     "Fluoro": "-F",
                                     "G0": "-OCH3", "G1": "-OCH2CH2OCH3",
                                     "G2": "-O(CH2CH2O)2CH3",
                                     "G3": "-O(CH2CH2O)3CH3",
                                     "Hydroxyl": "-OH", "Methoxyl": "-OCH3",
                                     "Methyl": "-CH3",
                                     "Methylamine": "-NHCH3", "Nitro": "-NO2",
                                     "S2": "-OCH2OCH3", "S3": "-O(CH2O)2CH3",
                                     "Tribromomethyl": "-CBr3",
                                     "Tricholoromethyl": "-CCl3",
                                     "Triflouromethyl": "-CF3",
                                     "Vinyl": "-C=CH2", "butyl": "-C4H9"}


class MoleculesBuilder(eg_shared.ParallelBuilder):
    """Build derived 'molecules' collection.
    """

    # absolute electrode potentials of some metals/molecules
    ref_potentials = {
        'hydrogen': 4.44,
        'magnesium': 2.07,
        'lithium': 1.40}

    def __init__(self, collections, **kwargs):
        """Create new molecules builder.

        Args:
            collections: Set of connected DB collections
                Type: eg_shared.Collections
        """
        eg_shared.ParallelBuilder.__init__(self, **kwargs)
        self._c = collections
        self._c.molecules.remove()
        self.ref_charge = 0
        self.ref_charge_range = (-1, 0, 1)
        logging.basicConfig(level=logging.INFO)
        _log.setLevel(logging.INFO)
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setLevel(getattr(logging, 'INFO'))
        _log.addHandler(sh)

    def run(self):
        """Run the builder.
        """
        sss = []
        for ch in self.ref_charge_range:
            self.ref_charge = ch
            charge_state = QChemFireWorkCreator.get_state_name(ch, 1).split()[
                1]
            _log.info(
                "Getting distinct root INCHIs for {}s".format(charge_state))
            if ch == 0:
                spec = {"$or": [{"user_tags.initial_charge": 0},
                                {"user_tags.initial_charge": {
                                    "$exists": False}}]}
            else:
                spec = {"user_tags.initial_charge": ch}
            inchi_root = list(self._c.tasks.find(filter=spec,
                                                 projection='inchi_root')
                              .distinct('inchi_root'))
            _log.info(
                "There are total {} unique INCHIs".format(len(inchi_root)))
            list(map(self.add_item, inchi_root))
            _log.info("Beginning analysis")
            states = self.run_parallel()
            sss.extend(states)
            self._build_indexes()
        return self.combine_status(sss)

    def process_item(self, inchi_root):
        """Create and add material for a given grouping identifer.
        """
        query = {'state': 'successful', 'inchi_root': inchi_root,
                 'task_type': "single point energy"}
        solvents = self._c.tasks.find(filter=query,
                                      projection=TaskKeys.fields).distinct(
            "implicit_solvent.solvent_name"
        )
        # solvent_model = "ief-pcm"
        solvent_models = self._c.tasks.find(filter=query,
                                            projection=TaskKeys.fields).distinct(
            "implicit_solvent.model")
        molecule = dict()
        molecule['charge'] = self.ref_charge
        docs_available = False
        molecule['solvated_properties'] = dict()
        for solvent_model in solvent_models:
            for solvent in solvents:
                query['implicit_solvent.solvent_name'] = solvent
                query['implicit_solvent.model'] = solvent_model
                docs = list(self._c.tasks.find(filter=query,
                                               projection=TaskKeys.fields))
                if docs:
                    docs_available = True
                d = self.build_molecule_solvated_properties(docs)
                if d and len(d) > 0:
                    solvent_key = "{}_{}".format(solvent,
                                                 solvent_model).replace(".",
                                                                        "_")
                    molecule['solvated_properties'][solvent_key] = d
        if not docs_available:
            return 1
        if len(molecule['solvated_properties']) == 0:
            return 2
        del query['implicit_solvent.solvent_name']
        d = self.build_molecule_vacuum_properties(copy.deepcopy(query))
        if d and len(d) > 0:
            molecule['vacuum_properties'] = d
        else:
            return 2
        query['charge'] = self.ref_charge
        docs = self._c.tasks.find_one(filter=query, projection=TaskKeys.fields)
        if not docs:
            return 1
        d = self.build_molecule_common_properties(docs)
        if d and len(d) > 0:
            molecule.update(d)
        else:
            return 2
        d = self.build_molecule_structure_properties(docs)
        if d and len(d) > 0:
            molecule.update(d)
        molecule['created_at'] = datetime.datetime.now()
        molecule['updated_at'] = datetime.datetime.now()
        self._insert_molecule(molecule)
        return 0

    def build_molecule_ipea(self, docs, molecule, solution_phase=True):
        molecule["task_id"] = dict()
        molecule["task_id_deprecated"] = dict()
        molecule["snlgroup_id_final"] = dict()
        molecule["charge"] = dict()
        molecule["spin_multiplicity"] = dict()
        molecule["snl_final"] = dict()
        molecule["molecule"] = dict()
        molecule["xyz"] = dict()
        molecule["inchi"] = dict()
        molecule["can"] = dict()
        molecule["smiles"] = dict()
        for k in docs.keys():
            molecule["task_id"][k] = docs[k]["task_id"]
            molecule["task_id_deprecated"][k] = docs[k]["task_id_deprecated"]
            molecule["snlgroup_id_final"][k] = docs[k]["snlgroup_id_final"]
            molecule["charge"][k] = docs[k]["charge"]
            molecule["spin_multiplicity"][k] = docs[k]["spin_multiplicity"]
            molecule["snl_final"][k] = docs[k]["snl_final"]
            molecule["molecule"][k] = docs[k]["molecule_final"]
            molecule["xyz"][k] = docs[k]["xyz"]
            molecule["inchi"][k] = docs[k]["inchi_final"]
            molecule["can"][k] = docs[k]["can"]
            molecule["smiles"][k] = docs[k]["smiles"]
        if solution_phase:
            # get the solution phase scf key name, scf_pcm, scf_sm12mk, etc.
            scf_all = set(docs["neutral"]["calculations"].keys())
            scf_all.remove('scf')
            scf_name = scf_all.pop()
        else:
            scf_name = 'scf'
        if "cation" in docs:
            molecule["IP"] = \
                docs["cation"]["calculations"][scf_name]["energies"][-1][-1] \
                - \
                docs["neutral"]["calculations"][scf_name]["energies"][-1][-1]
        if "anion" in docs:
            molecule["EA"] = \
                docs["neutral"]["calculations"][scf_name]["energies"][-1][-1] \
                - \
                docs["anion"]["calculations"][scf_name]["energies"][-1][-1]
        if "IP" in molecule and "EA" in molecule:
            molecule["electrochemical_window_width"] = molecule["IP"] - \
                                                       molecule["EA"]
        molecule['electrode_potentials'] = dict()
        if solution_phase:
            if 'IP' in molecule:
                molecule['electrode_potentials']['oxidation'] = dict()
                for electrode in self.ref_potentials.keys():
                    molecule['electrode_potentials']['oxidation'][electrode] \
                        = molecule['IP'] - self.ref_potentials[electrode]
            if 'EA' in molecule:
                molecule['electrode_potentials']['reduction'] = dict()
                for electrode in self.ref_potentials.keys():
                    molecule['electrode_potentials']['reduction'][electrode] \
                        = molecule['EA'] - self.ref_potentials[electrode]
        return scf_name

    def build_molecule_solvated_properties(self, taskdocs):
        docs = dict()
        for td in taskdocs:
            if td["charge"] == self.ref_charge:
                if "neutral" not in docs:
                    docs["neutral"] = td
                else:
                    if td["last_updated"] > docs["neutral"]["last_updated"]:
                        docs["neutral"] = td
                    else:
                        continue
            elif td["charge"] == self.ref_charge + 1:
                if "cation" not in docs:
                    docs["cation"] = td
                else:
                    if td["last_updated"] > docs["cation"]["last_updated"]:
                        docs["cation"] = td
                    else:
                        continue
            elif td["charge"] == self.ref_charge - 1:
                if "anion" not in docs:
                    docs["anion"] = td
                else:
                    if td["last_updated"] > docs["anion"]["last_updated"]:
                        docs["anion"] = td
                    else:
                        continue
        if len(docs) < 2 or ("neutral" not in docs):
            return None
        molecule = dict()
        scf = self.build_molecule_ipea(docs, molecule, solution_phase=True)
        molecule['solvation_energy'] = docs["neutral"]["calculations"]["scf"][
                                           "energies"][-1][-1] - \
                                       docs["neutral"]["calculations"][scf][
                                           "energies"][-1][-1]
        molecule["implicit_solvent"] = copy.deepcopy(docs['neutral'][
                                                         "implicit_solvent"])
        return molecule

    def build_molecule_vacuum_properties(self, query):
        docs = dict()
        for c, i in zip(["anion", "neutral", "cation"], [-1, 0, 1]):
            query['charge'] = self.ref_charge + i
            taskdocs = self._c.tasks.find(filter=query,
                                          projection=TaskKeys.fields,
                                          sort=[("_id", pymongo.DESCENDING)])
            if taskdocs.count() == 0:
                continue
            docs[c] = taskdocs[0]
        if len(docs) < 2 or ("neutral" not in docs):
            return None
        molecule = dict()
        self.build_molecule_ipea(docs, molecule, solution_phase=False)
        return molecule

    @staticmethod
    def build_molecule_common_properties(docs):
        """Transforms task document to molecules document.
        """
        molecule = dict()
        molecule["inchi_root"] = docs["inchi_root"]
        molecule["elements"] = copy.deepcopy(docs["elements"])
        molecule["nelements"] = docs["nelements"]
        molecule["user_tags"] = copy.deepcopy(docs["user_tags"])
        molecule["run_tags"] = copy.deepcopy(docs["run_tags"])
        molecule["reduced_cell_formula_abc"] = docs["reduced_cell_formula_abc"]
        molecule["pretty_formula"] = docs["pretty_formula"]
        molecule["formula"] = docs["formula"]
        molecule["pointgroup"] = docs["pointgroup"]
        molecule["svg"] = docs["svg"]
        return molecule

    @staticmethod
    def parse_derivation_name(derivation_name, molname):
        group_name_texts, base_mol_name = derivation_name.split(';')
        if base_mol_name == "thiophene" and "thiophen" not in molname:
            base_mol_name = "tetrathiafulvulene"
        functional_groups = []
        if base_mol_name != 'dbbb':
            for g in group_name_texts.split(','):
                position, gn = g.split('-')
                functional_groups.append(
                    TaskKeys.literal_to_formula_group_name[gn])
        else:
            pos_text, group1_text, group2_text = group_name_texts.split(',')
            group1_name = group1_text.split('-')[1]
            group2_name = group2_text.split('-')[1]
            functional_groups.append(
                TaskKeys.literal_to_formula_group_name[group1_name])
            functional_groups.append(
                TaskKeys.literal_to_formula_group_name[group2_name])
        return functional_groups, base_mol_name

    @staticmethod
    def build_molecule_structure_properties(docs):
        molecule = dict()
        if "molname" in docs["user_tags"]:
            if "derivation_name" in docs["user_tags"] and \
                            ";" in docs["user_tags"]["derivation_name"]:
                functional_groups, base_mol_name = MoleculesBuilder.parse_derivation_name(
                    docs["user_tags"]["derivation_name"],
                    docs["user_tags"]['molname'])
                molecule["base_molecule"] = base_mol_name
                molecule["functional_groups"] = sorted(functional_groups)
            else:
                molname = docs["user_tags"]['molname']
                m = TaskKeys.lei_1_group_pattern.search(molname)
                if m:
                    base_mol = m.group("base_mol")
                    literal_group_name = m.group("group_name")
                    if literal_group_name in TaskKeys.literal_to_formula_group_name:
                        formula_group_name = \
                            TaskKeys.literal_to_formula_group_name[
                                literal_group_name]
                        molecule["base_molecule"] = base_mol
                        molecule["functional_groups"] = [formula_group_name]
                m = TaskKeys.lei_2_group_pattern.search(molname)
                if m:
                    base_mol = m.group("base_mol")
                    literal_group1 = m.group("group1")
                    literal_group2 = m.group("group2")
                    if literal_group1 in TaskKeys.literal_to_formula_group_name \
                            and literal_group2 in TaskKeys.literal_to_formula_group_name:
                        formula_group1 = \
                            TaskKeys.literal_to_formula_group_name[
                                literal_group1]
                        formula_group2 = \
                            TaskKeys.literal_to_formula_group_name[
                                literal_group2]
                        molecule["base_molecule"] = base_mol
                        molecule["functional_groups"] = sorted(
                            [formula_group1, formula_group2])
        return molecule

    def _build_indexes(self):
        self._c.molecules.ensure_index(
            [('inchi_root', ASCENDING), ('charge', ASCENDING)], unique=True)
        for key in ['inchi_root', 'charge', 'nelements', 'elements',
                    'reduced_cell_formula', 'pretty_formula']:
            _log.info("Building {} index".format(key))
            self._c.molecules.ensure_index(key)
        _log.info("Building nelements and elements compound index")
        compound_index = [('nelements', ASCENDING), ('elements', ASCENDING)]
        self._c.molecules.ensure_index(compound_index)

    def _insert_molecule(self, doc):
        """All database insertion should be done from this method
        """
        _log.info("Inserting Material with InChI {i}, ".
                  format(i=str(doc['inchi_root'])))
        self._c.molecules.insert(doc)
