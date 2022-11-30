#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2021 Rui Ribeiro
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
This python wrap-up of the GOMoDO webserver
"""
__author__ = "Rui Ribeiro"
__email__ = "rui.ribeiro@univr.it"
import warnings, os , sys
warnings.filterwarnings('ignore')
homedirectory=os.path.dirname(__file__)
sys.path.append(homedirectory)

class Tethering:
    """Tethered minimization on the MCS of ligands and a reference 3D molecular structure"""
    def MolAlignment(**kwargs):
        """:parameter reference: (str) reference 3D molecule
            :parameter target_name: (str) name of ligand
            :parameter target_smile: (str) smiles string of ligand"""
        refFile = kwargs.pop('ref')
        targetMolID = kwargs.pop('target_name')
        targetMolSMI = kwargs.pop('target_smiles')
        
        from openbabel import pybel
        from rdkit import Chem
        from rdkit.Chem import rdFMCS, AllChem 
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')

        refmol_str = next(pybel.readfile(os.path.splitext(refFile)[1].strip('.'), refFile)).write('pdb')
        refmol = Chem.MolFromPDBBlock(refmol_str, removeHs=True)
        if 'ratioThreshold' in kwargs: ratioThreshold = kwargs.pop('ratioThreshold')
        else: ratioThreshold=0.20

        targetmol = Chem.MolFromSmiles(targetMolSMI)
        targetmolHs = Chem.AddHs(targetmol, addCoords=True)
        mcs = rdFMCS.FindMCS([refmol, targetmolHs], threshold=0.9, completeRingsOnly=True)

        if mcs.smartsString and len(mcs.smartsString)>0:
            patt = Chem.MolFromSmarts(mcs.smartsString, mergeHs=True)
        else: raise ValueError('MCS not fund!')
    
        refmolcore = Chem.AllChem.ReplaceSidechains(refmol,patt)

        if refmolcore:
            core=AllChem.DeleteSubstructs(refmolcore, Chem.MolFromSmiles('*'))
            core.UpdatePropertyCache()
        else: raise AttributeError('Reference mcs core not found.')
        

        GetFF=lambda x,confId=-1:AllChem.MMFFGetMoleculeForceField(x,AllChem.MMFFGetMoleculeProperties(x),confId=confId)
        AllChem.ConstrainedEmbed(targetmolHs,core, getForceField=GetFF, useTethers=True)
        matchratio = float(core.GetNumAtoms())/float(refmolcore.GetNumAtoms())
        tethered_atom_ids=targetmolHs.GetSubstructMatches(patt)[0]

        if tethered_atom_ids and matchratio>ratioThreshold:
            atoms = map(lambda x:x+1, list(tethered_atom_ids))
            atoms_string = ','.join(str(el) for el in atoms)
        else: raise AttributeError('ratioThreshold error.')

        targetmolHs.SetProp('TETHERED ATOMS',atoms_string)
        print('Ligand: {}; Tethered atoms: {}'.format(targetMolID, atoms_string))
        
        w=Chem.SDWriter('{}_tethered.sd'.format(targetMolID))
        w.write(targetmolHs)
        w.flush()
        print('Output:', '{}_tethered.sd'.format(targetMolID), '\n')
        return '{}_tethered.sd'.format(targetMolID)

    def MolsAlignment(**kwargs):
        """:parameter reference: (str) reference 3D molecule
            :parameter target_name: (list of str) list of names of ligands
            :parameter target_smile: (list of str) list of smiles strings of ligands"""
        refFile = kwargs.pop('reference')
        targetMolID = kwargs.pop('targets_names')
        targetMolSMI = kwargs.pop('targets_smiles')
        
        if isinstance(targetMolID, list) and isinstance(targetMolSMI, list): 
            if len(targetMolID) == len(targetMolSMI):pass
            else: raise AttributeError('"target_name" and "target_smiles" have different leghts.')
        else: raise AttributeError('"target_name" and "target_smiles" must be a list')

        tetheredFiles = [Tethering.MolAlignment(ref=refFile, target_name=targetMolID[idx],target_smiles=targetMolSMI[idx]) for idx in range(len(targetMolID))]
        return tetheredFiles
