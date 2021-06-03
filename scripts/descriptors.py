# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 16:38:19 2021

@author: Li-Cheng Xu
"""
from process_utils import process_desc,maxminscale,zscorescale
from rdkit import Chem
from mordred import Calculator, descriptors
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect

def get_mordred_map(defined_chemical_space,scaler='maxmin'):
    '''
    desc_map: {category1:{smiles11:desc11,
                          smiles12:desc12},
               category2:{smiles21:desc21,
                          smiles22:desc22}}
    '''
    scaler = scaler.lower()
    assert scaler == 'maxmin' or scaler == 'z-score', 'scaler only support maxmin and z-score currently.'
    calc = Calculator(descriptors, ignore_3D=True)
    desc_map = {}
    if scaler == 'maxmin':
        for tmp_key in defined_chemical_space:
            smiles = defined_chemical_space[tmp_key]
            mols = [Chem.MolFromSmiles(tmp_smi) for tmp_smi in smiles]
            df = calc.pandas(mols).dropna(axis=1)
            desc = maxminscale(process_desc(df.to_numpy()))
            smi_desc_map = {tmp_smi:tmp_desc for tmp_smi,tmp_desc in zip(smiles,desc)}
            desc_map[tmp_key] = smi_desc_map
    elif scaler == 'z-score':
        for tmp_key in defined_chemical_space:
            smiles = defined_chemical_space[tmp_key]
            mols = [Chem.MolFromSmiles(tmp_smi) for tmp_smi in smiles]
            df = calc.pandas(mols).dropna(axis=1)
            desc = zscorescale(process_desc(df.to_numpy()))
            smi_desc_map = {tmp_smi:tmp_desc for tmp_smi,tmp_desc in zip(smiles,desc)}
            desc_map[tmp_key] = smi_desc_map
    return desc_map

def get_mf_map(defined_chemical_space,scaler='maxmin',radius=4,nBits=2048,useChirality=True):
    '''
    desc_map: {category1:{smiles11:desc11,
                          smiles12:desc12},
               category2:{smiles21:desc21,
                          smiles22:desc22}}
    '''
    scaler = scaler.lower()
    assert scaler == 'maxmin' or scaler == 'z-score', 'scaler only support maxmin and z-score currently.'
    desc_map = {}
    if scaler == 'maxmin':
        for tmp_key in defined_chemical_space:
            smiles = defined_chemical_space[tmp_key]
            mols = [Chem.MolFromSmiles(tmp_smi) for tmp_smi in smiles]
            desc = np.array([[eval(tmp_item) for tmp_item in Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(mol,
                              radius=radius,nBits=nBits,useChirality=useChirality).ToBitString()] for mol in mols])

            desc = maxminscale(process_desc(desc))
            smi_desc_map = {tmp_smi:tmp_desc for tmp_smi,tmp_desc in zip(smiles,desc)}
            desc_map[tmp_key] = smi_desc_map
    elif scaler == 'z-score':
        for tmp_key in defined_chemical_space:
            smiles = defined_chemical_space[tmp_key]
            mols = [Chem.MolFromSmiles(tmp_smi) for tmp_smi in smiles]
            desc = np.array([[eval(tmp_item) for tmp_item in Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(mol,
                              radius=radius,nBits=nBits,useChirality=useChirality).ToBitString()] for mol in mols])

            desc = zscorescale(process_desc(desc))
            smi_desc_map = {tmp_smi:tmp_desc for tmp_smi,tmp_desc in zip(smiles,desc)}
            desc_map[tmp_key] = smi_desc_map
    return desc_map

def get_rdkit_map(defined_chemical_space,scaler='maxmin'):
    '''
    desc_map: {category1:{smiles11:desc11,
                          smiles12:desc12},
               category2:{smiles21:desc21,
                          smiles22:desc22}}
    '''
    scaler = scaler.lower()
    assert scaler == 'maxmin' or scaler == 'z-score', 'scaler only support maxmin and z-score currently.'
    descs = [desc_name[0] for desc_name in Descriptors._descList]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(descs)
    desc_map = {}
    if scaler == 'maxmin':
        for tmp_key in defined_chemical_space:
            smiles = defined_chemical_space[tmp_key]
            mols = [Chem.MolFromSmiles(tmp_smi) for tmp_smi in smiles]
            desc = np.array([calc.CalcDescriptors(mol) for mol in mols])
            desc = maxminscale(process_desc(desc))
            smi_desc_map = {tmp_smi:tmp_desc for tmp_smi,tmp_desc in zip(smiles,desc)}
            desc_map[tmp_key] = smi_desc_map
    elif scaler == 'z-score':
        for tmp_key in defined_chemical_space:
            smiles = defined_chemical_space[tmp_key]
            mols = [Chem.MolFromSmiles(tmp_smi) for tmp_smi in smiles]
            desc = np.array([calc.CalcDescriptors(mol) for mol in mols])
            desc = zscorescale(process_desc(desc))
            smi_desc_map = {tmp_smi:tmp_desc for tmp_smi,tmp_desc in zip(smiles,desc)}
            desc_map[tmp_key] = smi_desc_map
    return desc_map
