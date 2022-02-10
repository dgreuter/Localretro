import os, sys, re
import pandas as pd
from argparse import ArgumentParser

import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions
RDLogger.DisableLog('rdApp.*')
sys.path.append('../')

def matchwithtemp(template1, template2, replacement_dict):
    template1 = template1.split('>>')[1]
    template2 = template2.split('>>')[1]
    matched_idx = {}
    for char1, char2 in zip(template1, template2):
        if char1 == char2:
            pass
        else:
            matched_idx[char1] = char2

    new_replacement_dict = {}
    for k in replacement_dict.keys():
        if replacement_dict[k] not in matched_idx.keys():
            new_replacement_dict[k] = replacement_dict[k]
        else:
            new_replacement_dict[k] = matched_idx[replacement_dict[k]]
    return new_replacement_dict

def match_num(edit_idx, replace_dict):
    new_edit_idx = ''
    idx = ''
    for s in str(edit_idx):
        if s.isdigit():
            idx += s
        else:
            if idx.isdigit():
                new_edit_idx += replace_dict[idx]
            new_edit_idx += s
            idx = ''
    return new_edit_idx

def get_idx_map(product_smiles, replace_dict):
    mol = Chem.MolFromSmiles(product_smiles)
    idx_map = {}
    for atom in mol.GetAtoms():
        idx = str(atom.GetIdx())
        atom_map = str(atom.GetAtomMapNum())
        idx_map[atom_map] = idx
    return {i:idx_map[k] for k, i in replace_dict.items()}

def get_edit_site(smiles):
    mol = Chem.MolFromSmiles(smiles)
    A = [a for a in range(mol.GetNumAtoms())]
    B = []
    for atom in mol.GetAtoms():
        others = []
        bonds = atom.GetBonds()
        for bond in bonds:
            atoms = [bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()]
            other = [a for a in atoms if a != atom.GetIdx()][0]
            others.append(other)
        b = [(atom.GetIdx(), other) for other in sorted(others)]
        B += b
    return A, B
