from collections import defaultdict
import pandas as pd
import sys, os, re

import rdkit
from rdkit import Chem, RDLogger 
from rdkit.Chem import rdChemReactions

RDLogger.DisableLog('rdApp.*')

from LocalTemplate.template_extractor import extract_from_reaction

def remove_reagents(rxn):
    reactant = rxn.split('>>')[0]
    product = rxn.split('>>')[1]
    rs = reactant.split('.')
    ps = product.split('.')
    remove_frags = []
    for frag in ps:
        if frag in rs:
            remove_frags.append(frag)
    for frag in remove_frags:
        ps.remove(frag)
        rs.remove(frag)   
    return  '.'.join(sorted(rs)) + '>>' + '.'.join(sorted(ps))

def get_reaction_template(rxn, _id = 0):
    rxn = remove_reagents(rxn)
    rxn = {'reactants': rxn.split('>>')[0], 'products': rxn.split('>>')[1], '_id': _id}
    result = extract_from_reaction(rxn)
    return rxn, result

def dearomatic(template):
    for s in ['[c;', '[o;', '[n;', '[s;', '[c@']:
        template = template.replace(s, s.upper())
    return template

def fix_arom(mol):
    for atom in mol.GetAtoms():
        if not (atom.IsInRingSize(5) or atom.IsInRingSize(6)):
            atom.SetIsAromatic(False)
    return mol

def destereo(template):
    return template.replace('@', '')
        
def clean_smarts(sma):
    mol = fix_arom(Chem.MolFromSmarts(sma.replace(';', '')))
    smi = Chem.MolToSmiles(mol)
    try:
        smi = demap(smi)
        return smi
    except:
        return sma
    
def demap(smi):
    mol = Chem.MolFromSmiles(smi)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)

def reduce_template(temp):
    before_transform_smi = destereo(clean_smarts(temp.split('>>')[0]))
    before_transform_sma = destereo(temp.split('>>')[0])
    after_transform = temp.split('>>')[1]
    smi_template = before_transform_smi + '>>' + after_transform
    sma_template = before_transform_sma + '>>' + after_transform
    return smi_template, sma_template
