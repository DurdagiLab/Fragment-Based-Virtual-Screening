#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:24:37 2024

@author: Erenulug
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import BRICS
import types
import re
import logging


def read_smiles_from_excel(file_path):
    df = pd.read_excel(file_path)
    df = df[df['SMILES'].notna()]
    return [(row['Name'], row['SMILES']) for index, row in df.iterrows()] #If your smiles are missing after advanced search, it drops them for you.
#If your smiles are not missing, you can use these lines as comment lines. but not necessary


import logging

# Setting up logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def save_fragments_to_excel(fragments, output_file):
    df = pd.DataFrame(columns=['Name', 'Original SMILES', 'Fragment'])
    for fragment in fragments:
        name = fragment[0]
        original_smiles = fragment[1]
        fragments_smiles = fragment[2].split('.')
        for smile in fragments_smiles:
            # Remove patterns like "[1*]"
            cleaned_smile = re.sub(r'\[\d*\*\]', '', smile)
            if cleaned_smile:  # Only add non-empty fragments after cleaning
                df = df.append({'Name': name, 'Original SMILES': original_smiles, 'Fragment': cleaned_smile}, ignore_index=True)
    logging.info(f"Total fragments processed: {len(df)}")
    df.to_excel(output_file, index=False)
    logging.info(f"Output saved to {output_file}")

# Also, to make sure no errors are silently ignored, let's modify the `generate_fragments` function to log errors:

def generate_fragments(smiles_name_list):
    all_fragments = []
    for name, smiles in smiles_name_list:
        fragments_for_smiles = []
        if isinstance(smiles, float):
            logging.warning(f"Incorrect SMILES data: {smiles}")
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            try:
                frags = BRICS.BreakBRICSBonds(mol)
                if isinstance(frags, types.GeneratorType):
                    fragment_smiles = [(name, smiles, Chem.MolToSmiles(frag)) for frag in frags if frag is not None]
                else:
                    fragment_smiles = [(name, smiles, Chem.MolToSmiles(frags))] if frags is not None else []
            except Exception as e:
                logging.error(f"Error processing SMILES '{smiles}': {e}")
                continue
            fragments_for_smiles.extend(fragment_smiles)
        all_fragments.extend(fragments_for_smiles)
    return all_fragments


def save_fragments_to_excel(fragments, output_file):
    df = pd.DataFrame(columns=['Name', 'Original SMILES', 'Fragment'])
    for fragment in fragments:
        name = fragment[0]
        original_smiles = fragment[1]
        fragments_smiles = fragment[2].split('.')
        for smile in fragments_smiles:
            df = df.append({'Name': name, 'Original SMILES': original_smiles, 'Fragment': smile}, ignore_index=True)
    df.to_excel(output_file, index=False)

file_path = 'TNBC_advanced_metacore.xlsx' #Enter the xlsx file you have prepared or received from #metacore advanced search or search.
output_file = 'output_fragments.xlsx'

smiles_name_list = read_smiles_from_excel(file_path)
fragments = generate_fragments(smiles_name_list)
save_fragments_to_excel(fragments, output_file)

