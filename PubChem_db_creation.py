#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 20:52:45 2019

@author: alexorlov
"""

import sqlite3
import csv
from io import open

from rdkit import RDLogger
from rdkit.Chem import MolFromInchi, MolToSmiles, MolFromSmarts
from rdkit.Chem.Descriptors import MolWt, ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from tqdm import tqdm
import datetime
now = datetime.datetime.now()

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

###Setting the path to PubChem Database using CID-InChI-Key file available via PubChem ftp site ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz
pubchem='/Kostyukevich_project/Databases/PubChem/CID-InChI-Key' 

###Function for creation of sqlite version of PubChem Database
def create_database():
    conn = sqlite3.connect('/Kostyukevich_project/Databases/PubChem/PubChem_{}.db'.format(now.strftime("%Y_%m_%d")))
    c = conn.cursor()
    c.execute("CREATE TABLE IF NOT EXISTS molecules ("
              "id integer primary key, "
              "smiles BLOB, "
              "mass  REAL, "
              "exactmass REAL,"
              "formula TEXT, "
              "CID, "
              "key BLOB"
              ")")
    c.execute("CREATE INDEX massIndex ON molecules(exactmass)")
    c.execute("CREATE INDEX formulaIndex ON molecules(formula)")
    return c,conn

def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b


##Main body of the script. At first, connection to sqlite database (created on the previous stage) is setting up. Then, structures of molecules are reading from Inchi format into RDKit mol format using RDKit python library.  Molecular descriptors (mass, exactmass) and smiles notations are calculating through RDKit python library. The results of the calculations are putting into database.

with open(pubchem, "r", encoding="utf-8", errors='ignore') as f:
    nlines  = sum(bl.count("\n") for bl in blocks(f))
print(nlines)

c,conn = create_database()

with open(pubchem, newline='') as csvfile:
     spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
     idx = 0
     for row in tqdm(spamreader,total=nlines):
            cid = row[0]
            inchi = row[1]
            key = row[2]
            try:
                mol = MolFromInchi(inchi)
                mass = float(MolWt(mol))
                exactmass = float(ExactMolWt(mol))
                formula = CalcMolFormula(mol)
                if mass > 1500.: continue
                if len(mol.GetSubstructMatches(MolFromSmarts('[C]'))) < 2.: continue
                smiles = MolToSmiles(mol).encode('ascii')
                key = key.encode('ascii')
                c.execute('INSERT INTO molecules VALUES (?,?,?,?,?,?,?)', (idx,smiles,mass,exactmass,formula,cid,key))
                idx += 1
                if idx % 50000 ==0:
                    conn.commit()
            except:
                pass
     conn.commit()
