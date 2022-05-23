# -*- coding:utf-8 -*-
import argparse
import csv
from rdkit import Chem

def get_mols_without_h(smiles, mols):
    for s in smiles.split("."):
        m = Chem.MolFromSmiles(s)
        for a in m.GetAtoms():
            a.SetAtomMapNum(0)
            a.SetNoImplicit(False)
        Chem.SanitizeMol(m)
        mols.append(Chem.MolToSmiles(m))

def generate_flaver_db(file_name):
    db = set()
    for l in open(file_name):
        row = l.rstrip()
        m = Chem.MolFromSmiles(row)
        db.add(Chem.MolToSmiles(m))
    return db
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Preprocess Energy Activation Data')
    parser.add_argument('-i', '--input', required=True,type=str)
    parser.add_argument('-f', '--flavor', required=True,type=str)
    parser.add_argument('-o', '--output', required=True,type=str)
    args = parser.parse_args()
    print(args)
    db = generate_flaver_db(args.flavor)
    csvreader = csv.reader(open(args.input))
    header = next(csvreader)
    print(header)
    ofile = csv.writer(open(args.output, "w"))
    header = ["rsmi","psmi", "ea", "dh"]
    ofile.writerow(header)
    for row in csvreader:
        mols = []
        get_mols_without_h(row[1], mols)
        get_mols_without_h(row[2], mols)
        flag = False
        for ms in mols:
            if ms in db:
                flag = True
                break
        if flag:
            continue
        ofile.writerow(row[1:])
