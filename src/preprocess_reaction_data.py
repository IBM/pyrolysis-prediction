# -*- coding:utf-8 -*-
import argparse
import csv
from rdkit import Chem
import preprocess_ea_data as ped

def process_reaction(line, ofile, flavor_db):
    row = line.split(" ")
    smiles = row[0].split(">")
    flag = False
    for s in smiles:
        mols = []
        try:
            ped.get_mols_without_h(s, mols)
        except:
            flag = True
            break
        for m in mols:
            if m in flavor_db:
                flag = True
                break
        if flag:
            break
    if not flag:
        ofile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process Chemical Reaction Data')
    parser.add_argument('-i', '--input', required=True,type=str)
    parser.add_argument('-f', '--flavor', required=True,type=str)
    parser.add_argument('-o', '--output', required=True,type=str)
    args = parser.parse_args()
    print(args)
    db = ped.generate_flaver_db(args.flavor)
    ofile = open(args.output, mode="w")
    for line in open(args.input):
        process_reaction(line, ofile, db)
            
        
            
        
        
        
