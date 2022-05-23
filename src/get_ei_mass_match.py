import argparse
import csv
import pickle
import jcamp
import numpy as np
import re
import requests
from bs4 import BeautifulSoup
from rdkit import Chem
from rdkit.Chem import Descriptors

def download_ei_mass(website):
    response = requests.get(website, stream=False)
    if response.ok:
        soup = BeautifulSoup(response.text, 'html.parser')
        prop_a = soup.findAll(['p'])
        for each_a in prop_a:
            text = each_a.get_text()
            if "in JCAMP-DX" in text:
                web = "https://webbook.nist.gov" + each_a.a.get("href")
                response = requests.get(web, stream=False)
                if response.ok:
                    if "Spectrum not found" in response.text:
                        print("Spectrum not found with c= {0}".format(website))
                        return (None, None)
                content = response.text.splitlines()
                results = jcamp.jcamp_read(content)
                return results['x'], [v / 100.0 for v in results['y']]
    else:
        print("Website not found: {0}".format(website))
    return None

def get_ei_mass_value(p, eim):
    value = 0.0
    w = int(Descriptors.ExactMolWt(Chem.MolFromSmiles(p))) - 1
    for x, y in zip(eim[0], eim[1]):
        if x == w:
            value = max([value, y])
    return value


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract EI mass from NIST Chemistry Webbook')
    parser.add_argument('-w', '--web', required=True, type=str)
    parser.add_argument('-p', '--products', required=True, type=str)
    parser.add_argument('-t', '--threshold', required=True, type=float)
    args = parser.parse_args()
    eim = download_ei_mass(args.web)
    if eim is not None:
        print("EI mass matches")
        for smi in open(args.products):
            smi = smi.rstrip()
            value = get_ei_mass_value(smi, eim)
            if value >= args.threshold:
                print("{0}\t{1}".format(smi, value))
