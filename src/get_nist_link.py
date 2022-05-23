import argparse
import pickle
import jcamp
import numpy as np
import pubchempy as pcp
import re
import requests
from bs4 import BeautifulSoup
from rdkit import Chem

nist_id_rexp = re.compile(r'[0-9]+\-[0-9]+\-[0-9]+')

def process_nist_id(response, results, site_list, site_addr, inchikey):
    soup = BeautifulSoup(response.text, 'html.parser')
    prop_li = soup.findAll('li')
    for each_li in prop_li:
        text = each_li.get_text()
        if 'IUPAC Standard InChIKey:' in text:
            if not inchikey in text:
                return
    for each_li in prop_li:
        text = each_li.get_text()
        if 'CAS Registry Number:' in text:
            results.append(nist_id_rexp.search(text).group())
            site_list.append(site_addr)


def get_nist_id_from_inchi(smi):
    results = []
    sites = []
    inchikey = Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(smi))
    nist_site = "https://webbook.nist.gov/cgi/cbook.cgi?InChI=%s&Units=SI&Mask=200" % inchikey
    response = requests.get(nist_site, stream=False)
    if response.ok:
        if "The following matching species were found for the submitted" in response.text:
            soup = BeautifulSoup(response.text, 'html.parser')
            prop_a = soup.findAll(['a','h2'])
            for each_a in prop_a:
                if each_a.name == 'h2':
                    if "Similar" in each_a.string:
                        return results, sites
                    continue
                text = each_a.get("href")
                if "/cgi/cbook.cgi?" in text:
                    nist_site = "https://webbook.nist.gov" + text
                    response = requests.get(nist_site, stream=False)
                    if response.ok:
                        process_nist_id(response, results, sites, nist_site, inchikey)
        else:
            process_nist_id(response, results, sites, nist_site, inchikey)
    else:
        print("Cannot retrieve information with inchikey=%s" % inchikey)
    return results, sites

def get_mass_spectrum(nid):
    s = nid.replace('-','')
    nist_site ="https://webbook.nist.gov/cgi/cbook.cgi?JCAMP=C%s&Index=0&Type=Mass" % s
    response = requests.get(nist_site, stream=False)
    if response.ok:
        if "Spectrum not found" in response.text:
            return {}
        content = response.text.splitlines()
        results = jcamp.jcamp_read(content)
        return results
    return {}

def get_all_mass_spectrums(smi): 
    nid, sites = get_nist_id_from_inchi(smi)
    vals = []
    if not nid:
        return vals
    for c, s in zip(nid, sites):
        results = get_mass_spectrum(c)
        if results: 
            vals.append((s, results))
    return vals

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get links to flavors in the NIST Chemistry Webbook')
    parser.add_argument('-f', '--flavor', required=True,type=str)
    args = parser.parse_args()
    lno = 0
    all_info = []
    for line in open(args.flavor, "r"):
        line = line.strip()
        lno += 1
        val = get_all_mass_spectrums(line)
        sub_id = 0
        vl = len(val)
        if vl == 0:
            print("{0}\t{1}\tNone".format(lno, line, site))            
        if vl > 1:
            sub_id += 1
        for site, ei in val:
            if sub_id == 0:
                print("{0}\t{1}\t{2}".format(lno, line, site))
            else:
                print("{0}-{1}\t{2}\t{3}".format(lno, sub_id, line, site))
                sub_id += 1
