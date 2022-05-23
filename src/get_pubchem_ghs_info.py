import argparse
import chardet
from rdkit import Chem
import pickle
import pubchempy as pcp
import urllib.request, json

def extract_ghs_from_pictograms(sec, hazard):
    for s in sec:
        if "Name" in s.keys():
            if s["Name"] == "Pictogram(s)":
                for items in s["Value"]["StringWithMarkup"]:
                    for v in items["Markup"]:
                        if "Extra" in v:
                            hazard.add(v["Extra"])
            break


def check_compound_ghs(smi):
    inchikey = Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(smi))
    compounds = pcp.get_compounds(inchikey, 'inchikey')
    pair = []
    if not compounds:
        hazard = set()
        hazard.add("unknown")
        r = [(hazard, "not found")]
        return r
    for c in compounds:
        hazard = set()
        base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/" + str(c.cid)
        readable_url = "https://pubchem.ncbi.nlm.nih.gov/compound/" + str(c.cid)
        url = base +  "/JSON"
        req = urllib.request.Request(url)
        res = urllib.request.urlopen(req).read()
        ### cp932
        success_flag = True
        try:
            sections = json.loads(res.decode('utf-8'))['Record']['Section']
        except:
            success_flag = False
        if not success_flag:
            success_flag = True
            try:
                sections = json.loads(res.decode('cp932'))['Record']['Section']
            except:
                success_flag = False
        if not success_flag:
            result = chardet.detect(res)
            charenc = result['encoding']
            sections = json.loads(res.decode(charenc))['Record']['Section']
        for sec in sections:
            if sec["TOCHeading"] == "Safety and Hazards":
                for subsec in sec["Section"]:
                    if subsec["TOCHeading"] == "Hazards Identification":
                        for subsubsec in subsec["Section"]:
                            if subsubsec["TOCHeading"] == "GHS Classification":
                                extract_ghs_from_pictograms(subsubsec["Information"], hazard)
                                break
                        break
                break
        if hazard:
            pair.append((hazard, readable_url))
        else:
            hazard.add("none")
            pair.append((hazard, readable_url))
    return pair


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate GHS status for each product extracted from pubchem')
    parser.add_argument('-p', '--products', required=True,type=str)
    args = parser.parse_args()
    print(args)
    for p in open(args.products):
        p = p.rstrip()
        results = check_compound_ghs(p)
        for pair in results:
            print("{0}\t{1}\t{2}".format(p, pair[0], pair[1]))

