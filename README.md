# pyrolysis-prediction
Sample code to generate data and extract information on pyolysis

## Required packages

You will need to install RDKit, pubchempy and jcamp using pip. 

## Script to remove test flavor compounds from the training data (pyrolysis prediction)

Usage: ```python ./preprocess_reaction_data.py -i TRAINING_DATA -f FLAVOR_DATA -o NEW_TRAINING_DATA```

Each line in the training data should be formatted as reaction SMILES which includes molecules only in the reactant and product parts. 
Each line in the flavor data includes a flavor compound in SMILES format, such as:

```
CCC(=O)C(=O)C
CCCCCO
COC1=CC=C(C=C1)OC
CC(C(CCC)=O)=O
CC(=O)C1=CC=CO1
...
```

## Script to remove test flavor compounds from the training data (pyrolysis prediction)

Usage: ```python ./preprocess_ea_data.py -i TRAINING_DATA -f FLAVOR_DATA -o NEW_TRAINING_DATA```

The training data should be in CSV format which is acceptable for the reaction branch of chemprop 
The flavor data has the same format described above.


## Getting links to the flavors in the NIST database

Usage: ```python ./get_nist_link.py -f FLAVOR_DATA```


The flavor data has the same format described above.

## Finding EI mass spectrum matches to the products

Usage: ```python -w "NIST_WEBSITE" -p PRODUCT_FILE -t THRESHOLD```

NIST_WEBSITE is a link to a website obtained by `get_nist_link.py`, (e.g., https://webbook.nist.gov/cgi/cbook.cgi?InChI=TZMFJUDUGYTVRY-UHFFFAOYSA-N&Units=SI&Mask=200 for 2,3-Pentanedione)

PRODUCT_FILE includes a list of products generated from a flavor. Each line of the file consists of a molecule in SMILES format. 

THRESHOLD is a threshold of the EI mass spectrum values (set to 5.0 in our paper). 

## Finding GHS Classification of molecules in the PubChem database

Usage: ```python ./get_pubchem_ghs_info.py -p PRODUCT_FILE```

PRODUCT_FILE includes a list of products generated from a flavor. Each line of the file consists of a molecule in SMILES format. 
