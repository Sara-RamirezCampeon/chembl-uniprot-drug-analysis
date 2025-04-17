
from chembl_webresource_client.new_client import new_client
import requests
import pandas as pd
from collections import defaultdict
import time

# Step 1: Retrieve all approved drugs from ChEMBL
drug_query = new_client.molecule.filter(max_phase=4).filter(atc_class__isnull=False)  # max_phase=4 filters approved drugs
approved_drugs = []

for drug in drug_query:
    if drug.get('molecule_chembl_id') and drug.get('first_approval'):
        approved_drugs.append({
            'name': drug.get('pref_name'),
            'chembl_id': drug.get('molecule_chembl_id'),
            'approval_year': drug.get('first_approval')
        })

# Convert to DataFrame and sort
df_drugs = pd.DataFrame(approved_drugs)
df_drugs = df_drugs.dropna(subset=['name']).sort_values(by=['approval_year', 'name'])

# Step 2: Get targets for approved drugs since 2019
recent_drugs = df_drugs[df_drugs['approval_year'] >= 2019]
print(f"Found {len(recent_drugs)} drugs approved since 2019.")

chembl_to_uniprot = defaultdict(list)

for _, row in recent_drugs.iterrows():
    chembl_id = row['chembl_id']
    target_query = new_client.target.filter(target_type="SINGLE PROTEIN", molecule_chembl_id=chembl_id)
    for target in target_query:
        for comp in target.get('target_components', []):
            for xref in comp.get('target_component_xrefs', []):
                if xref['xref_src_db'] == 'UniProt':
                    chembl_to_uniprot[chembl_id].append(xref['xref_id'])

# Step 3: Retrieve UniProt keywords for each UniProt accession
uniprot_accessions = set([acc for lst in chembl_to_uniprot.values() for acc in lst])
print(f"Fetching keywords for {len(uniprot_accessions)} UniProt entries...")

uniprot_keywords = {}

for accession in uniprot_accessions:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        keywords = [kw['id'] for kw in data.get('keywords', [])]
        uniprot_keywords[accession] = keywords
    else:
        print(f"Failed to fetch {accession} (status {response.status_code})")
    time.sleep(0.2)  # to avoid hitting API rate limits

# Optional: Print results for inspection
for chembl_id, accessions in chembl_to_uniprot.items():
    print(f"\nDrug: {chembl_id}")
    for acc in accessions:
        print(f"  UniProt: {acc} - Keywords: {uniprot_keywords.get(acc, [])}")
