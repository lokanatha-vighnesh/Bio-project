from Bio import Entrez, SeqIO
import pandas as pd
from io import StringIO
import time
import ssl

# Temporary SSL bypass (since you had issues)
ssl._create_default_https_context = ssl._create_unverified_context

Entrez.email = "lokanathavighnesh@gmail.com"

# ---------------------------------------
# 1️⃣ SEARCH PROTEIN DATABASE
# ---------------------------------------

search_handle = Entrez.esearch(
    db="protein",
    term="EGFR AND Homo sapiens NOT lung cancer",
    retmax=1000
)

search_results = Entrez.read(search_handle)
search_handle.close()

id_list = search_results["IdList"]

print("Total Protein IDs retrieved:", len(id_list))

# ---------------------------------------
# 2️⃣ FETCH PROTEIN FASTA SEQUENCES
# ---------------------------------------

protein_data = []
batch_size = 100

for start in range(0, len(id_list), batch_size):
    end = min(start + batch_size, len(id_list))
    batch_ids = id_list[start:end]

    print(f"Fetching proteins {start} to {end}")

    fetch_handle = Entrez.efetch(
        db="protein",
        id=",".join(batch_ids),
        rettype="fasta",
        retmode="text"
    )

    fasta_data = fetch_handle.read()
    fetch_handle.close()

    records = SeqIO.parse(StringIO(fasta_data), "fasta")

    for record in records:
        protein_data.append({
            "Accession_ID": record.id,
            "Description": record.description,
            "Protein_Sequence": str(record.seq)
        })

    time.sleep(1)

# ---------------------------------------
# 3️⃣ SAVE TO EXCEL
# ---------------------------------------

df = pd.DataFrame(protein_data)
df.to_excel("lung_protein_sequences.xlsx", index=False)

print("Protein Excel file created successfully!")
