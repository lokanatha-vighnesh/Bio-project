from Bio import Entrez, SeqIO
import pandas as pd
from io import StringIO
import time
import certifi
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

Entrez.email = "lokanathavighnesh@gmail.com"

# ---------------------------------------
# 1️⃣ SEARCH FOR 1000 mRNA SEQUENCES
# ---------------------------------------

search_handle = Entrez.esearch(
    db="nucleotide",
    term="EGFR AND Homo sapiens AND lung cancer",
    retmax=1000
)

search_results = Entrez.read(search_handle)
search_handle.close()

id_list = search_results["IdList"]

print("Total IDs retrieved:", len(id_list))

# ---------------------------------------
# 2️⃣ FETCH SEQUENCES IN BATCHES
# ---------------------------------------

sequence_data = []

batch_size = 100

for start in range(0, len(id_list), batch_size):
    end = min(start + batch_size, len(id_list))
    batch_ids = id_list[start:end]

    print(f"Fetching records {start} to {end}")

    fetch_handle = Entrez.efetch(
        db="nucleotide",
        id=",".join(batch_ids),
        rettype="fasta",
        retmode="text"
    )

    fasta_data = fetch_handle.read()
    fetch_handle.close()

    records = SeqIO.parse(StringIO(fasta_data), "fasta")

    for record in records:
        sequence_data.append({
            "Accession_ID": record.id,
            "Description": record.description,
            "Sequence": str(record.seq)
        })

    time.sleep(1)  # NCBI rate limit safety

# ---------------------------------------
# 3️⃣ SAVE TO EXCEL
# ---------------------------------------

df = pd.DataFrame(sequence_data)

df.to_excel("lung_cancer_1000_sequences.xlsx", index=False)

print("Excel file created: lung_cancer_1000_sequences.xlsx")
