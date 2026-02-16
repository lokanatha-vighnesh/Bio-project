from Bio import Entrez
import pandas as pd
import ssl

ssl._create_default_https_context = ssl._create_unverified_context
Entrez.email = "lokanathavighnesh@gmail.com"

search = Entrez.esearch(
    db="protein",
    term="lung[All Fields] AND Homo sapiens[Organism] AND srcdb_refseq[PROP]",
    retmax=1
)

result = Entrez.read(search)
search.close()

protein_id = result["IdList"][0]

fetch = Entrez.efetch(
    db="protein",
    id=protein_id,
    rettype="fasta",
    retmode="text"
)

wild_type_sequence = fetch.read()

fetch.close()
print(wild_type_sequence)
df = pd.DataFrame({protein_id : [wild_type_sequence]})

df.to_excel("lung_normal_wildseq.xlsx",index=False)

