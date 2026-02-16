from Bio import pairwise2
import pandas as pd

# -----------------------------------
# 1️⃣ LOAD FILES
# -----------------------------------

all_sequences_path = "combined_file.xlsx"
wildtype_path = "lung_cancer_wildseq.xlsx"

all_sequences_df = pd.read_excel(all_sequences_path)
wild_df = pd.read_excel(wildtype_path)

sequence_column = "Sequence"

all_sequences = all_sequences_df[sequence_column].dropna().astype(str)
wild_type_sequence = str(wild_df.iloc[0, 0])

print("Total sequences loaded:", len(all_sequences))

# -----------------------------------
# 2️⃣ GLOBAL ALIGNMENT WITH SCORING
# -----------------------------------

# Scoring scheme

def compute_score_percentage(seq1, seq2):
    alignment = pairwise2.align.globalxx(seq1, seq2)[0]
    
    alignment_score = alignment.score
    
    max_possible_score = len(seq1) * 2
    
    percent_score = (alignment_score / max_possible_score) * 100
    
    return percent_score

# -----------------------------------
# 3️⃣ FILTER USING 80% SCORE
# -----------------------------------

filtered_data = []

for idx, seq in enumerate(all_sequences):
    score_percent = compute_score_percentage(wild_type_sequence, seq)
    
    if score_percent <= 20:
        filtered_data.append({
            "Protein_Sequence": seq,
            "Score_Percentage": score_percent
        })
    
    if idx % 50 == 0:
        print(f"Processed {idx} sequences")

# -----------------------------------
# 4️⃣ SAVE OUTPUT
# -----------------------------------

filtered_df = pd.DataFrame(filtered_data)

output_path = "filtered_sequences.xlsx"
filtered_df.to_excel(output_path, index=False)

print("Filtering complete.")
print("Sequences retained:", len(filtered_df))
print("Saved to:", output_path)