import pandas as pd
import numpy as np
from collections import Counter
from itertools import product
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ==========================================
# CONFIGURATION
# ==========================================
FILE_PATH = "filtered_sequences.xlsx"
OUTPUT_FILE = "protein_features_output.xlsx"

# Standard amino acids
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
DIPEPTIDES = [''.join(p) for p in product(AMINO_ACIDS, repeat=2)]

# Amino acid property groups (for CTD)
HYDROPHOBIC = set("AILMFWV")
POLAR = set("STNQYC")
CHARGED = set("KRHDE")

# ==========================================
# LOAD DATA SAFELY
# ==========================================
df = pd.read_excel(FILE_PATH)

# Clean column names
df.columns = df.columns.str.strip()

# Auto-detect sequence column
sequence_columns = [col for col in df.columns 
                    if "seq" in col.lower() or "protein" in col.lower()]

if not sequence_columns:
    raise ValueError("No protein sequence column found in the Excel file.")

SEQ_COLUMN = sequence_columns[0]
print("Using sequence column:", SEQ_COLUMN)

# ==========================================
# FEATURE EXTRACTION FUNCTION
# ==========================================
def extract_protein_features(sequence):
    
    features = {}
    
    # Clean sequence
    seq = str(sequence).upper()
    seq = ''.join([aa for aa in seq if aa in AMINO_ACIDS])
    
    if len(seq) == 0:
        return None
    
    analysis = ProteinAnalysis(seq)
    
    # -----------------------------
    # Basic Physicochemical Features
    # -----------------------------
    features["Length"] = len(seq)
    features["Molecular_Weight"] = analysis.molecular_weight()
    features["Isoelectric_Point"] = analysis.isoelectric_point()
    features["Net_Charge_pH7"] = analysis.charge_at_pH(7.0)
    features["Aromaticity"] = analysis.aromaticity()
    features["Instability_Index"] = analysis.instability_index()
    features["GRAVY"] = analysis.gravy()
    
    # Aliphatic Index
    aa_percent = analysis.get_amino_acids_percent()
    features["Aliphatic_Index"] = (
        aa_percent['A'] +
        2.9 * aa_percent['V'] +
        3.9 * (aa_percent['I'] + aa_percent['L'])
    ) * 100
    
    # -----------------------------
    # Amino Acid Composition (20)
    # -----------------------------
    aa_counts = Counter(seq)
    for aa in AMINO_ACIDS:
        features[f"AAC_{aa}"] = aa_counts.get(aa, 0) / len(seq)
    
    # -----------------------------
    # Dipeptide Composition (400)
    # -----------------------------
    dipep_counts = Counter([seq[i:i+2] for i in range(len(seq)-1)])
    total_dipeptides = len(seq) - 1
    
    for dp in DIPEPTIDES:
        features[f"Dipep_{dp}"] = (
            dipep_counts.get(dp, 0) / total_dipeptides
            if total_dipeptides > 0 else 0
        )
    
    # -----------------------------
    # Hydrophobicity Statistics
    # -----------------------------
    hydropathy_scale = {
        'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,
        'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,'M':1.9,'N':-3.5,
        'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,
        'W':-0.9,'Y':-1.3
    }
    
    hydro_values = [hydropathy_scale[aa] for aa in seq]
    features["Hydrophobicity_Mean"] = np.mean(hydro_values)
    features["Hydrophobicity_Std"] = np.std(hydro_values)
    
    # -----------------------------
    # CTD Composition (Simplified)
    # -----------------------------
    features["CTD_Hydrophobic"] = sum(aa in HYDROPHOBIC for aa in seq) / len(seq)
    features["CTD_Polar"] = sum(aa in POLAR for aa in seq) / len(seq)
    features["CTD_Charged"] = sum(aa in CHARGED for aa in seq) / len(seq)
    
    return features

# ==========================================
# APPLY FEATURE EXTRACTION
# ==========================================
feature_list = []

for idx, row in df.iterrows():
    feats = extract_protein_features(row[SEQ_COLUMN])
    feature_list.append(feats if feats is not None else {})

feature_df = pd.DataFrame(feature_list)

# ==========================================
# PRESERVE SCORE FUNCTION (UNCHANGED)
# ==========================================
if "Score_Percentage" in df.columns:
    feature_df["Score_Percentage"] = df["Score_Percentage"]

# Merge original identifiers if present
final_df = pd.concat([df, feature_df], axis=1)

# ==========================================
# SAVE OUTPUT
# ==========================================
final_df.to_excel(OUTPUT_FILE, index=False)

print("Feature extraction completed successfully.")
print("Total features generated:", feature_df.shape[1])
print("Output file saved as:", OUTPUT_FILE)