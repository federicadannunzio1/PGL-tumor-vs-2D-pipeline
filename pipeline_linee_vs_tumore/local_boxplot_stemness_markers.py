#!/usr/bin/env python3
# =============================================================================
# local_boxplot_stemness_markers.py
# Boxplot facettati (scala libera) per marker stemness + neuroendocrini
# Confronto tumore primario vs linee 2D
#
# Eseguire in locale (Mac):
#   python3 local_boxplot_stemness_markers.py
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# ---- Paths ----
PIPELINE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR  = os.path.join(PIPELINE_DIR, "results")
TPM_FILE     = os.path.join(RESULTS_DIR, "bulk_preprocessing", "tpm_matrix.csv")
META_FILE    = os.path.join(RESULTS_DIR, "bulk_preprocessing", "sample_metadata_ordered.csv")
DEGS_FILE    = os.path.join(RESULTS_DIR, "bulk_degs", "degs_full_table.csv")
OUTPUT_DIR   = os.path.join(RESULTS_DIR, "boxplot_markers")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---- Target genes ----
TARGET_GENES = [
    "SOX2", "PROM1", "ALDH1A1", "LGR5", "EPCAM", "CD44", "KLF4",
    "BMI1", "EZH2", "MYCN", "HMGA2", "IGF2BP2", "IGF2BP3",
    "CHGA", "CHGB", "INSM1", "SYP", "ASCL1", "NEUROD1",
    "DLL3", "CD276"  # DLL3 = DDL3 corretto; CD276 = B7-H3
]

COLORS = {"tumor": "#E64B35", "2D": "#4DBBD5"}

# ---- Load data ----
print("Caricamento dati...")
tpm = pd.read_csv(TPM_FILE, index_col=0)
meta = pd.read_csv(META_FILE)
degs = pd.read_csv(DEGS_FILE)

# Strip ENSEMBL version numbers
tpm.index = tpm.index.str.replace(r'\.\d+$', '', regex=True)

# Build ENSEMBL -> symbol mapping from DEG table
ens2sym = dict(zip(degs['gene_id'], degs['gene_symbol']))
sym2ens = {}
for ens, sym in ens2sym.items():
    if pd.notna(sym):
        sym2ens[sym] = ens

# ---- Find available genes ----
found_genes = []
missing_genes = []
for g in TARGET_GENES:
    if g in sym2ens and sym2ens[g] in tpm.index:
        found_genes.append(g)
    else:
        missing_genes.append(g)

print(f"Geni trovati: {len(found_genes)}/{len(TARGET_GENES)}")
if missing_genes:
    print(f"Geni mancanti: {', '.join(missing_genes)}")

# ---- Build long-format data ----
np.random.seed(42)
rows = []
for gene in found_genes:
    ens_id = sym2ens[gene]
    for _, row in meta.iterrows():
        sample = row['sample_id']
        cond = row['condition']
        tpm_val = tpm.loc[ens_id, sample]
        log2_val = np.log2(tpm_val + 1)
        rows.append({
            'gene': gene, 'condition': cond,
            'log2TPM': log2_val, 'sample': sample
        })

df = pd.DataFrame(rows)

# ---- Plot ----
n_genes = len(found_genes)
ncols = 4
nrows = int(np.ceil(n_genes / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
fig.suptitle(
    'Stemness & Neuroendocrine markers: Tumor vs 2D cell lines\n'
    'log\u2082(TPM+1), free Y-axis per gene',
    fontsize=13, fontweight='bold', y=1.01
)

axes_flat = axes.flatten() if n_genes > 1 else [axes]

for idx, gene in enumerate(found_genes):
    ax = axes_flat[idx]
    sub = df[df['gene'] == gene]

    for i, cond in enumerate(['tumor', '2D']):
        vals = sub[sub['condition'] == cond]['log2TPM'].values
        bp = ax.boxplot(
            vals, positions=[i], widths=0.4,
            patch_artist=True,
            medianprops=dict(color='black', linewidth=2),
            whiskerprops=dict(linewidth=1.2),
            capprops=dict(linewidth=1.2),
            flierprops=dict(marker='')
        )
        bp['boxes'][0].set_facecolor(COLORS[cond])
        bp['boxes'][0].set_alpha(0.35)
        x_jit = np.random.normal(i, 0.06, size=len(vals))
        ax.scatter(x_jit, vals, color=COLORS[cond], s=40, zorder=3,
                   alpha=0.9, edgecolors='white', linewidths=0.3)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Tumor', '2D'], fontsize=9)
    ax.set_ylabel('log\u2082(TPM+1)', fontsize=9)
    ax.set_title(gene, fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # y_min a 0 per chiarire che valori > 0 = gene espresso
    ax.set_ylim(bottom=-0.3)
    ax.axhline(y=0, color='grey', linestyle=':', linewidth=0.8, alpha=0.5)

# Hide empty panels
for idx in range(n_genes, len(axes_flat)):
    axes_flat[idx].set_visible(False)

plt.tight_layout()

out_pdf = os.path.join(OUTPUT_DIR, "boxplot_stemness_neuroendocrine_markers.pdf")
out_png = os.path.join(OUTPUT_DIR, "boxplot_stemness_neuroendocrine_markers.png")
plt.savefig(out_pdf, bbox_inches='tight')
plt.savefig(out_png, bbox_inches='tight', dpi=150)
print(f"\nSalvato:\n  {out_pdf}\n  {out_png}")

# ---- Summary table ----
summary_rows = []
for gene in found_genes:
    ens_id = sym2ens[gene]
    sub = df[df['gene'] == gene]
    mean_tumor = sub[sub['condition'] == 'tumor']['log2TPM'].mean()
    mean_2d    = sub[sub['condition'] == '2D']['log2TPM'].mean()

    # Get DESeq2 stats if available
    deg_row = degs[degs['gene_symbol'] == gene]
    lfc = deg_row['log2FoldChange'].values[0] if len(deg_row) > 0 else np.nan
    padj = deg_row['padj'].values[0] if len(deg_row) > 0 else np.nan

    summary_rows.append({
        'gene': gene,
        'mean_log2TPM_tumor': round(mean_tumor, 3),
        'mean_log2TPM_2D': round(mean_2d, 3),
        'DESeq2_log2FC': round(lfc, 3) if not np.isnan(lfc) else 'NA',
        'DESeq2_padj': f'{padj:.2e}' if not np.isnan(padj) else 'NA'
    })

summary_df = pd.DataFrame(summary_rows)
summary_csv = os.path.join(OUTPUT_DIR, "stemness_neuroendocrine_markers_summary.csv")
summary_df.to_csv(summary_csv, index=False)
print(f"  {summary_csv}")
