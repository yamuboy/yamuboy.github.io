"""
PERM-to-Compustat Matching with LinkTransformer
=================================================
ECON 5293 — Nishant Yamujala

Run this LOCALLY (not in Claude.ai) — it requires:
  pip install linktransformer pandas

LinkTransformer uses sentence-transformer embeddings + cosine similarity
for semantic matching, which handles cases like:
  "NXP SEMICONDUCTORS NV" <-> "NXP SEMICONDUCTORS USA INC"
  "CIRRUS LOGIC INTERNATIONAL" <-> "CIRRUS LOGIC INC"
much better than character-level fuzzy matching (rapidfuzz/fuzzywuzzy).

INPUTS:
  - final_h1b_compustat_panel.csv  (your existing panel)
  - perm_panel.csv                 (output from harmonize_perm.py)

OUTPUTS:
  - lt_perm_compustat_crosswalk.csv  (matched employer names with scores)
  - h1b_perm_lt_matched_panel.csv    (enriched panel with PERM + CEM strata)
"""

import pandas as pd
import numpy as np
import linktransformer as lt
import os

# ============================================================
# CONFIG — adjust paths to where your files are
# ============================================================
PANEL_PATH = "final_h1b_compustat_panel.csv"
PERM_PATH = "perm_panel.csv"
OUTPUT_DIR = "."

# LinkTransformer model — all-MiniLM-L6-v2 is fast and good for English
# For even better company name matching, try "dell-research-harvard/lt-wikidata-comp-en"
LT_MODEL = "sentence-transformers/all-MiniLM-L6-v2"

# Match threshold — lower = more permissive, higher = stricter
# 0.5 is a reasonable starting point; review matches below this manually
SCORE_THRESHOLD = 0.5


# ============================================================
# STEP 1: Load data
# ============================================================
print("Loading data...")
panel = pd.read_csv(PANEL_PATH, low_memory=False)
perm = pd.read_csv(PERM_PATH, low_memory=False)

print(f"  Panel: {len(panel):,} rows, {panel['gvkey'].nunique()} firms")
print(f"  PERM:  {len(perm):,} rows")


# ============================================================
# STEP 2: Prepare unique employer name DataFrames
# ============================================================
print("\nPreparing employer name lists...")

# Compustat unique employers
compustat_names = (panel[['conm']].drop_duplicates()
                   .rename(columns={'conm': 'employer_name'})
                   .assign(employer_name=lambda x: x['employer_name'].str.upper().str.strip()))

# PERM unique employers
perm_names = (perm[['employer_name']].drop_duplicates()
              .assign(employer_name=lambda x: x['employer_name'].str.upper().str.strip()))

print(f"  Compustat unique names: {len(compustat_names)}")
print(f"  PERM unique names:     {len(perm_names):,}")


# ============================================================
# STEP 3: LinkTransformer merge
# ============================================================
print(f"\nRunning LinkTransformer merge (model: {LT_MODEL})...")
print("  This may take a few minutes on first run (downloads model)...")

# lt.merge does a 1-nearest-neighbor semantic merge
# It embeds all names from both sides, then finds the closest PERM name
# for each Compustat name using cosine similarity
matched = lt.merge(
    compustat_names,
    perm_names,
    merge_type='1:1',          # each Compustat firm gets 1 best match
    on="employer_name",         # column to match on
    model=LT_MODEL,
    suffixes=('_compustat', '_perm'),
)

# Rename score column
if 'score' not in matched.columns:
    # Some versions use different column names
    score_cols = [c for c in matched.columns if 'score' in c.lower() or 'sim' in c.lower()]
    if score_cols:
        matched = matched.rename(columns={score_cols[0]: 'score'})

print(f"\n  Matched {len(matched)} Compustat firms")
print(f"  Score distribution:")
print(f"    Mean:   {matched['score'].mean():.3f}")
print(f"    Median: {matched['score'].median():.3f}")
print(f"    Min:    {matched['score'].min():.3f}")
print(f"    > 0.9:  {(matched['score'] > 0.9).sum()}")
print(f"    > 0.7:  {(matched['score'] > 0.7).sum()}")
print(f"    > 0.5:  {(matched['score'] > 0.5).sum()}")


# ============================================================
# STEP 4: Review and filter matches
# ============================================================
print(f"\nHigh-confidence matches (score > 0.9):")
high = matched[matched['score'] > 0.9].sort_values('score', ascending=False)
for _, row in high.head(20).iterrows():
    print(f"  {row['employer_name_compustat'][:40]:<40} -> "
          f"{row['employer_name_perm'][:40]:<40} ({row['score']:.3f})")

print(f"\nBorderline matches (score {SCORE_THRESHOLD}-0.7) — REVIEW THESE:")
borderline = matched[matched['score'].between(SCORE_THRESHOLD, 0.7)].sort_values('score')
for _, row in borderline.iterrows():
    print(f"  {row['employer_name_compustat'][:40]:<40} -> "
          f"{row['employer_name_perm'][:40]:<40} ({row['score']:.3f})")

print(f"\nLow-confidence matches (score < {SCORE_THRESHOLD}) — LIKELY FALSE:")
low = matched[matched['score'] < SCORE_THRESHOLD].sort_values('score')
for _, row in low.head(10).iterrows():
    print(f"  {row['employer_name_compustat'][:40]:<40} -> "
          f"{row['employer_name_perm'][:40]:<40} ({row['score']:.3f})")

# Filter to threshold
crosswalk = matched[matched['score'] >= SCORE_THRESHOLD].copy()
crosswalk['match_quality'] = pd.cut(
    crosswalk['score'],
    bins=[0, 0.7, 0.9, 1.01],
    labels=['review', 'good', 'excellent']
)

print(f"\nKept {len(crosswalk)} matches above threshold {SCORE_THRESHOLD}")

# Save crosswalk
xwalk_path = os.path.join(OUTPUT_DIR, "lt_perm_compustat_crosswalk.csv")
crosswalk.to_csv(xwalk_path, index=False)
print(f"Saved crosswalk: {xwalk_path}")


# ============================================================
# STEP 5: Collapse PERM to employer-year and merge
# ============================================================
print("\nCollapsing PERM to employer-year...")

perm['employer_clean'] = perm['employer_name'].str.upper().str.strip()

perm_agg = perm.groupby(['employer_clean', 'fiscal_year']).agg(
    perm_total_filings=('case_number', 'count'),
    perm_certified=('case_status', lambda x: (x == 'CERTIFIED').sum()),
    perm_certified_expired=('case_status', lambda x: (x == 'CERTIFIED_EXPIRED').sum()),
    perm_denied=('case_status', lambda x: (x == 'DENIED').sum()),
    perm_from_h1b=('from_h1b', 'sum'),
    perm_mean_wage=('wage_offer_annual', 'mean'),
    perm_median_wage=('wage_offer_annual', 'median'),
).reset_index()

# Map Compustat names -> PERM names via crosswalk
name_map = dict(zip(
    crosswalk['employer_name_compustat'],
    crosswalk['employer_name_perm']
))

panel['employer_clean_comp'] = panel['conm'].str.upper().str.strip()
panel['perm_employer'] = panel['employer_clean_comp'].map(name_map)

# Merge
panel = panel.merge(
    perm_agg.rename(columns={'employer_clean': 'perm_employer',
                              'fiscal_year': 'fyear'}),
    on=['perm_employer', 'fyear'],
    how='left'
)

# Fill NaN with 0
for col in ['perm_total_filings', 'perm_certified', 'perm_certified_expired',
            'perm_denied', 'perm_from_h1b']:
    panel[col] = panel[col].fillna(0)

# Derived variables
panel['perm_intensity'] = np.where(
    panel['emp'] > 0, panel['perm_total_filings'] / (panel['emp'] * 1000), 0)

pre_perm = panel[(panel['fyear'] <= 2016) & (panel['perm_total_filings'] > 0)]['gvkey'].unique()
panel['perm_user_pre'] = panel['gvkey'].isin(pre_perm).astype(int)
panel['h1b_perm_pipeline'] = ((panel['h1b_user_pre'] == 1) &
                               (panel['perm_user_pre'] == 1)).astype(int)


# ============================================================
# STEP 6: Coarsened Exact Matching (NAICS3 × mktcap quartile)
# ============================================================
print("\nBuilding CEM strata...")

pre = panel[panel['fyear'].between(2014, 2016)]
firm_chars = pre.groupby('gvkey').agg(
    naics3=('naicsh_str', lambda x: str(x.iloc[0])[:3]),
    mean_mktcap=('mktcap', 'mean'),
    h1b_user_pre=('h1b_user_pre', 'first'),
).reset_index()

firm_chars['mktcap_qtile'] = firm_chars.groupby('naics3')['mean_mktcap'].transform(
    lambda x: pd.qcut(x, q=4, labels=[1,2,3,4], duplicates='drop')
)
overall_q = pd.qcut(firm_chars['mean_mktcap'], q=4, labels=[1,2,3,4], duplicates='drop')
firm_chars['mktcap_qtile'] = firm_chars['mktcap_qtile'].fillna(overall_q)

firm_chars['stratum'] = firm_chars['naics3'] + '_Q' + firm_chars['mktcap_qtile'].astype(str)

# Valid strata: has both T and C
strata_counts = firm_chars.groupby('stratum').agg(
    n_t=('h1b_user_pre', 'sum'),
    n_c=('h1b_user_pre', lambda x: (x==0).sum())
).reset_index()
valid = strata_counts[(strata_counts['n_t'] >= 1) & (strata_counts['n_c'] >= 1)]['stratum']

firm_chars['matched'] = firm_chars['stratum'].isin(valid).astype(int)
panel = panel.merge(firm_chars[['gvkey','stratum','naics3','mktcap_qtile','matched']],
                     on='gvkey', how='left')

n_matched = panel[panel['matched']==1]['gvkey'].nunique()
print(f"  Matched firms: {n_matched}")


# ============================================================
# STEP 7: Save
# ============================================================
out_path = os.path.join(OUTPUT_DIR, "h1b_perm_lt_matched_panel.csv")
panel.to_csv(out_path, index=False)
print(f"\nSaved enriched panel: {out_path}")
print(f"  {len(panel):,} rows, {len(panel.columns)} cols")
print(f"  PERM-matched firms: {panel['perm_employer'].notna().sum()} firm-years")
print(f"  H-1B + PERM pipeline: {panel[panel['h1b_perm_pipeline']==1]['gvkey'].nunique()} firms")
print(f"\nDone! Review lt_perm_compustat_crosswalk.csv for match quality.")
