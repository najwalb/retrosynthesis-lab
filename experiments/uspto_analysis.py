"""Mode/diversity analysis on the 5000-product USPTO-50k test set samples.

Reads the 10 sample shards in
  /home/laabidn1/laabidn1/DiffAlign/experiments/align_absorbing_20260407_144212/
each containing 500 conditions x 100 samples generated at 100 diffusion steps
from the deployed checkpoint (epoch760).

Outputs (next to this script, under experiments/uspto_analysis_out/):
  per_condition.csv      — one row per of the 5000 products
  mode_reuse.csv         — for each canonical reactant SMILES, count of distinct
                            products it served as the top mode for
  bucket_summary.txt     — quick text summary of the inpainting-relevance buckets
"""
import os, re, gzip, glob
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, RDConfig
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

ART_DIR = Path('/home/laabidn1/laabidn1/DiffAlign/experiments/align_absorbing_20260407_144212')
OUT_DIR = Path(__file__).parent / 'uspto_analysis_out'
OUT_DIR.mkdir(exist_ok=True)

SAMPLE_FILES = sorted(ART_DIR.glob('samples_epoch760_steps100_cond500_sampercond100_s*.txt'))
print(f'Found {len(SAMPLE_FILES)} sample files.')

# ── Parse ──────────────────────────────────────────────────────────────────
HEADER_RE = re.compile(r'^\(cond (\d+)\) (.+):\s*$')

def parse_shard(path):
    """Yield (cond_idx, true_rxn, [sample_rxn, ...])."""
    cond_idx = None
    true_rxn = None
    samples = []
    with open(path) as f:
        for raw in f:
            m = HEADER_RE.match(raw)
            if m:
                if cond_idx is not None:
                    yield cond_idx, true_rxn, samples
                cond_idx = int(m.group(1))
                true_rxn = m.group(2)
                samples = []
            elif raw.startswith('\t'):
                samples.append(raw.strip())
    if cond_idx is not None:
        yield cond_idx, true_rxn, samples

def split_rxn(rxn):
    """rxn -> (reactant_str, product_str). Returns (None, None) on malformed."""
    if '>>' not in rxn:
        return None, None
    rcts, prod = rxn.split('>>', 1)
    return rcts, prod

def canon(smi):
    """Sort fragments + canonicalise. Returns None if any fragment fails."""
    parts = smi.split('.')
    mols = [Chem.MolFromSmiles(p) for p in parts]
    if any(m is None for m in mols):
        return None
    return '.'.join(sorted(Chem.MolToSmiles(m, canonical=True) for m in mols))

def morgan_fp(smi, radius=2, nbits=2048):
    m = Chem.MolFromSmiles(smi)
    return AllChem.GetMorganFingerprintAsBitVect(m, radius, nbits) if m is not None else None

def pairwise_tanimoto(smiles_iter):
    fps = [morgan_fp(s) for s in smiles_iter]
    fps = [fp for fp in fps if fp is not None]
    if len(fps) < 2:
        return np.nan
    sims = []
    for i in range(len(fps)):
        for j in range(i + 1, len(fps)):
            sims.append(DataStructs.TanimotoSimilarity(fps[i], fps[j]))
    return float(np.mean(sims))

def product_complexity(prod_smi):
    m = Chem.MolFromSmiles(prod_smi)
    if m is None:
        return {'heavy_atoms': np.nan, 'rings': np.nan, 'stereocenters': np.nan}
    return {
        'heavy_atoms': m.GetNumHeavyAtoms(),
        'rings': m.GetRingInfo().NumRings(),
        'stereocenters': len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)),
    }

# ── Per-condition pass ─────────────────────────────────────────────────────
rows = []
top_modes = []                   # (canon_reactants_of_top_mode, product_canon)
total_processed = 0

for shard in SAMPLE_FILES:
    print(f'shard {shard.name} ...', flush=True)
    for cond_idx, true_rxn, samples in parse_shard(shard):
        true_rcts_raw, true_prod_raw = split_rxn(true_rxn)
        true_rcts_canon = canon(true_rcts_raw) if true_rcts_raw else None
        true_prod_canon = canon(true_prod_raw) if true_prod_raw else None

        # parse + canonicalise every sample
        sample_rcts_canon = []
        n_invalid = 0
        n_identity = 0
        for s in samples:
            rcts_raw, prod_raw = split_rxn(s)
            if rcts_raw is None:
                n_invalid += 1
                continue
            rcts_c = canon(rcts_raw)
            if rcts_c is None:
                n_invalid += 1
                continue
            sample_rcts_canon.append(rcts_c)
            # Identity rxn: reactants == product after canonicalisation
            prod_c = canon(prod_raw)
            if prod_c is not None and rcts_c == prod_c:
                n_identity += 1

        n_valid = len(sample_rcts_canon)
        counter = Counter(sample_rcts_canon)
        unique_modes = list(counter.keys())
        n_unique = len(unique_modes)
        top_mode = counter.most_common(1)[0] if counter else (None, 0)

        # diversity over unique modes (capped to avoid quadratic blow-up if huge)
        if n_unique > 1:
            mean_tan = pairwise_tanimoto(unique_modes[:50])
        else:
            mean_tan = np.nan

        # truth match
        truth_in_modes = bool(true_rcts_canon and true_rcts_canon in counter)
        truth_is_top = bool(true_rcts_canon and top_mode[0] == true_rcts_canon)
        truth_count = counter.get(true_rcts_canon, 0) if true_rcts_canon else 0

        cx = product_complexity(true_prod_canon) if true_prod_canon else {'heavy_atoms': np.nan, 'rings': np.nan, 'stereocenters': np.nan}

        # globally-unique condition id = shard_start + cond_idx
        shard_start = int(re.search(r'_s(\d+)\.txt$', shard.name).group(1))
        global_cond = shard_start + cond_idx

        rows.append({
            'cond':              global_cond,
            'true_product':      true_prod_canon,
            'true_reactants':    true_rcts_canon,
            'n_samples':         len(samples),
            'n_valid':           n_valid,
            'n_invalid':         n_invalid,
            'n_unique':          n_unique,
            'top_mode_count':    top_mode[1],
            'top_mode_smiles':   top_mode[0],
            'mean_tanimoto':     mean_tan,
            'identity_count':    n_identity,
            'truth_in_modes':    truth_in_modes,
            'truth_is_top':      truth_is_top,
            'truth_count':       truth_count,
            **cx,
        })

        if top_mode[0] is not None and true_prod_canon is not None:
            top_modes.append((top_mode[0], true_prod_canon))

        total_processed += 1

print(f'Processed {total_processed} conditions.')

df = pd.DataFrame(rows)
df.to_csv(OUT_DIR / 'per_condition.csv', index=False)
print(f'Wrote {OUT_DIR / "per_condition.csv"}')

# ── Mode reuse ─────────────────────────────────────────────────────────────
mode_to_products = {}
for top_smi, prod in top_modes:
    mode_to_products.setdefault(top_smi, set()).add(prod)
mode_reuse_rows = sorted(
    ({'top_mode_smiles': k, 'n_distinct_products': len(v)} for k, v in mode_to_products.items()),
    key=lambda r: r['n_distinct_products'], reverse=True,
)
mr_df = pd.DataFrame(mode_reuse_rows)
mr_df.to_csv(OUT_DIR / 'mode_reuse.csv', index=False)
print(f'Wrote {OUT_DIR / "mode_reuse.csv"}  (top mode reuse count = {mr_df["n_distinct_products"].max()})')

# ── Bucket summary ─────────────────────────────────────────────────────────
def bucket(n):
    if n <= 1:
        return '1 (collapsed)'
    if n <= 4:
        return '2-4 (low)'
    if n <= 20:
        return '5-20 (sweet spot for inpainting)'
    if n <= 50:
        return '21-50 (high)'
    return '51+ (very high)'

df['bucket'] = df['n_unique'].apply(bucket)
bucket_counts = df['bucket'].value_counts().reindex(
    ['1 (collapsed)', '2-4 (low)', '5-20 (sweet spot for inpainting)',
     '21-50 (high)', '51+ (very high)'], fill_value=0,
)

summary_lines = []
summary_lines.append(f'== USPTO-50k test (5000 products × 100 samples, 100 diffusion steps) ==')
summary_lines.append('')
summary_lines.append('-- Diversity / mode-collapse buckets (n_unique per condition) --')
for k, v in bucket_counts.items():
    summary_lines.append(f'  {k:<40} {v:>5}  ({v/len(df)*100:5.1f}%)')
summary_lines.append('')
summary_lines.append('-- Truth recovery --')
summary_lines.append(f'  truth in any mode      : {int(df["truth_in_modes"].sum()):>5}  ({df["truth_in_modes"].mean()*100:5.1f}%)')
summary_lines.append(f'  truth is the top mode  : {int(df["truth_is_top"].sum()):>5}  ({df["truth_is_top"].mean()*100:5.1f}%)')
summary_lines.append(f'  mean truth_count / 100 : {df["truth_count"].mean():>5.1f}')
summary_lines.append('')
summary_lines.append('-- Identity (reactants == product) --')
summary_lines.append(f'  mean identity_count / 100              : {df["identity_count"].mean():>5.1f}')
summary_lines.append(f'  conditions with >50 identity samples   : {int((df["identity_count"] > 50).sum())}')
summary_lines.append('')
summary_lines.append('-- Sample validity --')
summary_lines.append(f'  mean valid_count / 100   : {df["n_valid"].mean():>5.1f}')
summary_lines.append(f'  mean invalid_count / 100 : {df["n_invalid"].mean():>5.1f}')
summary_lines.append('')
summary_lines.append('-- Mode reuse (top mode appearing across distinct products) --')
summary_lines.append(f'  modes used as top for >=2 products : {int((mr_df["n_distinct_products"] >= 2).sum())}')
summary_lines.append(f'  modes used as top for >=10 products: {int((mr_df["n_distinct_products"] >= 10).sum())}')
summary_lines.append(f'  most-reused top mode               : {mr_df.iloc[0]["n_distinct_products"]} products')
summary_lines.append(f'    SMILES: {mr_df.iloc[0]["top_mode_smiles"][:120]}')
summary_lines.append('')
summary_lines.append('-- Diversity correlations --')
for col in ['heavy_atoms', 'rings', 'stereocenters']:
    sub = df[[col, 'n_unique']].dropna()
    if len(sub) > 1 and sub[col].std() > 0:
        corr = sub[col].corr(sub['n_unique'])
        summary_lines.append(f'  corr(n_unique, {col:<14}) = {corr:+.3f}  (n={len(sub)})')

text = '\n'.join(summary_lines)
print(text)
(OUT_DIR / 'bucket_summary.txt').write_text(text)
print(f'\nWrote {OUT_DIR / "bucket_summary.txt"}')
