"""Drill into the per_condition CSV: histogram of n_unique, and example products
in each bucket so we can eyeball what mode collapse looks like vs. the healthy middle."""
import pandas as pd
from pathlib import Path

OUT = Path(__file__).parent / 'uspto_analysis_out'
df = pd.read_csv(OUT / 'per_condition.csv')

print('-- n_unique histogram (binned) --')
bins = [0, 1, 2, 4, 10, 20, 30, 50, 80, 101]
labels = ['1', '2', '3-4', '5-10', '11-20', '21-30', '31-50', '51-80', '81-100']
hist = pd.cut(df['n_unique'], bins=bins, labels=labels, right=True).value_counts().reindex(labels, fill_value=0)
for lab, c in hist.items():
    bar = '#' * int(c / hist.max() * 50)
    print(f'  n_unique={lab:<7} {c:>5}  {bar}')

print('\n-- Truth-recovery vs. n_unique --')
for lab in ['1', '2', '3-4', '5-10', '11-20', '21-30', '31-50', '51-80', '81-100']:
    sub = df[pd.cut(df['n_unique'], bins=bins, labels=labels, right=True) == lab]
    if len(sub) == 0:
        continue
    print(f'  bucket {lab:<7} n={len(sub):>4}  truth_in_modes={sub["truth_in_modes"].mean()*100:5.1f}%  '
          f'truth_is_top={sub["truth_is_top"].mean()*100:5.1f}%  mean_tan={sub["mean_tanimoto"].mean():.3f}')

print('\n-- COLLAPSED group (n_unique == 1) --')
collapsed = df[df['n_unique'] == 1].head(8)
for _, r in collapsed.iterrows():
    truth_match = '✓ TRUTH' if r['truth_is_top'] else '✗ NOT truth'
    print(f"  cond {int(r['cond'])}: heavy={int(r['heavy_atoms']) if pd.notna(r['heavy_atoms']) else '?'} {truth_match}")
    print(f"    product: {r['true_product'][:90]}")
    print(f"    sole mode: {r['top_mode_smiles'][:90]}")
    print(f"    truth:     {r['true_reactants'][:90] if pd.notna(r['true_reactants']) else '(?)'}")

print('\n-- SWEET SPOT example (n_unique in [5,20]) --')
sweet = df[(df['n_unique'].between(5, 20)) & (df['truth_is_top'])].iloc[3]
print(f"  cond {int(sweet['cond'])}: n_unique={int(sweet['n_unique'])}, top_count={int(sweet['top_mode_count'])}, mean_tan={sweet['mean_tanimoto']:.2f}")
print(f"    product: {sweet['true_product'][:100]}")
print(f"    top mode (truth): {sweet['top_mode_smiles'][:100]}")

print('\n-- VERY HIGH diversity (n_unique >= 80) --')
chaos = df[df['n_unique'] >= 80].head(5)
for _, r in chaos.iterrows():
    print(f"  cond {int(r['cond'])}: n_unique={int(r['n_unique'])}, top_count={int(r['top_mode_count'])}, "
          f"mean_tan={r['mean_tanimoto']:.2f}, truth_in_modes={r['truth_in_modes']}, heavy={int(r['heavy_atoms']) if pd.notna(r['heavy_atoms']) else '?'}")

print('\n-- Top-mode count distribution --')
# How dominant is the top mode? top_count = 100 means full collapse.
for thr in [100, 90, 75, 50, 25, 10]:
    n = (df['top_mode_count'] >= thr).sum()
    print(f'  top_mode_count >= {thr:>3}: {n:>5}  ({n/len(df)*100:5.1f}%)')

print('\n-- Identity reactions per condition --')
print(f'  conditions with 0 identity     : {int((df["identity_count"] == 0).sum())} ({(df["identity_count"] == 0).mean()*100:.1f}%)')
print(f'  conditions with >=10 identity  : {int((df["identity_count"] >= 10).sum())}')
print(f'  conditions with >=50 identity  : {int((df["identity_count"] >= 50).sum())}')

print('\n-- Validity per condition --')
print(f'  mean valid    : {df["n_valid"].mean():.1f} / 100')
print(f'  conditions with <80 valid : {int((df["n_valid"] < 80).sum())}')
print(f'  conditions with <50 valid : {int((df["n_valid"] < 50).sum())}')
