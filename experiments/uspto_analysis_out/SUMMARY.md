# USPTO-50k mode/diversity analysis

**Source:** 4,951 conditions × 100 samples each from
`align_absorbing_20260407_144212/samples_epoch760_steps100_cond500_sampercond100_s*.txt`
(epoch 760 checkpoint, **100 diffusion steps**).

## Headline

The deployed checkpoint is **not mode-collapsed** at 100 steps. The "low diversity" experience users see in the app comes from `diffusion_steps_eval=1`, not from the checkpoint. With proper sampling, the model produces healthy mode distributions and is genuinely product-conditional.

## Diversity buckets

| n_unique per condition | count | % | inpainting interest |
|---|---:|---:|---|
| 1 (collapsed) | 49 | 1.0% | none — but 92% of these *are* the ground truth (confidence, not failure) |
| 2-4 (low) | 559 | 11.3% | thin |
| **5-20 (sweet spot)** | **2,848** | **57.5%** | **inpainting actually offers alternatives** |
| 21-50 (high) | 1,359 | 27.4% | inpainting useful as a narrowing tool |
| 51+ (very high) | 136 | 2.7% | model is guessing — diversity is noise |

Histogram peak: **n_unique = 11–20** (1,450 conditions). Median product behaves like "model has ~10–15 plausible disconnections in mind, with one preferred."

## Truth recovery is excellent

- **truth in any mode: 85.6%** (4,237 / 4,951)
- **truth is the top mode: 56.9%** (2,819 / 4,951)
- mean truth_count / 100 = **41.7** (when truth is recovered, it dominates)

Inversely correlated with diversity bucket — collapsed cases are 92% truth-correct; very-high-diversity cases drop to 67% truth-in-modes / 24% truth-as-top. **Diversity is appropriately calibrated to model uncertainty**, not random noise.

## Top-mode dominance

| top_mode_count ≥ | conditions | % |
|---:|---:|---:|
| 100 (full collapse) | 28 | 0.6% |
| 90 | 709 | 14.3% |
| 75 | 1,592 | 32.2% |
| 50 (plurality) | 2,937 | 59.3% |
| 25 | 4,327 | 87.4% |
| 10 | 4,866 | 98.3% |

Even when "diverse," the top mode is rarely below 10% mass — there's structure, not chaos.

## Healthy signals

- **Mode reuse = 0.** No top-mode SMILES appears as the top mode for two distinct products. The model genuinely conditions on input — no generic-fallback shortcuts.
- **Identity reactions (reactants == product) = 0.8 / 100 mean.** 75.9% of conditions have **zero** identity samples. The "model just copies the product" failure is essentially absent.
- **Validity = 93.3 / 100 mean.** Only 13 conditions (0.3%) have <50% valid samples.
- **Diversity correlations weak** (`r ≈ 0.06–0.14` with heavy_atoms / rings / stereocenters) — diversity isn't pathologically tied to product size or complexity.

## What this means for the small-panel notebook results

In our 6-product CPU sweep, aspirin/menthol/ibuprofen collapsed to 1–4 modes even at 50 steps. The USPTO-50k results show:

- These products are **smaller and simpler than the typical USPTO-50k training example**. The model's behaviour on them is "I know one good disconnection, here it is."
- That's the same pattern as the **49 collapsed USPTO conditions** — and 92% of those also picked the *correct* truth as their sole mode. So this isn't a model-quality bug; it's the model being confident on easy/well-known transformations.
- **Menthol** (single mode across all our cells) is the more concerning failure — but it's also the only stereo-rich molecule we tested, and stereocenters has the strongest (still weak) positive correlation with diversity. Worth a separate stereo-specific investigation.

## Should you retrain on a harder dataset?

**Not for the diversity issue.** The data doesn't support it:
- 56.9% top-1 truth recovery on USPTO-50k test — solid retrosynthesis number.
- Mode distributions are healthy at 100 steps.
- The deployed app's diversity complaint is fully explained by `diffusion_steps_eval=1`.

**Maybe, for other reasons:**
- If you care about **small-molecule** chemistry beyond USPTO-50k's pharma-like distribution, a broader corpus (USPTO-MIT, USPTO-full) would expose the model to more variety — but small-molecule "modes" might still be 1 because the chemistry is genuinely deterministic.
- If **stereochemistry** is a priority, the menthol failure suggests the training data is stereo-poor or the loss treats stereo as low-weight. A stereo-aware dataset/loss change is a more targeted fix than a wholesale retrain.

## Recommended app changes (in priority order)

1. **Raise `diffusion_steps_eval` default from 1 to 25–50** in `DiffAlign/diffalign/inference.py:65` and the HTML input default in `wsgi.py:539`. Wall-clock is ~15–35s for 20 samples on CPU (well under the 600s gunicorn timeout). Highest-ROI change by far.
2. Keep the inpainting constraint check we just shipped — at 25–50 steps the "atoms marked to change but stuck" case will be rarer, but the empty-state UX is still important for the cases where it does happen.
3. Don't ship temperature scaling as a UI knob. The sweep showed no clean win and the convention is inverted from standard usage (confusing).

## Files

- `per_condition.csv` — one row per condition, all metrics.
- `mode_reuse.csv` — top-mode SMILES + how many products it served as top for (max = 1).
- `bucket_summary.txt` — text version of the bucket counts.
