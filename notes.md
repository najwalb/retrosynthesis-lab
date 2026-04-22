# Inpainting
(no active issues)

# TODO [DONE]
- get diffalign samples to display on wsgi (get rid of dummy mols)
- integrate DiffAlign with syntheseus and add evaluation code with rxn-insight
- add options to: 1) apply full evaluation pipeline, 2) display precursors with levels of filtering, 3) add additional info next to the precursors
- selection highlights: now highlight bonds between selected atoms too, so
  implicit-carbon fragments show up visually
- selection summary: shows real fragment strings like "CC, C=O" per component
  instead of "N/M atoms of FULL_SMILES"
- fix /api/inpaint IndexError when sample_for_condition returns already-collapsed X