"""Thin re-export of PubChem compound lookup helpers."""
from evaluation.pubchem_lookup import get_compound_profile, lookup_all_compounds

__all__ = ['get_compound_profile', 'lookup_all_compounds']
