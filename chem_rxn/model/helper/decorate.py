from flask import current_app

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor


def decorate_svg(outcomes):
    results = []
    rank = 1

    for smis in outcomes:
        if rank > current_app.output_count:
            break

        smis = smis.replace(' ', '.')
        smis = [smi for smi in smis.split('.') if smi]
        ms = [Chem.MolFromSmiles(smi) for smi in smis]
        if None in ms: # filter invalid SMILES syntax
            continue

        svg = Draw.MolsToGridImage(ms, useSVG=True).replace('svg:','')
        results.append({
            'smiles': smis,
            'rank': rank,
            'svg': svg,
        })
        rank += 1

    return results
