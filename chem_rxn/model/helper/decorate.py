from flask import current_app

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor


def decorate_svg(outcomes):
    result = []
    rank = 1
    for smis in outcomes:
        if rank > current_app.output_count:
            break

        ms = [Chem.MolFromSmiles(smi) for smi in smis.split('.')]
        svg = Draw.MolsToGridImage(ms, useSVG=True).replace('svg:','')
        result.append({
            'smiles': smis.split('.'),
            'rank': rank,
            'svg': svg,
        })
        rank += 1
    return result
