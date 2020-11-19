from __future__ import print_function
from fairseq.models.transformer import TransformerModel

from pathlib import Path



forward_root_dir = f'{Path().absolute()}/chem_rxn/nlp-rxn-main/forward'

def load_forward_model():
    forwardModel = TransformerModel.from_pretrained(
        str(Path().absolute()),
        checkpoint_file=f'{forward_root_dir}/checkpoints/checkpoint_best.pt',
        data_name_or_path=f'{forward_root_dir}/data-bin/uspto.rct-prd/',
        bpe='subword_nmt',
        bpe_codes=f'{forward_root_dir}/preprocess/uspto/code'
    )

    forwardModel.eval()

    print('Load FORWARD model OK!')
    return forwardModel



retro_root_dir = f'{Path().absolute()}/chem_rxn/nlp-rxn-main/retro'

def load_retro_model():
    retroModel = TransformerModel.from_pretrained(
        str(Path().absolute()),
        checkpoint_file=f'{retro_root_dir}/checkpoints/checkpoint_best.pt',
        data_name_or_path=f'{retro_root_dir}/data-bin/uspto.prd-rct/',
        bpe='subword_nmt',
        bpe_codes=f'{retro_root_dir}/preprocess/uspto/code'
    )

    retroModel.eval()

    print('Load RETRO model OK!')
    return retroModel
