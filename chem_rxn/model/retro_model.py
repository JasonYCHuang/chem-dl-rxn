from flask import current_app

from chem_rxn.model.helper.topk import translate_top_k
from chem_rxn.model.helper.decorate import decorate_svg


class RetroModel:
    @classmethod
    def infer(cls, products):
        ans = translate_top_k(
            current_app.retro_model, products, k=current_app.output_count
        )
        return [decorate_svg(a) for a in ans]
