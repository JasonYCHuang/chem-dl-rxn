from flask import current_app

from chem_rxn.model.helper.topk import translate_top_k
from chem_rxn.model.helper.decorate import decorate_svg


class ForwardModel:
    @classmethod
    def infer(cls, reactants):
        ans = translate_top_k(
            current_app.forward_model, reactants, k=5
        )
        return [decorate_svg(a) for a in ans]
