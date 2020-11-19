import io
import numpy as np
import json
from flask import (
    Blueprint, request, jsonify
)
from chem_rxn.controller.helper.settings import get_ip_white_list
from chem_rxn.model.forward_model import ForwardModel
from chem_rxn.model.retro_model import RetroModel


ctrl = Blueprint('api', __name__)


@ctrl.before_app_request
def filter_remote_ip():
    trusted_servers = get_ip_white_list()
    # if request.remote_addr not in trusted_servers:
    #     abort(403)


@ctrl.route('/forward', methods=['POST'])
def forward_prediction():
    payload = request.json['reactants']
    rsp = ForwardModel.infer(payload)
    return jsonify(rsp)


@ctrl.route('/retro', methods=['POST'])
def retro_prediction():
    payload = request.json['products']
    rsp = RetroModel.infer(payload)
    return jsonify(rsp)
