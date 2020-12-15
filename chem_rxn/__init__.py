import os

from flask import Flask

import chem_rxn.load_models as lm


OUTPUT_COUNT = 10


def create_app(test_config=None):
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
        IP_WHITE_LIST='',
    )

    if test_config is None:
        app.config.from_pyfile('config.py', silent=True)
    else:
        app.config.from_mapping(test_config)

    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    forward_model = lm.load_forward_model()
    app.forward_model = forward_model
    retro_model = lm.load_retro_model()
    app.retro_model = retro_model

    app.output_count = OUTPUT_COUNT

    @app.route('/ping')
    def ping():
        return 'pong'


    from chem_rxn.controller.api import ctrl
    app.register_blueprint(ctrl)


    return app
