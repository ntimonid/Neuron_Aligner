import pytest
from src import config

def test_config_port():
    assert config.PORT == 8000

def test_config_my_region():
    assert config.MY_REGION == 'VPM'

def test_config_soma_thr():
    assert config.SOMA_THR == 30

def test_config_cpd_params():
    assert 'max_it' in config.CPD_PARAMS
    assert config.CPD_PARAMS['max_it'] == 3
