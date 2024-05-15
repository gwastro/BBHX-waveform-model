import numpy as np
from pycbc.waveform import (
    get_fd_det_waveform,
    get_fd_det_waveform_sequence,
    get_waveform_filter_length_in_time,
)
import pytest


@pytest.fixture(params=["BBHX_PhenomD", "BBHX_PhenomHM"])
def approximant(request):
    return request.param


@pytest.fixture(params=["LISA", "SSB"])
def ref_frame(request):
    return request.param


@pytest.fixture(params=["1.5", "2.0"])
def tdi(request):
    return request.param


@pytest.fixture()
def params():
    params = {}
    params["approximant"] = "BBHX_PhenomHM"
    params["tdi"] = 1.5
    params["ref_frame"] = "LISA"
    params["ifos"] = ["LISA_A", "LISA_E", "LISA_T"]
    params["coa_phase"] = 0.0
    params["mass1"] = 1e6
    params["mass2"] = 8e5
    params["spin1z"] = 0.0
    params["spin2z"] = 0.0
    params["distance"] = 410
    params["inclination"] = np.pi / 2
    params["t_obs_start"] = 31536000
    params["delta_f"] = 1.0 / params["t_obs_start"]
    params["f_lower"] = 1e-4
    params["f_ref"] = 8e-4
    params["f_final"] = 0.1
    params["delta_t"] = 1 / 0.2
    params["t_offset"] = 9206958.120016199
    params["tc"] = 4799624.274911478
    params["eclipticlongitude"] = 0.5
    params["eclipticlatitude"] = 0.23
    params["polarization"] = 0.1
    return params


def test_get_fd_det_waveform(params, ref_frame, approximant, tdi):
    params["tdi"] = tdi
    params["ref_frame"] = ref_frame
    params["approximant"] = approximant
    wf = get_fd_det_waveform(**params)
    # Check all 3 ifos are returned
    assert len(wf) == 3


def test_get_fd_det_waveform_sequence(params, approximant, tdi):
    params["ifos"] = "LISA_A"
    params["tdi"] = tdi
    params["approximant"] = approximant
    freqs = np.array([1-4, 1e-3, 1e-2])
    wf = get_fd_det_waveform_sequence(sample_points=freqs, **params)
    # Check all 3 ifos are returned
    assert len(wf) == 1
    assert len(wf["LISA_A"]) == len(freqs)


@pytest.mark.parametrize("mode_array", [None, [(3, 3)], [(2, 2), (3, 3)]])
def test_phenomhm_mode_array(params, mode_array):
    params["approximant"] = "BBHX_PhenomHM"
    params["mode_array"] = mode_array
    wf = get_fd_det_waveform(**params)
    assert len(wf) == 3


def test_length_in_time(params, approximant):
    params["approximant"] = approximant
    # print(params)
    duration = get_waveform_filter_length_in_time(**params)
    assert duration > 0
