import numpy as np
from pycbc.waveform import get_fd_det_waveform
import pytest


@pytest.mark.parametrize("ref_frame", ["SSB", "LISA"])
@pytest.mark.parametrize("tdi", ["1.5", "2.5"])
def test_get_fd_det_waveform(ref_frame):
    params = {}
    params["tdi"] = tdi
    params["ref_frame"] = ref_frame
    params["approximant"] = "BBHX_PhenomD"
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
    wf = get_fd_det_waveform(ifos=["LISA_A", "LISA_E", "LISA_T"], **params)

    # Check all 3 ifos are returned
    assert len(wf) == 3
