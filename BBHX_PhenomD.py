import functools
import math
import numpy as np
from scipy.interpolate import interp1d
from bbhx.utils.constants import MTSUN_SI, YRSID_SI
from bbhx.waveformbuild import BBHWaveformFD
from bbhx.utils.transform import LISA_to_SSB


@functools.lru_cache(maxsize=128)
def get_waveform_genner(log_mf_min, run_phenomd=True):
    # See below where this function is called for description of how we handle
    # log_mf_min.
    mf_min = math.exp(log_mf_min/25.)
    wave_gen = BBHWaveformFD(amp_phase_kwargs=dict(run_phenomd=run_phenomd, mf_min=mf_min))
    return wave_gen

@functools.lru_cache(maxsize=10)
def cached_arange(start, stop, spacing):
    return np.arange(start, stop, spacing)

def chirptime(m1, m2, f_lower):
    from pycbc.waveform.spa_tmplt import findchirp_chirptime

    duration = findchirp_chirptime(m1=m1, m2=m2, fLower=f_lower, porder=7)
    return duration

def imr_duration(m1, m2, s1z, s2z, f_lower):
    # More accurate duration of the waveform, including merge, ringdown,
    # and aligned spin effects.
    import lal
    from pycbc import libutils
    lalsimulation = libutils.import_optional('lalsimulation')
    time_length = lalsimulation.SimIMRPhenomDChirpTime(
                        m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, s1z, s2z, f_lower)
    return time_length * 1.1

def interpolated_tf(m1, m2):
    # Using findchirp_chirptime in PyCBC to calculate 
    # the time-frequency track of dominant mode to get
    # the corresponding `f_min` for `t_obs_start`.
    freq_array = np.logspace(-4, 0, num=10)
    t_array = np.zeros(len(freq_array))
    for i in range(len(freq_array)):
        t_array[i] = chirptime(m1=m1, m2=m2, f_lower=freq_array[i])
    tf_track = interp1d(t_array, freq_array)
    return tf_track

def bbhx_fd(ifos=None, run_phenomd=True,
            ref_frame='LISA', sample_points=None, **params):

    if ifos is None:
        raise Exception("Must define data streams to compute")

    from pycbc.types import FrequencySeries, Array
    from pycbc import pnutils

    m1 = params['mass1']
    m2 = params['mass2']
    a1 = params['spin1z']
    a2 = params['spin2z']
    dist = pnutils.megaparsecs_to_meters(params['distance'])
    phi_ref = params['coa_phase']
    f_ref = 0 # This is now NOT standard LAL convention!
    inc = params['inclination']
    lam = params['eclipticlongitude']
    beta = params['eclipticlatitude']
    psi = params['polarization']
    t_ref = params['tc']
    t_obs_start = params['t_obs_start'] # in seconds

    if ('f_lower' not in params) or (params['f_lower'] < 0):
        # the default value of 'f_lower' in PyCBC is -1.
        tf_track = interpolated_tf(m1, m2)
        t_max = chirptime(m1=m1, m2=m2, f_lower=1e-4)
        if t_obs_start > t_max:
            # Avoid "above the interpolation range" issue.
            f_min = 1e-4
        else:
            f_min = tf_track(t_obs_start) # in Hz
    else:
        f_min = params['f_lower'] # in Hz
        tf_track = interpolated_tf(m1, m2)
        t_max = chirptime(m1=m1, m2=m2, f_lower=1e-4)
        if t_obs_start > t_max:
            f_min_tobs = 1e-4
        else:
            f_min_tobs = tf_track(t_obs_start) # in Hz
        if f_min < f_min_tobs:
            err_msg = f"Input 'f_lower' is lower than the value calculated from 't_obs_start'."

    # We want to cache the waveform generator, but as it takes a mass dependent
    # start frequency as input this is hard.
    # To solve this we *round* the *logarithm* of this mass-dependent start
    # frequency. The factor of 25 ensures reasonable spacing while doing this.
    # So we round down to the nearest 1/25 of the logarithm of the frequency
    log_mf_min = int(math.log(f_min*MTSUN_SI*(m1+m2)) * 25)
    wave_gen = get_waveform_genner(log_mf_min, run_phenomd=run_phenomd)

    if ref_frame == 'LISA':
        # Transform to SSB frame
        t_ref, lam, beta, psi = LISA_to_SSB(
            t_ref,
            lam,
            beta,
            psi
        )
    elif ref_frame == 'SSB':
        # Don't need to update variable names
        pass
    else:
        err_msg = f"Don't recognise reference frame {ref_frame}. "
        err_msg = f"Known frames are 'LISA' and 'SSB'."

    if sample_points is None:
        # It's likely this will be called repeatedly with the same values
        # in many applications.
        if 'f_final' in params and params['f_final'] != 0:
            freqs = cached_arange(0, params['f_final'], 1/t_obs_start)
        else:
            raise Exception("Please set 'f_final' in **params.")
    else:
        freqs = sample_points

    modes = [(2,2)] # More modes if not phenomd
    direct = False # See the BBHX documentation
    fill = True # See the BBHX documentation
    squeeze = True # See the BBHX documentation
    length = 1024 # An internal generation parameter, not an output parameter
    shift_t_limits = False # Times are relative to merger
    t_obs_end = 0.0 # Generates ringdown as well!

    wave = wave_gen(m1, m2, a1, a2,
                    dist, phi_ref, f_ref, inc, lam,
                    beta, psi, t_ref, freqs=freqs,
                    modes=modes, direct=direct, fill=fill, squeeze=squeeze,
                    length=length, t_obs_start=t_obs_start/YRSID_SI,
                    t_obs_end=t_obs_end,
                    shift_t_limits=shift_t_limits)[0]

    wanted = {}

    if 'LISA_A' in ifos:
        wanted['LISA_A'] = 0
    if 'LISA_E' in ifos:
        wanted['LISA_E'] = 1
    if 'LISA_T' in ifos:
        wanted['LISA_T'] = 2

    output = {}
    # Convert outputs to PyCBC arrays
    if sample_points is None:
        length_of_wave = t_obs_start
        loc_of_signal_merger_within_wave = t_ref % length_of_wave
        df = 1 / t_obs_start

        for channel, tdi_num in wanted.items():
            output[channel] = FrequencySeries(
                wave[tdi_num],
                delta_f=df,
                epoch=params['tc'] - loc_of_signal_merger_within_wave,
                copy=False
            )
    else:
        for channel, tdi_num in wanted.items():
            output[channel] = Array(wave[tdi_num], copy=False)

    return output
