import functools
import math
import numpy as np
from scipy.interpolate import interp1d
from bbhx.utils.constants import MTSUN_SI, YRSID_SI
from bbhx.waveformbuild import BBHWaveformFD
from bbhx.utils.transform import LISA_to_SSB, SSB_to_LISA


@functools.lru_cache(maxsize=128)
def get_waveform_genner(log_mf_min, run_phenomd=True, use_gpu=False):
    # See below where this function is called for description of how we handle
    # log_mf_min.
    mf_min = math.exp(log_mf_min/25.)
    wave_gen = BBHWaveformFD(amp_phase_kwargs=dict(run_phenomd=run_phenomd,
                                                   mf_min=mf_min),
                             use_gpu=use_gpu)
    return wave_gen

@functools.lru_cache(maxsize=10)
def cached_arange(start, stop, spacing):
    return np.arange(start, stop, spacing)

def chirptime(m1, m2, f_lower):
    from pycbc.waveform.spa_tmplt import findchirp_chirptime

    duration = findchirp_chirptime(m1=m1, m2=m2, fLower=f_lower, porder=7)
    return duration

def imr_duration(**params):
    # More accurate duration (within LISA frequency band) of the waveform,
    # including merge, ringdown, and aligned spin effects.
    # This is used in the time-domain signal injection in PyCBC.
    import warnings
    from pycbc.waveform.waveform import imrphenomd_length_in_time

    nparams = {'mass1':params['mass1'], 'mass2':params['mass2'],
               'spin1z':params['spin1z'], 'spin2z':params['spin2z'],
               'f_lower':params['f_lower']}
    time_length = np.float64(imrphenomd_length_in_time(**nparams))
    if time_length < 2678400:
        warnings.warn("Waveform duration is too short! Setting it to 1 month (2678400 s).")
        time_length = 2678400
    if time_length >= params['t_obs_start']:
        warnings.warn("Waveform duration is longer than data length! Setting it to `t_obs_start`.")
        time_length = params['t_obs_start']
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



# For now, always run phenom=False
def bbhx_fd(ifos=None, run_phenomd=False, use_gpu=False,
            ref_frame='LISA', sample_points=None, **params):

    if ifos is None:
        raise Exception("Must define data streams to compute")

    # If asking for anything but the (2,2) mode, run PhenomHM
    # if params['modes'] != [(2, 2)]:
    #     run_phenomd = False

    from pycbc.types import FrequencySeries, Array
    from pycbc import pnutils

    m1 = np.float64(params['mass1'])
    m2 = np.float64(params['mass2'])
    a1 = np.float64(params['spin1z'])
    a2 = np.float64(params['spin2z'])
    dist = np.float64(pnutils.megaparsecs_to_meters(params['distance']))
    phi_ref = np.float64(params['coa_phase'])
    f_ref = 0.0 # This is now NOT standard LAL convention!
    inc = np.float64(params['inclination'])
    lam = np.float64(params['eclipticlongitude'])
    beta = np.float64(params['eclipticlatitude'])
    psi = np.float64(params['polarization'])
    t_ref = np.float64(params['tc'])
    t_obs_start = np.float64(params['t_obs_start']) # in seconds

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
    wave_gen = get_waveform_genner(log_mf_min, run_phenomd=run_phenomd,
                                   use_gpu=use_gpu)

    if ref_frame == 'LISA':
        t_ref_lisa = t_ref
        # Transform to SSB frame
        t_ref, lam, beta, psi = LISA_to_SSB(
            t_ref_lisa,
            lam,
            beta,
            psi
        )
    elif ref_frame == 'SSB':
        # Don't need to update variable names
        t_ref_lisa, _, _, _ = SSB_to_LISA(
            t_ref,
            lam,
            beta,
            psi
        )
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

    # freqs = np.load('/home/connor/main/higher_modes_dev/likelihood_testing/freq.npy')

    if params['modes'] == 22.0:
        params['modes'] = [(2,2)]
    elif params['modes'] == 33.0:
        params['modes'] = [(3,3)]

    # If creating injection of many modes, or just single, compress = True
    # will do the same thing.
    compress = True #  If True, combine harmonics into single channel waveforms. (Default: True)
    # Need to give length if direct = False.
    direct = False # If True, directly compute the waveform without interpolation. (Default: False)
    fill = True # See the BBHX documentation
    squeeze = True # See the BBHX documentation
    length = 1024 # An internal generation parameter, not an output parameter
    shift_t_limits = False # Times are relative to merger
    t_obs_end = 0.0 # Generates ringdown as well!
    modes = params['modes'] # More modes if not phenomd


    # NOTE: This does not allow for the seperation of multiple modes into
    # their own streams.
    wave = wave_gen(m1, m2, a1, a2,
                    dist, phi_ref, f_ref, inc, lam,
                    beta, psi, t_ref, freqs=freqs,
                    modes=modes, direct=direct, fill=fill, squeeze=squeeze,
                    length=length, t_obs_start=t_obs_start/YRSID_SI,
                    t_obs_end=t_obs_end, compress=compress,
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
        loc_of_signal_merger_within_wave = t_ref_lisa % length_of_wave
        df = 1 / t_obs_start

        for channel, tdi_num in wanted.items():
            output[channel] = FrequencySeries(
                wave[tdi_num],
                delta_f=df,
                epoch=t_ref_lisa - loc_of_signal_merger_within_wave,
                copy=False
            )
    else:
        for channel, tdi_num in wanted.items():
            output[channel] = Array(wave[tdi_num], copy=False)

    return output
