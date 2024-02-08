import functools
import math
import numpy as np
from scipy.interpolate import interp1d
from bbhx.utils.constants import MTSUN_SI, YRSID_SI
from bbhx.waveformbuild import BBHWaveformFD
from pycbc.coordinates import TIME_OFFSET_20_DEGREES, lisa_to_ssb, ssb_to_lisa


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


# Value to be used if no f_lower information is provided, and a limit for
# things like interpolation
_MINIMUM_FLOWER = 1E-6
_LOG_MINIMUM_FLOWER = -6


def interpolated_tf(m1, m2):
    # Using findchirp_chirptime in PyCBC to calculate 
    # the time-frequency track of dominant mode to get
    # the corresponding `f_min` for `t_obs_start`.
    freq_array = np.logspace(_LOG_MINIMUM_FLOWER, 0, num=200)
    t_array = chirptime(m1=m1, m2=m2, f_lower=freq_array)
    tf_track = interp1d(t_array, freq_array)
    return tf_track


def bbhx_fd(ifos=None, run_phenomd=True,
            ref_frame='LISA', sample_points=None, **params):

    if ifos is None:
        raise Exception("Must define data streams to compute")

    from pycbc.types import FrequencySeries, Array
    from pycbc import pnutils

    m1 = np.float64(params['mass1'])
    m2 = np.float64(params['mass2'])
    a1 = np.float64(params['spin1z'])
    a2 = np.float64(params['spin2z'])
    inc = np.float64(params['inclination'])
    dist = np.float64(pnutils.megaparsecs_to_meters(params['distance']))
    f_ref = np.float64(params['f_ref'])
    phi_ref = np.float64(params['coa_phase']) # phase at f_ref
    if 't_offset' in params:
        if params['t_offset'] == 'TIME_OFFSET_20_DEGREES':
            t_offset = TIME_OFFSET_20_DEGREES
        else:
            t_offset = np.float64(params['t_offset']) # in seconds
    else:
        raise Exception("Must set `t_offset`, if you don't have a preferred value, \
                        please set it to be the default value %f, which will put LISA behind \
                        the Earth by ~20 degrees." % TIME_OFFSET_20_DEGREES)
    t_obs_start = np.float64(params['t_obs_start']) # in seconds

    if ref_frame == 'LISA':
        t_ref_lisa = np.float64(params['tc']) + t_offset
        lam = np.float64(params['eclipticlongitude'])
        beta = np.float64(params['eclipticlatitude'])
        psi = np.float64(params['polarization'])
        # Transform to SSB frame
        t_ref, lam, beta, psi = lisa_to_ssb(
            t_lisa=t_ref_lisa,
            longitude_lisa=lam,
            latitude_lisa=beta,
            polarization_lisa=psi,
            t0=0
        )
    elif ref_frame == 'SSB':
        t_ref = np.float64(params['tc']) + t_offset
        lam = np.float64(params['eclipticlongitude'])
        beta = np.float64(params['eclipticlatitude'])
        psi = np.float64(params['polarization'])
        # Don't need to update variable names,
        # because wave_gen receives parameters in SSB frame.
        t_ref_lisa, _, _, _ = ssb_to_lisa(
            t_ssb=t_ref,
            longitude_ssb=lam,
            latitude_ssb=beta,
            polarization_ssb=psi,
            t0=0
        )
    else:
        err_msg = f"Don't recognise reference frame {ref_frame}. "
        err_msg = f"Known frames are 'LISA' and 'SSB'."

    if ('f_lower' not in params) or (params['f_lower'] < 0):
        # the default value of 'f_lower' in PyCBC is -1.
        tf_track = interpolated_tf(m1, m2)
        t_max = chirptime(m1=m1, m2=m2, f_lower=_MINIMUM_FLOWER)
        if t_obs_start > t_max:
            # Avoid "above the interpolation range" issue.
            f_min = _MINIMUM_FLOWER
        else:
            f_min = tf_track(t_obs_start) # in Hz
    else:
        f_min = np.float64(params['f_lower']) # in Hz
        tf_track = interpolated_tf(m1, m2)
        t_max = chirptime(m1=m1, m2=m2, f_lower=_MINIMUM_FLOWER)
        if t_obs_start > t_max:
            f_min_tobs = _MINIMUM_FLOWER
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

    if sample_points is None:
        if 'delta_f' in params and params['delta_f'] > 0:
            # In PyCBC, default value for `delta_f` is -1.
            # Besides, soget_td_waveform_from_fd or _base_get_td_waveform_from_fd
            # will set nparams['delta_f'] = 1.0 / fudge_duration, and this is not
            # same as 1 / t_obs_start.
            df = np.float64(params['delta_f'])
        else:
            raise Exception("Please set 'delta_f' in **params.")
        # It's likely this will be called repeatedly with the same values
        # in many applications.
        if 'f_final' in params and params['f_final'] != 0:
            freqs = cached_arange(0, np.float64(params['f_final']), df)
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
        length_of_wave = 1 / df
        loc_of_signal_merger_within_wave = t_ref_lisa % length_of_wave

        for channel, tdi_num in wanted.items():
            output[channel] = FrequencySeries(
                wave[tdi_num],
                delta_f=df,
                epoch=t_ref_lisa - loc_of_signal_merger_within_wave,
                copy=False
            )
            # move the merge to the end of the vector
            output[channel] = output[channel].cyclic_time_shift(
                length_of_wave - loc_of_signal_merger_within_wave)
            output[channel].start_time -= t_offset
    else:
        for channel, tdi_num in wanted.items():
            output[channel] = Array(wave[tdi_num], copy=False)
            if t_offset != 0:
                # subtract t_offset from FD waveform
                output[channel] *= np.exp(2j*np.pi*sample_points*t_offset)

    return output
