import functools
import math
import numpy as np
from scipy.interpolate import interp1d
from bbhx.utils.constants import MTSUN_SI, YRSID_SI, L_SI
from bbhx.waveformbuild import BBHWaveformFD
from pycbc.coordinates import TIME_OFFSET_20_DEGREES, lisa_to_ssb, ssb_to_lisa
from warnings import warn


def get_waveform_genner(log_mf_min=None, mf_min=None, run_phenomd=True):
    # See below where this function is called for description of how we handle
    # log_mf_min.
    if mf_min is None:
        mf_min = math.exp(log_mf_min/25.)
    wave_gen = BBHWaveformFD(
        amp_phase_kwargs=dict(
            run_phenomd=run_phenomd,
            mf_min=mf_min,
        ),
    )
    return wave_gen


@functools.lru_cache(maxsize=128)
def cached_get_waveform_genner(*args, **kwargs):
    """Cached version of get_waveform_genner"""
    return get_waveform_genner(*args, **kwargs)


@functools.lru_cache(maxsize=10)
def cached_arange(start, stop, spacing):
    return np.arange(start, stop, spacing)


@functools.lru_cache(maxsize=128)
def cached_freq_logspace(f_lower, num_interp):
    return np.logspace(math.log10(f_lower), 0, num=num_interp)


def chirptime(m1, m2, f_lower, m_mode=None):
    """Compute the chirptime.

    Defaults for the (2,2) mode.
    """
    from pycbc.waveform.spa_tmplt import findchirp_chirptime

    # Find the (2,2) mode duration
    duration = findchirp_chirptime(m1=m1, m2=m2, fLower=f_lower, porder=7)
    # If a mode is specified, convert to that mode
    if m_mode is not None:
        # See https://arxiv.org/abs/2005.08830, eqs. 2-3
        factor = (2 / (m_mode)) ** (-8 / 3)
        duration *= factor
    return duration


def validate_length_in_time(length_in_time, t_obs_start):
    """Validate the estimate length in time compared to t_obs_start.

    Parameters
    ----------
    length_in_time : float
        Estimated length in time.
    t_obs_start
        Observation start time.

    Returns
    -------
    float
        Valid length in time.
    """
    if length_in_time < 2678400:
        warn("Waveform duration is too short! Setting it to 1 month (2678400 s).")
        length_in_time = 2678400
    if length_in_time >= t_obs_start:
        warn("Waveform duration is longer than data length! Setting it to `t_obs_start`.")
        length_in_time = t_obs_start
    return length_in_time * 1.1


def bbhx_phenomd_length_in_time(**params):
    """Compute the length in for BBHx PhenomD"""
    from pycbc.waveform.waveform import get_imr_length

    time_length = np.float64(get_imr_length("IMRPhenomD", **params))
    return validate_length_in_time(time_length, params["t_obs_start"])


def bbhx_phenomhm_length_in_time(**params):
    """Compute the length in for BBHx PhenomHM"""
    from pycbc.waveform.waveform import get_hm_length_in_time

    # Max m mode is 4 for BBHx PhenomHM, use if the `mode_array` is not
    # specified.
    length_in_time = np.float64(
        get_hm_length_in_time("IMRPhenomD", 4, **params)
    )
    return validate_length_in_time(length_in_time, params["t_obs_start"])


def interpolated_tf(m1, m2, m_mode=None, num_interp=100, f_lower=1e-4):
    """Interpolate the time frequency-track.

    Defaults to the dominant (2,2) mode and uses :code:`chirptime` to compute
    the track.
    """
    # Using findchirp_chirptime in PyCBC to calculate
    # the time-frequency track of dominant mode to get
    # the corresponding `f_min` for `t_obs_start`.
    freq_array = cached_freq_logspace(f_lower, num_interp)
    t_array = chirptime(m1=m1, m2=m2, f_lower=freq_array, m_mode=m_mode)
    tf_track = interp1d(t_array, freq_array)
    return tf_track


def waveform_setup(**kwargs):
    if kwargs['approximant'] == "BBHX_PhenomD":
        if kwargs.get('mode_array') is not None  and len(kwargs['mode_array']) != 1:
            raise RuntimeError("BBHX_PhenomD only supports the (2,2) mode!")
        kwargs['mode_array'] = [(2, 2)]
        return _bbhx_fd(run_phenomd=True, **kwargs)
    elif kwargs['approximant'] == "BBHX_PhenomHM":
        if kwargs.get('mode_array') is None:
            kwargs['mode_array'] = [(2, 2), (2, 1), (3, 3), (3, 2), (4, 4), (4, 3)]
        return _bbhx_fd(run_phenomd=False, **kwargs)
    else:
        raise ValueError(f"Invalid approximant: {kwargs['approximant']}")


def _bbhx_fd(
    ifos=None,
    run_phenomd=True,
    ref_frame='LISA',
    tdi=None,
    sample_points=None,
    length=1024,
    direct=False,
    num_interp=100,
    interp_f_lower=1e-4,
    cache_generator=True,
    mf_min=None,
    **params
):
    
    """Function to generate frequency-domain waveforms using BBHx.
    
    Parameters
    ----------
    ifos : list
        List of interferometers
    run_phenomd : bool
        Flag passed to :code:`bbhx.waveformbuild.BBHWaveformFD` that determines
        if PhenomD or PhenomHM is used.
    tdi : {'1.5', '2.0}
        Version of TDI to use.
    ref_frame : {'LISA', 'SSB'}
        Reference frame.
    samples_points : numpy.ndarray, optional
        Array of frequencies for computing the waveform
    length : int
        Length parameter passed to BBHx. Must be specified if
        :code:`direct=False`. See BBHx documentation for more details.
    direct : bool
        See BBHx documentation.
    num_interp : int
        Number of interpolation points used for computing chirp time.
    interp_f_lower : float
        Lower frequency cutoff used for interpolation when computing the 
        chirp time.
    mf_min : float, optional
        Minimum frequency used by BBHx when performing interpolation. If
        not specified, the value will be set based on the total mass and
        minimum frequency.
    cache_generator : bool
        If true, the BBHx waveform generator is cached based on the computed
        ``mf_min``. Must be ``False`` if ``mf_min`` is specfied.
    
    Returns
    -------
    dict
        A dictionary containing the the waveforms for each interferometer.
    """

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
    mode_array = list(params["mode_array"])
    num_interp = int(num_interp)
    length = int(length) if length is not None else None
    interp_f_lower = float(interp_f_lower)
    

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

    # We follow the convention used in LAL and set the frequency based on the
    # highest m mode. This means that lower m modes will start at later times.
    max_m_mode = max([mode[1] for mode in mode_array])
    if ('f_lower' not in params) or (params['f_lower'] < 0):
        # the default value of 'f_lower' in PyCBC is -1.
        t_max = chirptime(
            m1=m1, m2=m2, f_lower=interp_f_lower, m_mode=max_m_mode
        )
        if t_obs_start > t_max:
            # Avoid "above the interpolation range" issue.
            f_min = interp_f_lower
        else:
            tf_track = interpolated_tf(
                m1,
                m2,
                m_mode=max_m_mode,
                num_interp=num_interp,
                f_lower=interp_f_lower,
            )
            f_min = tf_track(t_obs_start) # in Hz
    else:
        f_min = np.float64(params['f_lower']) # in Hz
        t_max = chirptime(
            m1=m1, m2=m2, f_lower=interp_f_lower, m_mode=max_m_mode
        )
        if t_obs_start > t_max:
            f_min_tobs = interp_f_lower
        else:
            tf_track = interpolated_tf(
                m1,
                m2,
                m_mode=max_m_mode,
                num_interp=num_interp,
                f_lower=interp_f_lower,
            )
            f_min_tobs = tf_track(t_obs_start) # in Hz
        if f_min < f_min_tobs:
            err_msg = f"Input 'f_lower' is lower than the value calculated from 't_obs_start'."
            warn(err_msg, RuntimeWarning)

    # We want to cache the waveform generator, but as it takes a mass dependent
    # start frequency as input this is hard.
    # To solve this we *round* the *logarithm* of this mass-dependent start
    # frequency. The factor of 25 ensures reasonable spacing while doing this.
    # So we round down to the nearest 1/25 of the logarithm of the frequency
    log_mf_min = math.log(f_min*MTSUN_SI*(m1+m2)) * 25
    if cache_generator:
        if mf_min is not None:
            raise RuntimeError(
                "Cannot use `cache_generator` when `mf_min` is specified"
            )
        # Use int to round down
        wave_gen = cached_get_waveform_genner(
            int(log_mf_min),
            mf_min=None,
            run_phenomd=run_phenomd,
        )
    else:
        wave_gen = get_waveform_genner(
            log_mf_min,
            mf_min=mf_min,
            run_phenomd=run_phenomd,
        )

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

    # If creating injection of many modes, or just single, compress = True
    # will do the same thing.
    compress = True #  If True, combine harmonics into single channel waveforms. (Default: True)
    fill = True # See the BBHX documentation
    squeeze = True # See the BBHX documentation
    shift_t_limits = False # Times are relative to merger
    t_obs_end = 0.0 # Generates ringdown as well!

    # NOTE: This does not allow for the separation of multiple modes into
    # their own streams. All modes requested are combined into one stream.
    wave = wave_gen(
        m1, m2, a1, a2,
        dist, phi_ref, f_ref, inc, lam,
        beta, psi, t_ref,
        freqs=freqs,
        modes=mode_array,
        direct=direct,
        fill=fill,
        squeeze=squeeze,
        t_obs_start=t_obs_start / YRSID_SI,
        t_obs_end=t_obs_end,
        compress=compress,
        length=length,
        shift_t_limits=shift_t_limits,
    )
    # For some reason, the shape is different depending on if direct is True
    # or False.
    if not direct:
        wave = wave[0]

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

    # convert TDI version
    if str(tdi) == '2.0':
        from pycbc.psd.analytical_space import omega_length
        if sample_points is None:
            # assume all channels share the same sample_frequencies
            omega_len = omega_length(f=output[channel].sample_frequencies, len_arm=L_SI)
        else:
            omega_len = omega_length(f=sample_points, len_arm=L_SI)
        rescale = 2j*np.sin(2*omega_len)*np.exp(-2j*omega_len)
        for key in output:
            output[key] *= rescale
    if str(tdi) not in ['1.5', '2.0']:
        raise ValueError("The TDI version only supports '1.5' and '2.0' for now.")

    return output
