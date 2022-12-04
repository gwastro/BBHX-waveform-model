import functools
from bbhx.utils.constants import *


def get_waveform_genner(mf_min, run_phenomd=True):
    from bbhx.waveformbuild import BBHWaveformFD

    wave_gen = BBHWaveformFD(amp_phase_kwargs=dict(run_phenomd=run_phenomd, mf_min=mf_min))
    return wave_gen

def chirptime(m1, m2, f_lower):
    from pycbc.waveform.spa_tmplt import findchirp_chirptime

    duration = findchirp_chirptime(m1=m1, m2=m2, fLower=f_lower, porder=7)
    return duration

def interpolated_tf(m1, m2):
    import numpy as np
    from scipy.interpolate import interp1d

    # Using findchirp_chirptime in PyCBC to calculate 
    # the time-frequency track of dominant mode to get
    # the corresponding `f_min` for `t_obs_start`.
    freq_array = np.logspace(-4, 0, num=100)
    t_array = []
    for freq in freq_array:
        t_array.append(chirptime(m1=m1, m2=m2, f_lower=freq))
    tf_track = interp1d(t_array, freq_array)
    return tf_track

def bbhx_fd(ifos=None, run_phenomd=True, f_final=0.1,
            ref_frame='LISA', sample_points=None, **params):

    if ifos is None:
        raise Exception("Must define data streams to compute")

    import numpy as np
    from pycbc.types import FrequencySeries, Array
    from pycbc import pnutils
    from bbhx.utils.transform import LISA_to_SSB

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

    if 'f_lower' in params and params['f_lower'] < 0:
        # the default value of 'f_lower' in PyCBC is -1.
        params.pop('f_lower')
    else: pass

    if 'f_lower' in params and 't_obs_start' not in params:
        f_min = params['f_lower'] # in Hz
        params['t_obs_start'] = chirptime(m1=m1, m2=m2, f_lower=f_min) # in seconds
    elif 'f_lower' not in params and 't_obs_start' in params:
        tf_track = interpolated_tf(m1, m2)
        t_max = chirptime(m1=m1, m2=m2, f_lower=1e-4)
        if params['t_obs_start'] > t_max:
            # Avoid "above the interpolation range" issue.
            f_min = 1e-4
        else:
            f_min = tf_track(t_obs_start) # in Hz
    elif 'f_lower' in params and 't_obs_start' in params:
        f_min_input = params['f_lower'] # in Hz
        tf_track = interpolated_tf(m1, m2)
        t_max = chirptime(m1=m1, m2=m2, f_lower=1e-4)
        if params['t_obs_start'] > t_max:
            f_min_tobs = 1e-4
        else:
            f_min_tobs = tf_track(t_obs_start) # in Hz
        if f_min_input <= f_min_tobs:
            f_min = f_min_input
        else:
            f_min = f_min_tobs        
    else:
        err_msg = f"Must provide 'f_lower', 't_obs_start', or both."

    wave_gen = get_waveform_genner(mf_min=f_min*MTSUN_SI*(m1+m2), run_phenomd=run_phenomd)

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
        freqs = np.arange(f_min, f_final, 1/params['t_obs_start'])
        print("Calculating values at frequencies (%d points): %s" % (len(freqs),freqs))
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
        length_of_wave = params['t_obs_start']
        loc_of_signal_merger_within_wave = t_ref % length_of_wave
        df = 1 / params['t_obs_start']
        padding = int(f_min/df)

        for channel, tdi_num in wanted.items():
            wave_padded = np.concatenate((np.zeros(padding, dtype=complex),wave[tdi_num]),axis=0)
            output[channel] = FrequencySeries(wave_padded, delta_f=df,
                                  epoch=params['tc'] - loc_of_signal_merger_within_wave)
    else:
        for channel, tdi_num in wanted.items():
            output[channel] = Array(wave[tdi_num])
    print(channel)
    print(f_min, f_final, freqs)
    print(output[channel].data[:10])
    return output
