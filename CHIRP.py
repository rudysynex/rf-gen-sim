#!/usr/bin/env python3
# Copyright: Rudraksha Dutta Majumdar, 2021
"""
Author: R.D Majumdar
Date: 2021Dec17
Brief: Functions to generate frequency swept CHIRP pulses
"""
import os
import json
import numpy as np
import matplotlib.pyplot as plt

PULSE_LENGTH = 5000 #us
NR_OF_POINTS = 256
SMOOTHING_PC = 10 # % smoothing of the edges
SWEEP_BW = 4000
SWEEP_DIRECTION = "descend" # descend, ascend
FREQ_OFFSET = 0
INIT_PHASE = 0

SAVE_PULSE = True
BASE_PATH = "C:/Users/RudrakshaMajumdar/Documents/GitHub/rf-bloch-simulator/saved_rf_pulses"
NAME_PULSE = "chirp_down_5ms_4kHz"


def chirp_rf(
        pulse_length = PULSE_LENGTH,
        shape_pts = NR_OF_POINTS,
        smoothing_pc = SMOOTHING_PC,
        sweep_bw=SWEEP_BW,
        sweep_dir = SWEEP_DIRECTION,
        freq_offset=FREQ_OFFSET,
        init_ph=INIT_PHASE
    ):
    """
    Function to generate a Chirp pulse with a linear frequency sweep,
    given the following parameters:
        pulse_length = duration of RF pulse in us
        shape_pts = No. of points to use for the shape
        smoothing_pc = Percent of total points to smooth
            at the start and the end
        sweep_bw = frequency sweep bandwidth in Hz
        sweep_dir = Sweep direction, "ascend"/"descend"
        freq_offset = Offset frequency of the pulse in Hz
        init_ph = Initial phase of the pulse
    
    """

    # Time resolution
    time_res = pulse_length/(shape_pts)
    # Actual time axis
    time = np.arange(0, pulse_length, time_res)

    # Define apparent time, from -1 to +1
    # t_vec = np.linspace(-1, 1, shape_pts)

    # Define the linear chirp 
    rf = np.linspace(1, 1,shape_pts)

    # Apply smoothing to the ends of the chirp
    pts_to_smooth = np.trunc(shape_pts * (smoothing_pc/100))
    
    for idx in range (shape_pts):
        
        if idx < pts_to_smooth:
            scale = np.sin((np.pi/2) * (idx/pts_to_smooth))
            rf[idx] = rf[idx] * scale
        
        if idx > shape_pts - pts_to_smooth:
            scale = np.sin((np.pi/2) * ((shape_pts - idx)/pts_to_smooth))
            rf[idx] = rf[idx] * scale

    # Frequency modulation
    if sweep_dir == "ascend":
        freq_mod = np.linspace(-sweep_bw/2, sweep_bw/2, shape_pts)
    elif sweep_dir == "descend":
        freq_mod = np.linspace(sweep_bw/2, -sweep_bw/2, shape_pts)

    # Phase modulation = time itegral of frequency
    phs_mod = abs(360 * np.cumsum(freq_mod * (time_res/1e6)))
    # Force the 1st pont in phs_mod to zero
    phs_mod = phs_mod - phs_mod[0]

    

    # Add frequency offset to freq_mod
    freq_mod = freq_mod + freq_offset

    # Add freq offset as a phase ramp to phase modulation
    phs_mod = phs_mod + 360 * freq_offset * (time/1e6)

    # Add initial phase offset
    phs_mod = phs_mod + init_ph

    return (rf, freq_mod, phs_mod, time)
if __name__ == "__main__":
    #Example plot
    rf, freq_mod, phs_mod, time = chirp_rf(pulse_length=5000, sweep_dir="descend")
    # Convert degrees to radians for complex plot
    phs_mod_rad = np.deg2rad(phs_mod)
    rf_cos_mod = rf * np.cos(phs_mod_rad)
    rf_sin_mod = rf * np.sin(phs_mod_rad)


    # Plot the shape
    fig = plt.figure(figsize=[8,8])
    plt_rfamp = fig.add_subplot(221)
    plt_freq_mod = fig.add_subplot(222)
    plt_phs_mod = fig.add_subplot(223)
    plt_complex = fig.add_subplot(224)

    plt_rfamp.plot(time, rf, 'b')
    plt_rfamp.set_ylabel("RF amplitude")
    plt_rfamp.set_xlabel("Time (us)")


    plt_freq_mod.plot(time, freq_mod, 'b')
    plt_freq_mod.set_ylabel("Frequency (Hz)")
    plt_freq_mod.set_xlabel("Time (us)")


    plt_phs_mod.plot(time, phs_mod, 'b')
    plt_phs_mod.set_ylabel("Phase (deg.)")
    plt_phs_mod.set_xlabel("Time (us)")

    plt_complex.plot(time, rf_cos_mod, 'b', label = "B1x")
    plt_complex.plot(time, rf_sin_mod, 'r', label = "B1y")
    plt_complex.set_ylabel("RF amplitude")
    plt_complex.set_xlabel("Time (us)")
    plt_complex.set_ylim(-1,1)
    plt_complex.legend()

    plt_rfamp.grid()
    plt_freq_mod.grid()
    plt_phs_mod.grid()
    plt_complex.grid()

    plt.tight_layout()
    plt.show()

    # Option to save the shape file
    if SAVE_PULSE is True:
        SAVE_PATH = os.path.join(BASE_PATH, NAME_PULSE)
        os.makedirs(SAVE_PATH, exist_ok=True)                

        # save the shape file
        file_name = os.path.join(SAVE_PATH, "rf_pulse_file")
        np.savez(file_name, rf, phs_mod, freq_mod, rf_cos_mod, rf_sin_mod, time)

        # save the pulse parameters
        _global_dict = globals()
        _glob_names = list(_global_dict.keys())
        global_vars = {}

        for name in _glob_names:
            if not name.startswith("_") and isinstance(
                _global_dict[name], (str, bool, float, int, list, dict)
            ):
                global_vars[name] = _global_dict[name]
        
        params_file_name = os.path.join(SAVE_PATH, "rf_pulse_parameters.json")
        with open(params_file_name, "w+") as param_file:
            json.dump(global_vars, param_file, indent=4)
