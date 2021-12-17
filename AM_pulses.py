#!/usr/bin/env python3
# Copyright: Rudraksha Dutta Majumdar, 2021
"""
Author: R.D Majumdar
Date: 2021Dec17
Brief: Functions to generate amplitude modulated pulses
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt

# Pulse Parameters
PULSE_LENGTH = 132 #us
NR_OF_POINTS = 256
PULSE_FUNCTION = "square"  # sinc, square, gaussian, hermite
SINC_LOBES = 3
GAUSSIAN_TRUNCATION = 1 # in %
HERMITE_NUT_ANGLE = 90
FREQ_OFFSET = 0
INIT_PHASE = 0

SAVE_PULSE = True
BASE_PATH = "C:/Users/RudrakshaMajumdar/Documents/GitHub/rf-bloch-simulator/saved_rf_pulses"
NAME_PULSE = "square_132us"

def AM_rf(
    func = PULSE_FUNCTION,
    pulse_length = PULSE_LENGTH,
    shape_pts = NR_OF_POINTS,
    lobes = SINC_LOBES,
    trunc_lvl = GAUSSIAN_TRUNCATION,
    nut_angle = HERMITE_NUT_ANGLE,
    freq_offset = FREQ_OFFSET,
    init_ph = INIT_PHASE
):
    """
    A function to generate an Amplitude Modulated (AM) 
    pulse given the following parameters:
        func = Modulation function
        pulse_length = duration of RF pulse in us
        shape_pts = No. of points to use for the shape
        lobes = Number of lobes in the sinc shape
        trunc_lvl = Truncation level % for Gaussian pulse. 
                    1 or 10 are standard
        nut_angle = 90 or 180. Nutation angle for Hermite pulse
        freq_offset = Offset frequency of the pulse in Hz
        init_ph = Initial phase of the pulse
    """
    
    # Time resolution
    time_res = pulse_length/(shape_pts)
    # Actual time axis
    time = np.arange(0, pulse_length, time_res)

    
    if func == "sinc":
    
        # Define apparent time, from -1 to +1
        t_vec = np.linspace(-1, 1, shape_pts)
        # Define the shape
        rf = np.sin(np.pi*lobes*t_vec)/t_vec

        # Normalize so that max val. is 1
        rf = rf/np.amax(rf)

        # Frequency modulation
        freq_mod = np.zeros(shape_pts)
        freq_mod[0:shape_pts] = freq_offset


        # Calculate phase mod from freq mod
        phs_mod = 360 * freq_offset * (time/1e6)

        # Add initial phase offset
        phs_mod = phs_mod + init_ph
    
    elif func == "gaussian":
        
        # Define apparent time, from -1 to +1
        t_vec = np.linspace(-1, 1, shape_pts)

        # Convert truncation level (%) to exponential scaling
        gaussian_scale = 1.0 * np.log(trunc_lvl/100.0)

        # Define the shape
        rf = np.exp(gaussian_scale * t_vec * t_vec)

        # Normalize so that max val. is 1
        rf = rf/np.amax(rf)

        # Frequency modulation
        freq_mod = np.zeros(shape_pts)
        freq_mod[0:shape_pts] = freq_offset


        # Calculate phase mod from freq mod
        phs_mod = 360 * freq_offset * (time/1e6)

        # Add initial phase offset
        phs_mod = phs_mod + init_ph
        
    elif func == "hermite":
        
        # Define apparent time, from -1 to +1
        t_vec = np.linspace(-1, 1, shape_pts)



        # Define the shape
        if nut_angle == 90:
            rf = (1.0 - 6.003 * t_vec * t_vec) * np.exp(-9.0 * t_vec * t_vec)
        elif nut_angle == 180:
            rf = (1.0 - 8.604 * t_vec * t_vec) * np.exp(-9.0 * t_vec * t_vec)
        else:
            print("nut_angle can only be 90 or 180")

        # Normalize so that max val. is 1
        rf = rf/np.amax(rf)

        # Frequency modulation
        freq_mod = np.zeros(shape_pts)
        freq_mod[0:shape_pts] = freq_offset


        # Calculate phase mod from freq mod
        phs_mod = 360 * freq_offset * (time/1e6)

        # Add initial phase offset
        phs_mod = phs_mod + init_ph
    
    elif func == "square":

        # Define apparent time, from -1 to +1
        t_vec = np.linspace(-1, 1, shape_pts)
      
        # Define the shape
        rf = np.linspace(1, 1, shape_pts)

        # Frequency modulation
        freq_mod = np.zeros(shape_pts)
        freq_mod[0:shape_pts] = freq_offset


        # Calculate phase mod from freq mod
        phs_mod = 360 * freq_offset * (time/1e6)

        # Add initial phase offset
        phs_mod = phs_mod + init_ph
        
    else:
        print("ERROR: Please check 'func' argument")
        
        
    return (rf, freq_mod, phs_mod, time)

if __name__ == "__main__":
    #Example plot
    rf, freq_mod, phs_mod, time = AM_rf()
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
