#!/usr/bin/env python3
# Copyright: Rudraksha Dutta Majumdar, 2021
"""
Author: R.D Majumdar
Date: 2021Dec17
Brief: Functions to generate adiabatic half-passage pulses
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt

# Pulse Parameters
PULSE_LENGTH = 5000 #us
NR_OF_POINTS = 256
PULSE_FUNCTION = "HSn"  # "HSn" OR "tanh/tan"
HYPSEC_ORDER = 4
SWEEP_BW = 4000
FREQ_OFFSET = 0
INIT_PHASE = 0

SAVE_PULSE = False
BASE_PATH = "C:/Users/RudrakshaMajumdar/Documents/GitHub/rf-bloch-simulator/saved_rf_pulses"
NAME_PULSE = "AHP_HS4_5ms_4kHz"

def arcsech(val):
    """
    Inverse hyperbolic secant
    """
    return np.arccosh(1. / val)

def sech(arg):
    """
    Hyperbolic secant
    """
    return 1 / np.cosh(arg)

def AHP_rf(
    pulse_length = PULSE_LENGTH,
    shape_pts = NR_OF_POINTS,
    func = PULSE_FUNCTION,
    beta = 5.3,
    n=HYPSEC_ORDER,
    k=1.52,
    sweep_bw=SWEEP_BW,
    freq_offset=FREQ_OFFSET,
    init_ph=INIT_PHASE
):
    """
    A function to generate an Adiabatic Half Passage pulses
    given the following parameters:
        func =  Modulation function, "HSn" or "tanh/tan"
                Defaults to an HS1 pulse
        pulse_length = duration of RF pulse in us
        shape_pts = No. of points to use for the shape
        beta = Rf modulation truncation factor, typically 5.3
        n = Power 'n' in Rf modulation. n=1 for HS1, n=8 for HS8 etc
        k = Freq modulation factor for tanh/tan, typically 1.52
        sweep_bw = frequency sweep bandwidth in Hz
        freq_offset = Offset frequency of the pulse in Hz
        init_ph = Initial phase of the pulse

    beta = Rf modulation truncation factor,
           typically 5.3 for HSn
                     10 for tanh/tan
    """

    # Time resolution
    time_res = pulse_length/(shape_pts)
    # Actual time axis
    time = np.arange(0, pulse_length, time_res)

    if func == "HSn":
        """
        Offset-independent adiabaticity (OIA) pulses based on
        hyperbolic secant functions raised to a power n.
        n = 1 leads to the original hyperbolic secant RF pulses.

        Original/relevant publications:
        1. M. S. Silver, R. I. Joseph, D. I. Hoult,
            J. Magn. Reson. 59, 347-351 (1984)
        2. A. Tannus, M. Garwood,
            J. Magn. Reson. 120A, 133-137 (1996)
        See also Chapter 5 in 'In Vivo NMR spectroscopy' book
        by Robin de Graaf
        """

        beta = 5.3

        # Define apparent time, from -1 to +1
        t_vec = np.linspace(0, 1, shape_pts)

        # Define the shape (first segment)
        rf = sech(beta*(1.0-t_vec)**n)

        # Normalize so that max val. is 1
        rf = rf/np.amax(rf)

        # Frequency modulation
        freq_mod = np.cumsum(rf * rf)
        # Normalize freq_mod to 0 to +1 range
        freq_mod = 1.0 - (freq_mod/np.amax(freq_mod))
        # Scale the freq_mod to full bandwidth
        freq_mod = (sweep_bw/2) * freq_mod

        # Phase modulation = time itegral of frequency
        phs_mod = 360 * np.cumsum(freq_mod * (time_res/1e6/2))
        # Force the 1st pont in phs_mod to zero
        phs_mod = phs_mod - phs_mod[0]


        # Add frequency offset to freq_mod
        freq_mod = freq_mod + freq_offset

        # Add freq offset as a phase ramp to phase modulation
        phs_mod = phs_mod + 360 * freq_offset * (time/1e6)

        # Add initial phase offset
        phs_mod = phs_mod + init_ph

    elif func == "tanh/tan":
        """
        Adiabatic half passage pulses based on tanh/tan modulation
        """

        beta = 10

        # Define apparent time, from 0 to +1
        t_vec = np.linspace(0, 1, shape_pts)

        # Define the shape (first segment)
        rf = np.tanh(beta * t_vec)

        # Normalize so that max val. is 1
        rf = rf/np.amax(rf)

        # Frequency modulation
        freq_mod = np.tan(k * (1.0 - t_vec))
        # Normalize freq_mod to 0 to +1 range
        freq_mod = freq_mod/np.amax(freq_mod)
        # Scale the freq_mod to full bandwidth
        freq_mod = (sweep_bw/2.0) * freq_mod

        # Phase modulation = time itegral of frequency
        phs_mod = 360 * np.cumsum(freq_mod * (time_res/1e6/2))
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
    rf, freq_mod, phs_mod, time = AHP_rf()
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
