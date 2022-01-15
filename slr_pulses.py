#!/usr/bin/env python3
# Copyright: Rudraksha Dutta Majumdar, 2021
"""
Author: R.D Majumdar
Date: 2021Dec17
Brief: Functions to generate SLR pulses
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import sigpy
import sigpy.mri.rf as rf

# SLR pulse parameters
N = 256  # Number of time points
TBW = 9.62  # time-bandwidth product
PULSE_LEN = 1  #ms
PULSE_BW = TBW/PULSE_LEN #kHz
PBR = 0.01  # pass-band ripple
SBR = 0.01  # stop-band ripple
PULSE_TYPE = "ex" # inv, se, ex, sat, st (small-tip)
FILTER_TYPE = "min" # 'pm', 'ls', 'ms', 'min', 'max'

# Root-flipped pulse parameters
ROOT_FLIP = False
if PULSE_TYPE == "ex" or PULSE_TYPE == "sat":
    ROOT_FLIP_ANGLE = 90
elif PULSE_TYPE == "inv" or PULSE_TYPE == "se":
    ROOT_FLIP_ANGLE = 180
else:
    ROOT_FLIP_ANGLE = 180

# Multiband pulse parameters
MULTI_BAND = False
N_BANDS = 2
PHS_TYPE = 'quad_mod' # for n_bands >= 3 only: phs_mod, amp_mod, 
                        # for all n_bands: quad_mod, or 'None'
BAND_SEP = 6*TBW # separated by BAND_SEP slice widths

# Saving
SAVE_PULSE = False
BASE_PATH = "C:/Users/RudrakshaMajumdar/Documents/GitHub/rf-bloch-simulator/saved_rf_pulses"
NAME_PULSE = "SLR_test"

def slr_rootflip(flip=ROOT_FLIP_ANGLE):
    """
    flip: Target flip angle in degrees
    """
    flip_rad = np.deg2rad(flip)
    # code extracted from dzrf():
    [bsf, d1, d2] = rf.slr.calc_ripples(PULSE_TYPE, PBR, SBR)
    b = rf.slr.dzmp(N, TBW, d1, d2)
    b = b[::-1]
    b = bsf*b

    # root flipping the pulse, using the b polynomial
    [am_rootflip, bRootFlipped] = rf.slr.root_flip(b, d1, flip_rad, TBW)

    return am_rootflip, bRootFlipped

def slr_pulse(
    num=N, time_bw=TBW, 
    ptype=PULSE_TYPE, ftype=FILTER_TYPE, 
    d_1=PBR, d_2=SBR, 
    root_flip=ROOT_FLIP, 
    multi_band = MULTI_BAND,
    n_bands = N_BANDS,
    phs_type = PHS_TYPE,
    band_sep = BAND_SEP
    ):
    """Use Shinnar-Le Roux algorithm to generate pulse"""
    if root_flip is False:
        complex_pulse = rf.dzrf(n=num, tb=time_bw, ptype=ptype, ftype=ftype, d1=d_1, d2=d_2)
        amp_arr = complex_pulse
        
    else:
        amp_arr, b_rootflip = slr_rootflip(ROOT_FLIP_ANGLE)

    phs_arr = np.zeros(num)
    for idx in range(num):
        if amp_arr[idx] < 0:
            phs_arr[idx] = 180
        else:
            phs_arr[idx] = 0
    
    if multi_band is True:
        amp_arr = rf.multiband.mb_rf(amp_arr, n_bands, band_sep, phs_type)
    
    # prepare pulse for instrument, which takes absolute only
    # cast negative values to positive
    amp_arr_abs = np.abs(amp_arr)
    # shift amplitude such that the lowest value is 0
    amp_arr_abs = amp_arr_abs - amp_arr_abs.min()
    # fold back phase when it exceeds 360
    phs_arr = phs_arr % 360
    freq_arr = (np.diff(phs_arr)/num)/360

    return amp_arr, freq_arr, phs_arr, amp_arr_abs

def ck_to_mag(alpha, beta, se=False):
    """Convert Shinnar-Le Roux parameters to magnetization"""
    if se is False:
        m_z = 1-2*np.abs(beta)**2
        m_xy = 2 * np.conj(alpha) * beta
    else:
        m_z = 1-2*np.abs(beta)**2
        m_xy = 1j * (np.conj(alpha) **2 + beta**2)

    
    return m_z, m_xy

if __name__ == "__main__":
    am, fm, pm, am_abs = slr_pulse()
    # simulations
    if MULTI_BAND is True:
        space = np.linspace(-20*TBW, 20*TBW, 2000)
    else:
        space = np.linspace(-2*TBW, 2*TBW, 200)
    alpha, beta = sigpy.mri.rf.sim.abrm(am, space, balanced=False)
    
    if PULSE_TYPE == "se":
        m_z, m_xy = ck_to_mag(alpha, beta, se=True)
    else:
        m_z, m_xy = ck_to_mag(alpha, beta, se=False)


    fig = plt.figure()
    time_ax = np.linspace(0, PULSE_LEN, N)
    plt_am = fig.add_subplot(221)
    plt_pm = fig.add_subplot(222)
    plt_mz = fig.add_subplot(223)
    plt_mxy = fig.add_subplot(224)

    plt_am.plot(time_ax, am.real, 'b', label = "B1x", alpha = 0.3)
    plt_am.plot(time_ax, am.imag, 'r', label = "B1y", alpha = 0.3)
    plt_am.plot(time_ax, np.abs(am), 'k', label = "Abs.")
    plt_am.set_ylabel("Amplitude")
    plt_am.set_xlabel("time (ms)")
    plt_am.set_title("Time-bandwidth Product = " + str(TBW))
    plt_am.legend()
    plt_am.grid()

    plt_pm.plot(time_ax, pm, 'b')
    plt_pm.set_ylabel("Phase (deg.)")
    plt_pm.set_xlabel("time (ms)")
    plt_pm.grid()

    plt_mz.plot(space, m_z.real, 'b', label="Real", alpha = 0.3)
    plt_mz.plot(space, m_z.imag, 'r', label="Imag.", alpha = 0.3)
    plt_mz.plot(space, np.abs(m_z), 'k', label="Abs.")
    plt_mz.set_ylabel("Mz")
    plt_mz.set_xlabel("Frequency/Slice")
    plt_mz.set_title("Z-magnetization Profile")
    plt_mz.grid()
    plt_mz.legend()

    
    plt_mxy.plot(space, m_xy.real, 'b', label="Real", alpha = 0.3)
    plt_mxy.plot(space, m_xy.imag, 'r', label="Imag.", alpha = 0.3)
    plt_mxy.plot(space, np.abs(m_xy), 'k', label="Abs.")
    plt_mxy.set_ylabel("Mxy")
    plt_mxy.legend()
   
    plt_mxy.set_xlabel("Frequency/Slice")
    plt_mxy.set_title("XY-magnetization Profile")
    plt_mxy.grid()

    # plt.plot(fm)

    plt.tight_layout()
    plt.show()

    # Option to save the shape file
    if SAVE_PULSE is True:
        SAVE_PATH = os.path.join(BASE_PATH, NAME_PULSE)
        os.makedirs(SAVE_PATH, exist_ok=True)                

        # save the shape file
        file_name = os.path.join(SAVE_PATH, "rf_pulse_file")
        np.savez(file_name, am_abs, pm, fm)

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
