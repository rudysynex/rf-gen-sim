#!/usr/bin/env python3
# Copyright: Rudraksha Dutta Majumdar, 2021
"""
Author: R.D Majumdar
Date: 2021Dec17
Brief: Functions to generate multiband pulses
        Requires a base shape generated by any of 
        the other pulse functions (except SLR)
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt

from AFP_all import sech, AFP_rf
from AHP_all import AHP_rf
from AM_pulses import AM_rf
from BIR4 import BIR4_rf
from CHIRP import chirp_rf

"""
Available shaped pulse functions:
Amplitude Modulated (AM_rf): 'sinc', 'gaussian', 'hermite'
Adiabatic Full/Half Passage (AFP_rf/AHP_rf): 'HSn', 'tanh/tan'
B1 Insensitive Rotation (BIR4): 'tanh/tan', 'sech/tanh'
CHIRP Linear frequency sweep: chirp_rf
"""
BASE_RF_FILE = "C:/Users/RudrakshaMajumdar/Documents/GitHub/rf-bloch-simulator/saved_rf_pulses/gaussian1_2ms/rf_pulse_file.npz"
BASE_RF_PULSE = np.load(BASE_RF_FILE)

PULSE_LENGTH = 2000 #us
NR_OF_POINTS = len(BASE_RF_PULSE['arr_0'])
TIME_REVERSAL = "on"
BLOCH_SIEGERT_COMP = "off"
FREQ_1 = 1000 #H
FREQ_2 = -1000
PHS_1 = 0
PHS_2 = 0
RF_AMPLITUDE = 1000 # Hz

SAVE_PULSE = True
BASE_PATH = "C:/Users/RudrakshaMajumdar/Documents/GitHub/rf-bloch-simulator/saved_rf_pulses"
NAME_PULSE = "multifreq_gauss_2ms"

def multi_freq(
    pulse_length = PULSE_LENGTH,
    shape_pts = NR_OF_POINTS,
    time_reversal = TIME_REVERSAL,
    bs_compensation = BLOCH_SIEGERT_COMP,
    freq_1 = FREQ_1,
    freq_2 = FREQ_2,
    phs_1 = PHS_1,
    phs_2 = PHS_2,
    rf_amplitude = RF_AMPLITUDE
):

    # Generate the base RF pulse shape to modulate
    # rf_amp, freq_mod, phs_mod, _ = AM_rf(func="gaussian",pulse_length=pulse_length,shape_pts=shape_pts,trunc_lvl=1)

    rf_amp = BASE_RF_PULSE['arr_0']
    phs_mod = BASE_RF_PULSE['arr_1']
    freq_mod = BASE_RF_PULSE['arr_2']

    rf_dwelltime = 1e-6 * (pulse_length/shape_pts)  # convert to seconds
    time = np.arange(0, pulse_length, pulse_length/shape_pts)

    # Convert RF phase from deg to rad
    phs_mod = phs_mod * (np.pi/180)

    # Time reversal of pulse 2 relative to pulse 1
    if time_reversal == "on":
        rf_amp_1 = rf_amp
        phs_mod_1 = phs_mod
        rf_amp_2 = np.flip(rf_amp)
        phs_mod_2 = np.flip(phs_mod)
    elif time_reversal == "off":
        rf_amp_1 = rf_amp
        phs_mod_1 = phs_mod
        rf_amp_2 = rf_amp
        phs_mod_2 = phs_mod
    else:
        rf_amp_1 = rf_amp
        phs_mod_1 = phs_mod
        rf_amp_2 = rf_amp
        phs_mod_2 = phs_mod
        print("Warning: No time reversal option specified. Continuing without time reversal")

    # Memory allocation
    Mx_final_1 = np.zeros(shape_pts)
    Mx_final_2 = np.zeros(shape_pts)
    My_final_1 = np.zeros(shape_pts)
    My_final_2 = np.zeros(shape_pts)

    M = np.zeros([3,1])
    # Initial Magnetization
    Mx_0 = 1.0
    My_0 = 0.0
    Mz_0 = 0.0

    if bs_compensation == "on":
        # Calculate the effect of pulse 1 with freq_1 at freq_2
        phs_1_freq_1 = phs_1 + 2 * np.pi * freq_1 * (time / 1e6)
        # Scale RF pulse over the entire RF range and convert from Hz to rad/s
        rf_amp = 2 * np.pi * rf_amplitude * rf_amp_1.reshape(-1,1)

        M[0,0] = Mx_0
        M[1,0] = My_0
        M[2,0] = Mz_0

        R = np.identity(3)

        # Convert frequency offset from hz to rad/s
        rf_offset = 2 * np.pi * freq_1

        for rf_pulse_counter in range(shape_pts):
            term_0 = rf_amp[rf_pulse_counter] ** 2
            term_1 = rf_offset ** 2

            #B_effective
            Be = np.sqrt(term_0 + term_1) * rf_dwelltime
            alpha = np.arctan2(rf_offset, rf_amp[rf_pulse_counter])

            cosBe = np.cos(Be)
            sinBe = np.sin(Be)

            cos1alpha = np.cos(alpha)
            cos2alpha = np.cos(alpha) * np.cos(alpha)
            sin1alpha = np.sin(alpha)
            sin2alpha = np.sin(alpha) * np.sin(alpha)

            cos1phi = np.cos(phs_1_freq_1[rf_pulse_counter])
            cos2phi = np.cos(phs_1_freq_1[rf_pulse_counter]) * np.cos(phs_1_freq_1[rf_pulse_counter])
            sin1phi = np.sin(phs_1_freq_1[rf_pulse_counter])
            sin2phi = np.sin(phs_1_freq_1[rf_pulse_counter]) * np.sin(phs_1_freq_1[rf_pulse_counter])

            # Construct the total rotation matrix
            R[0,0] = cos2phi*(cos2alpha + cosBe*sin2alpha) + cosBe*sin2phi
            R[1,0] = -sin1alpha*sinBe + sin1phi*cos1phi*cos2alpha*(1 - cosBe)
            R[2,0] = cos1alpha*(cos1phi*sin1alpha*(cosBe - 1) - sinBe*sin1phi)

            R[0,1] = sin1alpha*sinBe + sin1phi*cos1phi*cos2alpha*(1 - cosBe)
            R[1,1] = cosBe*cos2phi + (cos2alpha + cosBe*sin2alpha)*sin2phi
            R[2,1] = sin1alpha*cos1alpha*sin1phi*(cosBe - 1) + cos1alpha*cos1phi*sinBe

            R[0,2] = cos1alpha*(cos1phi*sin1alpha*(cosBe - 1) + sinBe*sin1phi)
            R[1,2] = sin1alpha*cos1alpha*sin1phi*(cosBe - 1) - cos1alpha*cos1phi*sinBe
            R[2,2] = cosBe*cos2alpha + sin2alpha

            M = R @ M

            Mx_final_1[rf_pulse_counter] = M[0,0]
            My_final_1[rf_pulse_counter] = M[1,0]
        
        # Calculate phase evolution during pulse 1
        phs_from_pulse1 = np.arctan2(My_final_1,Mx_final_1)
        phs_from_pulse1 = np.unwrap(phs_from_pulse1)
        phs_from_chem_shft1 = 2 * np.pi * freq_1 * (time / 1e6)
        phs_bs_1 = phs_from_pulse1 + phs_from_chem_shft1

        # Calculate the effect of pulse 2 with freq_2 at freq_1
        phs_2_freq_2 = phs_2 + 2 * np.pi * freq_2 * (time / 1e6)
        # Scale RF pulse over the entire RF range and convert from Hz to rad/s
        rf_amp = 2 * np.pi * rf_amplitude * rf_amp_2.reshape(-1,1)

        M[0,0] = Mx_0
        M[1,0] = My_0
        M[2,0] = Mz_0

        R = np.identity(3)

        # Convert frequency offset from hz to rad/s
        rf_offset = 2 * np.pi * freq_1

        
        for rf_pulse_counter in range(shape_pts):
            term_0 = rf_amp[rf_pulse_counter] ** 2
            term_1 = rf_offset ** 2

            #B_effective
            Be = np.sqrt(term_0 + term_1) * rf_dwelltime
            alpha = np.arctan2(rf_offset, rf_amp[rf_pulse_counter])

            cosBe = np.cos(Be)
            sinBe = np.sin(Be)

            cos1alpha = np.cos(alpha)
            cos2alpha = np.cos(alpha) * np.cos(alpha)
            sin1alpha = np.sin(alpha)
            sin2alpha = np.sin(alpha) * np.sin(alpha)

            cos1phi = np.cos(phs_2_freq_2[rf_pulse_counter])
            cos2phi = np.cos(phs_2_freq_2[rf_pulse_counter]) * np.cos(phs_2_freq_2[rf_pulse_counter])
            sin1phi = np.sin(phs_2_freq_2[rf_pulse_counter])
            sin2phi = np.sin(phs_2_freq_2[rf_pulse_counter]) * np.sin(phs_2_freq_2[rf_pulse_counter])

            # Construct the total rotation matrix
            R[0,0] = cos2phi*(cos2alpha + cosBe*sin2alpha) + cosBe*sin2phi
            R[1,0] = -sin1alpha*sinBe + sin1phi*cos1phi*cos2alpha*(1 - cosBe)
            R[2,0] = cos1alpha*(cos1phi*sin1alpha*(cosBe - 1) - sinBe*sin1phi)

            R[0,1] = sin1alpha*sinBe + sin1phi*cos1phi*cos2alpha*(1 - cosBe)
            R[1,1] = cosBe*cos2phi + (cos2alpha + cosBe*sin2alpha)*sin2phi
            R[2,1] = sin1alpha*cos1alpha*sin1phi*(cosBe - 1) + cos1alpha*cos1phi*sinBe

            R[0,2] = cos1alpha*(cos1phi*sin1alpha*(cosBe - 1) + sinBe*sin1phi)
            R[1,2] = sin1alpha*cos1alpha*sin1phi*(cosBe - 1) - cos1alpha*cos1phi*sinBe
            R[2,2] = cosBe*cos2alpha + sin2alpha

            M = R @ M

            Mx_final_2[rf_pulse_counter] = M[0,0]
            My_final_2[rf_pulse_counter] = M[1,0]
        
        # Calculate phase evolution during pulse 2
        phs_from_pulse2 = np.arctan2(My_final_2,Mx_final_2)
        phs_from_pulse2 = np.unwrap(phs_from_pulse2)
        phs_from_chem_shft2 = 2 * np.pi * freq_1 * (time / 1e6)
        phs_bs_2 = phs_from_pulse2 + phs_from_chem_shft2

    # Recalculate phase evolutions
    if bs_compensation == "on":
        phs_1_freq_1 = phs_1 + 2 * np.pi * freq_1 * (time / 1e6) - phs_bs_2
        phs_2_freq_2 = phs_2 + 2 * np.pi * freq_2 * (time / 1e6) - phs_bs_1
    elif bs_compensation == "off":
        phs_1_freq_1 = phs_1 + 2 * np.pi * freq_1 * (time / 1e6)
        phs_2_freq_2 = phs_2 + 2 * np.pi * freq_2 * (time / 1e6)
    else:
        phs_1_freq_1 = phs_1 + 2 * np.pi * freq_1 * (time / 1e6)
        phs_2_freq_2 = phs_2 + 2 * np.pi * freq_2 * (time / 1e6)


    # Add additional phase offset
    phs_1_freq_1 = phs_1_freq_1 + (np.pi/180) * phs_1
    phs_2_freq_2 = phs_2_freq_2 + (np.pi/180) * phs_2

    # Calculate real and imaginary RF amplitudes
    rf_real_1 = rf_amp_1 * np.cos(phs_1_freq_1)
    rf_img_1 = rf_amp_1 * np.sin(phs_1_freq_1)

    rf_real_2 = rf_amp_2 * np.cos(phs_2_freq_2)
    rf_img_2 = rf_amp_2 * np.sin(phs_2_freq_2)

    # Calculate final RF amplitude, phase and frequency
    rf = np.sqrt((rf_real_1 + rf_real_2)**2 + (rf_img_1 + rf_img_2)**2)
    rf = rf/np.amax(rf)

    phs_mod = np.arctan2((rf_real_1 + rf_real_2), (rf_img_1 + rf_img_2))
    phs_mod = (180/np.pi) * np.unwrap(phs_mod)
    phs_mod = phs_mod % 360

    freq_mod = np.diff(phs_mod, append=freq_mod[shape_pts-1])
    freq_mod[shape_pts-1] = 0

    return(rf, freq_mod, phs_mod, time)
if __name__ == "__main__":
    rf, freq_mod, phs_mod, time = multi_freq()
    # Plotting
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
