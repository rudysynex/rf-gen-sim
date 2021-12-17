# rf-gen-sim
A set of functions to generate various shaped RF pulses, including SLR optimized pulses, and simulate their frequency response

Bloch simulation programs for RF pulse analysis, including programs
for generating shaped RF pulses
Author: Rudraksha Dutta Majumdar, 2021

Many of these functions are Python versions of Robin de Graaf's Matlab tool 'PulseWizard' (except the SLR pulses).
The SLR tool uses parts of the '[sigpy.mri.rf](https://sigpy.readthedocs.io/en/latest/mri_rf.html#module-sigpy.mri.rf)' modules.

The 'Bloch_simulator.ipynb' contains the main
simulation program to simulate the following:
1. Magnetization VS RF amplitude
2. Magnetization VS RF Offset Frequency
3. Magnetization VS Time, including 3D trajectory on a Bloch sphere

**Do not use this tool for SLR pulses!**

'Bloch_simulator.ipynb' requires the 'qutip' module to be installed
(http://qutip.org) in addition to Numpy and Matplotlib.

'AM_pulses.py' contains functions to generate Amplitude Modulated pulses (AM).

'AFP_all.py' contains functions to generate Adiabatic Full Passage pulses,
which are Amplitude and Frequency Modulated (AM/FM).

'AHP_all.py' contains functions to generate Adiabatic Half Passage pulses,
which are Amplitude and Frequency Modulated (AM/FM).

'BIR4.py' contains functions to generate B1-Insenstive Rotation pulses,
which are Amplitude and Frequency Modulated (AM/FM).

'CHIRP.py' contains functions to generate a chirp pulse with a
linear frequency sweep (AM/FM)

'multi_freq.py' contains functions to generate a dual band pulse

These pulse functions are designed to generate amplitude, frequency,
and phase tables as Numpy arrays.

**SLR PULSES**

The 'slr_pulses.py' program is a standalone script to both generate and simulate 
Shinnar - Le Roux optimized pulses, including options for root-flipping and multibanding.
Do not use the SLR pulses with the 'Bloch_simulator.ipynb' notebook
