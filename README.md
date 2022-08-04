# gamma_calibration
Gamma calibration of TOFu system using Na-22 source.

Notes
-----
Data is taken using a Na-22 source which decays through either electron capture emitting a 1275 keV gamma or through positron emission resulting in the emission of 511 keV gammas when the positron annihilates. The Compton edges of the gamma spectrum can be used to calibrate the TOFOR detectors. The goal of these scripts is to find a offset (O) and a multiplier (M) to translate the acquired pulse waveform areas to the corresponding light yield produced in each TOFOR sub-detector. The equation is on the form
$$
E_{ly} = O + M\cdot A
$$
where $A$ is the area under the pulse waveform.