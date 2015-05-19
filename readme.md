# BDA scripts
- 01_sim.py: Simulation script which generates corrupted data.
  - OSKAR sim
    - Parameterise by ...

- 02_cor.py
  - Corrupt data by generating calibration table and applying in CASA?
  - Data products:
    - Corrupted and uncorrupted MS.

- 03_avg.py
  - Average MS using mstransform or other?

- 04_cal.py
  - Calibrate and apply calibrations
  - Data products:
    - Calibrated MS

- 05_img.py
  - Image measurement sets (Corrupted, uncorrupted and calibrated)
    - wsclean or CASA clean

  - Generate any metrics
