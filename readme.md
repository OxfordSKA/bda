# Baseline dependent averaging scripts

This is a set of python scripts for simulations of baseline dependent 
averaging using [CASA](https://casa.nrao.edu), 
[OSKAR](https://oerc.ox.ac.uk/~ska/oskar2/), 
and [WSClean](https://sourceforge.net/projects/wsclean/).

Test simulations:
* [Test simulation 1](https://github.com/OxfordSKA/bda/wiki/test_sim_001)

Some benchmark results (work in progress) can be found 
[here](https://github.com/OxfordSKA/bda/wiki/benchmarks).

-----------------------------

Timings [1Feb2016]

Total run time of 839.5 s
    - simulation (104s)
    - corrupt (239s)
    - bda (80.7s)
    - ~20s*6 calibrations = 120s
    - imaging > 250s (using lots of combinations **)
    
** imaging costs will drop considerably if we don't image everything!
    

- Simulation1 (model): 50 times (20 over-sample) => 1000 time oskar time steps.
    * 104.3 s
- Simulation2 (model ref):  50 times (no sky)
    * 5.7 s
- Corrupt MS (corrupt + applycal):
    * copy ms: 4.9 s
    * create MODEL_DATA and CORRECTED_DATA columns: 45.1 s 
    * copy DATA to MODEL_DATA: 2.4 s
    * Create calibration table for corruptions: 23.4 s
    * Fill calibration table: 14.9 s
    * Apply corruptions: 146.2 s
    * Update DATA column (copy CORRECTED_DATA to DATA?): 1.4 s
    * Total: 239.4 s
- Average (collapse over-sample)
    * ???
- Add noise
    * ???
- Average (BDA) (mstrasnform)
    * Model -> model_bda: 23.4s
    * Corrupted -> corrupted_bda: 28.6 s 
    * Corrupted_noisy -> corrupted_noisy_bda: 28.7 s
- Expand BDA prior to calibration
    * ???
- Calibrate (gaincal + applycal)
    * corrupted: 13.5 s + 7.3 s
    * corrupted_bda: 12.6 s + 4.2 s
    * corrupted_bda_expanded: 11.4 s + 7.2 s
    * corrupted noisy: 11.3 s + 7.3 s
    * corrupted_noisy_bda: 9.2 s + 4.2 s
    * corrupted_noisy_bda_expanded: 11.7 s + 7.2 s
- BDA (the expanded data sets)
    * calibrated_bda_expanded: 28.3 s 
    * calibrated_noisy_bda_expanded: 28.4 s
- Imaging (CASA imtool) (uniform, on target, off target, DATA, CORRECTED, MODEL):
    * calibrated: (3.7 + 7.3 + 3.5 + 7.2 + 3.6 + 7.3)
    * calibrated bda: (2.1 + 4.1 + 2.1 + 4.1 + 2.1 + 4.3)
    * calibrated bda expanded: (3.7 + 7.1 + 3.5 + 7.1 + 3.7 + 7.4)
    * calibrated bda expanded bda: (2.1 + 4.1 + 2.1 + 4.0 + 2.1 + 4.3)
    * calibrated noisy: 3.4 + 7.1 + 3.6 + 7.2 + 3.8 + 7.4
    * calibrated noisy bda: 2.1 + 4.1 + 2.0 + 4.0 + 2.1 + 4.2
    * calibrated noisy bda expanded: 3.6 + 7.3 + 3.7 + 7.2 + 3.7 + 7.4
    * calibrated noisy bda expanded bda: 2.0 + 4.1 + 2.0 + 4.1 + 2.1 + 4.3
    * corrupted: 3.6 + 7.2 (uniform DATA, on and off target)
    * corrupted_bda: 2.1 + 4.1 (uniform DATA, on and off target)
    * corrupted_bda_expanded: 3.6 + 7.2 (uniform DATA, on and off target)
    * corrupted_noisy: 3.6 + 7.2 (uniform DATA, on and off target)
    * corrupted_noisy_bda: 2.0 + 4.0 (uniform DATA, on and off target)
    * corrupted_noisy_bda_expanded: 3.6 + 7.3 (uniform DATA, on and off target)
    * model: 3.4 + 6.9 (uniform DATA, off + on target)
    * model_bda: 2.1 + 4.2 (uniform DATA, off + on target)
    * model_ref:  3.6 + 7.2 (uniform DATA, off + on target)
- Diffs
    * ??



