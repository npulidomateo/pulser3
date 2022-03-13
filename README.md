# Pulser 3

Create 2-qubit pulse sequences for *cycle benchmarking* and read sequences from the generated `.dc` files or from `openQASM 2.0` files.

## Prerequisites
* numpy
* matplotlib
* qiskit

## Installation

1. Create a virtual environment e.g.
  
`conda create -n pulser python=3`

2. Install with pip

`pip install -e .`

## TODO

- [x] visualize sequences of pulses ~~using `matplotlib`~~
- [x] Introduce MS gate
- [x] Allow swith between MS gate and Identity (for debugging)
- [x] Calculate the state
- [x] `back_to_z()`  (engineered pulse)
  - [x] Fix `back_to_z()` not giving the same gates if run from `cycle_benchmark()` and externally (see `delme.py`)
- [x] Debug the MS gate!!
- [x] Make engineered pulses use multiples of $\pi$ instead of long decimal numbers
- [x] `compiler` to produce hfgui / ~~artiq~~ compatible pulses
  - [x] write `theta` as a function of `scan_pi`
    - [x] hfgui wants "." after floating point numbers --> `put_the_dot()`
    - [x] different ions have different pitimes `scan_pi_0` / `scan_pi_1`
  - [x] add the ACZS correction to SIA pulses
  - [x] write pulses into files (use `pathlib`)
  - [ ] ~~compile single-ion circuits (global or MMSB pulses?)~~
- [x] Use multiple cores
- [x] `decompiler` to load the pulses into gates
- [ ] `decompiler` fix ion_0 ion_1 convention!

### Wishlist

- [ ] Arbitrary, user specified MS operator
- [ ] Function to transform `U3Gate`s in minimum set of`RGate`s
- [ ] Read circuits from files
- [ ] Write jupyter tutorial 
  - [ ] Include conda environment file


## Qiskit version

| module                    | version   |
|---------------------------|-----------|
| qiskit-terra | 0.19.2 |
| qiskit-aer | 0.10.3 |
| qiskit-ignis | 0.7.0 |
| qiskit-ibmq-provider | 0.18.3 |
| qiskit-aqua | 0.9.5 |
| qiskit | 0.34.2 |
| qiskit-nature | None |
| qiskit-finance | None |
| qiskit-optimization | None |
| qiskit-machine-learning | None |
