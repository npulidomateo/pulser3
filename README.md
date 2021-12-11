# Pulser 3

Create 2-qubit pulse sequences for *cycle benchmarking* and read sequences from the generated `.dc` files or from `openQASM 2.0` files.

## Prerequisites
* numpy
* matplotlib
* qiskit

## Installation

### Method 1
Copy the file `pulser3.py` in the same folder as your script.

### Method 2
Add the following lines to the beginning of your script:
``` python
import sys
import pathlib
installation_path = pathlib.Path('/absolute/path/to/installation/folder')
sys.path.append(installation_path)
```
### Method 3 (Linux only)
Add the absolute path of the installation folder to your `$PYTHONPATH` enviroment variable.

For example add to your `.bashrc`

``` bash
export PYTHONPATH="${PYTHONPATH}:/absolute/path/to/installation/folder"
```

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

- [ ] Function to transform `U3Gate`s in minimum set of`RGate`s
- [ ] Read circuits from files
- [ ] Write jupyter tutorial 
  - [ ] Include conda environment file


## Qiskit version

| module                    | version   |
|---------------------------|-----------|
| qiskit-aer                | 0.9.1     |
| qiskit-ignis              | 0.6.0     |
| qiskit-ibmq-provider      | 0.17.0    |
| qiskit-aqua               | 0.9.5     |
| qiskit                    | 0.31.0    |
| qiskit-nature             | None      |
| qiskit-finance            | None      |
| qiskit-optimization       | None      |
| qiskit-machine-learning   | None      |
