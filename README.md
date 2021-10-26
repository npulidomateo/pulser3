# Pulser 2

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
- [x] Calculate the state
- [x] `back_to_z()`  (engineered pulse)
- [ ] Fix `back_to_z()` not giving the same gates if run from `cycle_benchmark()` and externally (see `delme.py`)
- [ ] Debug the MS gate!!
- [ ] `compiler` to produce hfgui / artiq compatible pulses
- [ ] `decompiler` to load the pulses into gates
