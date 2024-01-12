## Cantera Adaptive Testing Documentation
This is the documentation used in generation of the results in multiple papers on the
generalized preconditioning approach for both gas phase and surface phase chemistry.
All study data and scripts can be found in the `study` directory.

### Getting Started
The easiest way to repeat or use my analyses is to install this package with

```sh
pip install -e .
```
from the github repository.

### Commands
A series of commands was developed and used for these studies, see `-h` for configuration options for each command.
- adaptive-testing: The primary interface for running problems.
- adaptive-utilities: A series of utilities for transforming and manipulating data.
- adaptive-yp: A plotting utilities to plot data in useful forms.
