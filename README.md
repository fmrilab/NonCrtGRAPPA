# Non-Cartesian GRAPPA

Reference implementation of:\
[A GRAPPA algorithm for arbitrary 2D/3D non-Cartesian sampling trajectories with rapid calibration](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.27801)\
doi: [10.1002/mrm.27801](https://doi.org/10.1002/mrm.27801)

cite as:
```
@article{Luo2019GRAPPA,
author = {Luo, Tianrui and Noll, Douglas C. and Fessler, Jeffrey A. and Nielsen, Jon-Fredrik},
title = {A GRAPPA algorithm for arbitrary 2D/3D non-Cartesian sampling trajectories with rapid calibration},
journal = {Magnetic Resonance in Medicine},
doi = {10.1002/mrm.27801},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.27801},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.27801}
}
``` 

setup_NonCrtGRAPPA.m does configurations.

Navigate to './GRAPPA/test' for running test file grappaTestScript.m

The algorithm is implemented as stand-alone, no dependency is required.\
However, [MIRT](http://web.eecs.umich.edu/~fessler/irt/fessler.tgz) box is needed for running grappaTestScript.m
