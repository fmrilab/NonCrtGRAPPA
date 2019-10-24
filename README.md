# Non-Cartesian GRAPPA

Reference implementation of:\
[A GRAPPA algorithm for arbitrary 2D/3D non-Cartesian sampling trajectories with rapid calibration](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.27801)\
doi: [10.1002/mrm.27801](https://doi.org/10.1002/mrm.27801)

Accecible [here](http://web.eecs.umich.edu/~fessler/papers/lists/files/jour/19/web/luo-19-aga.pdf).

cite as:

```bib
@article{Luo2019GRAPPA,
author = {Luo, Tianrui and Noll, Douglas C. and Fessler, Jeffrey A. and Nielsen, Jon-Fredrik},
title = {A GRAPPA algorithm for arbitrary 2D/3D non-Cartesian sampling trajectories with rapid calibration},
journal = {Magnetic Resonance in Medicine},
doi = {10.1002/mrm.27801},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.27801},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.27801}
}
```

## General comments

`setup_NonCrtGRAPPA.m` does configurations.\
Demos are provided in `./demo`.

``` md
# Parameter naming
P:    Image;
kP:   Fourier transform of P;
PS:   Coil images, i.e., P modulated with sensitivity maps;
kPS:  Fourier transform of PS, channel-wise;
ukPS: Under-sampled kPS, zero-filled;
fkPS: Full kPS, either fully sampled or reconstructed.
```

This repo has included binary test data files for basic accessibility in certain regions.\
Future binary data files will be available on <https://drive.google.com/drive/folders/1FxF5jcMhL8Z2IzB-i__mb4Q3wDOmi20R>.

## Dependencies

- MIRT: Github [link](https://github.com/JeffFessler/mirt), Package [link](http://web.eecs.umich.edu/~fessler/irt/fessler.tgz).

## Compatible MATLAB versions  
For MATLAB < R2018a, if LAPACK library is not correctly linked when building LS_fft_mex.c, one would need an additional flag 'lmwlapack' as: 

```matlab
>> mex LS_fft_mex.c -lmwlapack
```

