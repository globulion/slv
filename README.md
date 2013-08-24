slv
===

The Solvshift (SLV) project is designed to
perform computations of frequency shift
of a solute immerged in a solvent in terms of discrete and
implicit solvatochromic models of prof. Minhaeng Cho [1-4] and its extended versions.
Available models are:
  * Discrete electrostatic SolX models and its arbitrary contractions.
    The SolX models available are SolDMTP[5] and SolMMM[6]
  * Kirkwood-Onsager continuum model

To see the help box and the options type
```
slv --help
```
To install the package you have to type:
```
sudo python setup.py install
```

The [tutorial](https://github.com/globulion/slv/blob/master/USAGE.md "Title") is under preparation.

Good Luck!

References
----------

[1]: M. Cho, J. Chem. Phys. 130 (9), 094505 (2009)

[2]: M. Cho, J. Chem. Phys. 118 (8), 3480-3490 (2003)

[3]: H. Lee, J.-H. Choi and M. Cho, J. Chem. Phys. 137 (11), 114307 (2012)

[4]: B. Błasiak, H. Lee and M. Cho, J. Chem. Phys. 139, 044111 (2013)

[5]: Distributed Multipole Moments (Mulliken, Chelpg, DMA, CAMM etc.)

[6]: Molecular Multipole Moments (centered at molecular origin, eg. center of mass)
