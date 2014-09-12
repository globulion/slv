slv
===

The Solvshift (SLV) project is designed to
perform computations of frequency shift
of a solute immerged in a solvent in terms of discrete and
implicit solvatochromic models of prof. Minhaeng Cho [1-5] and its extended versions.
Available models are:
  * Solvatochromic Effective Fragment Potentials (SolEFP)[5]
  * Discrete electrostatic SolX models and its arbitrary contractions.
    The SolX models available are SolDMTP[4,6] and SolMMM[4,7]
  * Kirkwood-Onsager continuum solvatochromic model[4]


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

[5]: B. Błasiak and M. Cho, J. Chem. Phys. 140, 164107 (2014)

[6]: Distributed Multipole Moments (Mulliken, Chelpg, DMA, CAMM etc.)

[7]: Molecular Multipole Moments (centered at molecular origin, eg. center of mass)
