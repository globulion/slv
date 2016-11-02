Solvshift - Solvatochromic Shift Quantum Chemistry Program
==========================================================

Bartosz Błasiak, 2012-2016

Description
-----------

The Solvshift (SLV) project is designed to develop a tool that enables performing fast and accurate 
computations of the interaction-induced vibrational property fluctuations of a chosen solute's vibrational degree of freedom.  
At present, the available code allows for the vibrational frequency shift predictions
relative to the gas-phase state, in which IR active spectator is isolated from other substances.
Discrete solvatochromic models [1-5] and its extended versions [6-9] are currently developed.
In particular, Solvshift implements:
  * Solvatochromic Effective Fragment Potentials (SolEFP)[6-8]
  * Solvatochromic Shifts from Supermolecular Energy Decomposition Scheme (SolEDS) [6-9]
  * Discrete electrostatic, multipole-based solvatochromic models.
    Available are SolCAMM [5,8,9] models, their arbitrary contractions and SolMMM [5,8,9] models.
  * Kirkwood-Onsager continuum solvatochromic model [5]. This model is highly qualitative
    and is of predominantly didactic importance.

**Table 1.** Range of applications of various models implemented in Solvshift.

| Method      | Intermolecular Interaction | Accuracy              | Level of Theory  | Target Systems                | Purpose                              | 
| ----------- | -------------------------- | --------------------- | ---------------- | ----------------------------- | ------------------------------------ | 
| SolEDS      | Supermolecular Approach    | High (quantitative)   | HF and MP2       | Small clusters                | Validation of simplified models      |
| SolEFP      | Perturbation Theory        | Low (qualitative)     | HF               | Bulk solutions, proteins      | Simulations of vibrational spectra   |
| SolCAMM     | Multipole Expansion        | Low (qualitative)     | HF, MP2, CC, DFT | Model systems, bulk, proteins | Simulations of vibrational spectra   |
| SolMMM      | Multipole Expansion        | good only for molecular properties | HF, MP2, CC, DFT | Single molecule  | Electrostatic solvatochromic properties (e.g. Stark tuning rates) | 

The [tutorial](https://github.com/globulion/slv/blob/master/USAGE.md "Title") is under preparation.

Good Luck!

References
----------

[1]: [A. D. Buckingham, *Trans. Faraday Soc.* **1960**, 56, 753-760](http://pubs.rsc.org/en/content/articlepdf/1960/tf/tf9605600753 "Title")

[2]: M. Cho, *J. Chem. Phys.* **2009**, 130, 094505

[3]: M. Cho, *J. Chem. Phys.* **2003**, 118, 3480-3490

[4]: H. Lee, J.-H. Choi and M. Cho, *J. Chem. Phys.* **2012** 137, 114307

[5]: B. Błasiak, H. Lee and M. Cho, *J. Chem. Phys.* **2013** 139, 044111

[6]: B. Błasiak and M. Cho, *J. Chem. Phys.* **2014** 140, 164107

[7]: B. Błasiak and M. Cho, *J. Chem. Phys.* **2015**, 143, 164111

[8]: B. Błasiak, A. W. Ritchie, L. J. Webb and M. Cho, *Phys. Chem. Chem. Phys.* **2016** 18, 18094-18111

[9]: M. Maj, C. Ahn, B. Błasiak, K. Kwak, H. Han and M. Cho, *J. Phys. Chem. B* **2016**, 120, 10167-10180

Molecular Multipole Moments (centered at molecular origin, eg. center of mass)
