!
!     EFP PARAMETER FILE FOR --- LITHIUM CATION (+1) --- MOLECULE
!
!     Contains:
!     ---------
!
!         structure    : YES
!         electrst     : YES
!         polarization : YES
!         wave-func    : YES
!           localized  : NO
!           canonical  : YES
!         vibration    : NO
!         NLO          : NO
!
!     Notes:
!     ------
!
!        * DMTP type: CHARGE=+1 AU
!        * DPOL type: TOTAL POLARIZABILITY
!        * LMO  type: No localization (1 orbital)
!
!                                       11 Feb 2015
!
 [ molecule ]                            
   name       = Lithium cation (+1)
   shortname  = Li+
   basis      = RHF/6-311++G**
   natoms     = 1
   nbasis     = 23
   nmos       = 1
   ncmos      = 1
   ndma       = 1 
   npol       = 1
 
 [ Atomic coordinates ]                            N= 3
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00

 [ Atomic numbers ]                                N= 1
                  3

 [ Atomic masses ]                                 N= 1
    12760.2223

 [ DMTP centers ]                                  N= 3
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00

 [ DMTP charges ]                                  N= 1
    1.0000000000E+00

 [ DMTP dipoles ]                                  N= 3
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00

 [ DMTP quadrupoles ]                              N= 6
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00

 [ DMTP octupoles ]                                N= 10
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00

 [ Polarizable centers ]                           N= 3
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00

 [ Distributed polarizabilities ]                  N= 9
    1.38192740E-01   0.00000000E+00   0.00000000E+00    0.00000000E+00   1.38192740E-01
    0.00000000E+00   0.00000000E+00   0.00000000E+00    1.38192740E-01

 [ LMO centroids ]                                 N= 3
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00

 [ Fock matrix ]                                   N= 1
   -2.79237289E+00

 [ AO->LMO matrix ]                                N= 23
    6.18403405E-01    4.28879122E-01    6.42974743E-16    1.93810973E-16   -6.13602002E-16     
   -1.51025922E-03   -7.03950011E-16   -6.42573299E-17    8.39440070E-16    1.17617104E-03
    3.45971392E-16   -9.64027343E-17   -6.37154075E-16   -3.40445853E-04   -8.99681447E-17
    5.15557834E-17    2.19034405E-16   -7.04759737E-04   -7.04759737E-04   -7.04759737E-04
   -4.12329435E-18   -9.17772494E-18   -9.81437011E-18                                   

