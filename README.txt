Structure and chemical composition data in models 22L, 27LA, and 25M.

Data for "Three-dimensional Hydrodynamics Simulations of Pre-Collapse Shell
Burning in the Si and O-rich Layers",
by Yoshida, T., Takiwaki, T., Kotake, K., Takahashi K., Nakamura, K, Umeda, H.

Ver 2.0: 2020,12,17.
Ver 1.0: 2020, 4, 8.

- README.txt
  This file.


- sshctomach.f
  A program to make turbulent Mach number distribution from the SSH
  decomposition compnents a_nlm.
  Requested input data:
    - datafile2.txt or datafile3.txt or datafile4.txt in Yoshida et al. (2020)
    - The coordinate data of the convective region.
      The requested format is written in this program.


- mass300J1.txt
  300 nuclear species data.
  line 1: total species number = 300
------------------------------------------------------------------------
Byte-by-byte Description from line 2 of file: mass300J1.txt
------------------------------------------------------------------------
  Bytes Format Units   Explanations
------------------------------------------------------------------------
  3-  5 I4     ---     nuclear species number
  7- 11 A5     ---     nuclear species
 17- 23 F7.3   ---     mass number
 30- 31 I2     ---     proton number
 34- 35 I2     ---     neutron number
 39- 41 F3.1   ---     spin
 44- 51 F8.3   MeV     mass excess
 54- 57 A4     ---     Reference
------------------------------------------------------------------------


- strc22L directory contains
    strc_22L.txt, comp1_22L.txt, comp2_22L.txt, comp3_22L.txt.
- strc27LA directory contains
    strc_27LA.txt, comp1_27LA.txt, comp2_27LA.txt, comp3_27LA.txt.
- strc25M directory contains
    strc_25M.txt, comp1_25M.txt, comp2_25M.txt, comp3_25M.txt.


- strc(22L,27LA,25M).txt
  Progenitor structure data.
------------------------------------------------------------------------
Byte-by-byte Description from line 2 of file: strc(22L,27LA,25M).txt
------------------------------------------------------------------------
  Bytes Format Units            Explanations
------------------------------------------------------------------------
  3 - 6 I4     ---              mesh number
  7- 28 E22.14 Msun             Mass coordinate
 29- 50 E22.14 Msun             dMr
 51- 64 E14.6  ---              log(Radius/(Solar radius)): log(r/Rsun)
 65- 78 E14.6  log[g cm^{-3}]   log(Density)
 79- 92 E14.6  log[K]           log(Temperature)
 93-106 E14.6  log[dyn cm^{-2}] log(Pressure)
107-120 E14.6  cm s^{-1}        Radial velocity vr
121-134 E14.6  ---              Electron fraction (Ye)
------------------------------------------------------------------------
  

- comp1_(22L,27LA,25M).txt
  Progenitor composition data 1.
  There are 100 colums in each line.
  Line j, column k:
  Mass fraction of nuclei with species number k at mesh j.

- comp2_(22L,27LA,25M).txt
  Progenitor composition data 2.
  There are 100 colums in each line.
  Line j, column k:
  Mass fraction of nuclei with species number 100+k at mesh j.

- comp3_(22L,27LA,25M).txt
  Progenitor composition data 3.
  There are 100 colums in each line.
  Line j, column k:
  Mass fraction of nuclei with species number 200+k at mesh j.
