#!/bin/bash

echo "Running the ldos calculation for MgB2 ..."

# clean up, but do some checks that we are in the right directory.
if [[ ! -f ../../build_potcar.sh ]]; then
	echo "Are we in the right directory? Something is wrong ... exiting."
	exit 1
fi
rm -f CHG CHGCAR CONTCAR DOSCAR EIGENVAL IBZKPT INCAR infile KPOINTS ldos.dat OSZICAR OUTCAR PCDAT POSCAR POTCAR POTCAR.spec REPORT vasprun.xml WAVECAR XDATCAR vasp_run.log

echo Building VASP input ...""
# Build the POTCAR file
(
cat << POTCAR
B
Mg_pv
POTCAR
) > ./POTCAR.spec

echo "Building POTCAR according to POTCAR.spec ..."
../../build_potcar.sh POTCAR.spec
if [ $? -ne 0 ]; then
	echo "Problem building POTCAR. Exiting ..."
	exit 1;
fi

# Build the INCAR file ...
( 
cat << INCAR
ALGO = Fast
ECUT = 510
EDIFF = 1e-07
IBRION = -1
ICHARG = 1
ISIF = 3
ISMEAR = 1
ISPIN = 1
LREAL = False
LWAVE = True
PREC = Accurate
SIGMA = 0.05
NWRITE = 3
INCAR
) > ./INCAR

# ... the POSCAR file ...
(
cat << POSCAR
Mg1 B2
   1.0000000000000000
     3.0519418355658581    0.0000000000000000   -0.0000000000000000
    -1.5259709177829290    2.6430591965863179    0.0000000000000000
    -0.0000000000000000    0.0000000000000000    3.5056459752153359
   B    Mg
     2     1
Direct
  0.6666669999999968  0.3333330000000032  0.5000000000000000
  0.3333330000000032  0.6666669999999968  0.5000000000000000
 -0.0000000000000000 -0.0000000000000000  0.0000000000000000

  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00

POSCAR
) > ./POSCAR

# ... and the KPOINTS file.
(
cat << KPOINTS
elephon example file for MgB2
0
Gamma
7 7 6
KPOINTS
) > ./KPOINTS

echo "Running vasp ..."
../../run_vasp.sh
if [ $? -ne 0 ]; then
	echo "Problem running VASP. Exiting ..."
	exit 1;
fi

echo "Running elephon ldos ..."
echo "f_ldos = ldos.dat" > infile
../../../run_elephon.sh
