#!/bin/bash

# Copyright (C) 2009-2013 The PISM Authors

# PISM Greenland spinup using either constant present-day climate or modeled
# paleoclimate.  See README.md.

# Before using this script, run preprocess.sh to download and adjust metadata
# in the SeaRISE "Present Day Greenland" master dataset.

set -e  # exit on error

GRIDLIST="{20, 10, 5, 3, 2}"
CLIMLIST="{const, paleo}"
DYNALIST="{sia, hybrid}"

# preprocess.sh generates pism_*.nc files; run it first
PISM_DATAVERSION=1.1
PISM_DATANAME=pism_Greenland_5km_v$PISM_DATAVERSION.nc
PISM_TEMPSERIES=pism_dT.nc
PISM_SLSERIES=pism_dSL.nc

if [ $# -lt 5 ] ; then
  echo "spinup.sh ERROR: needs 5 or 6 or 7 positional arguments ... ENDING NOW"
  echo
  echo "usage:"
  echo
  echo "    spinup.sh PROCS CLIMATE DURATION GRID DYNAMICS [OUTFILE] [BOOTFILE]"
  echo
  echo "  where:"
  echo "    PROCS     = 1,2,3,... is number of MPI processes"
  echo "    CLIMATE   in $CLIMLIST"
  echo "    DURATION  = model run time in years; does '-ys -DURATION -ye 0'"
  echo "    GRID      in $GRIDLIST (km)"
  echo "    DYNAMICS  in $DYNALIST; sia is non-sliding"
  echo "    OUTFILE   optional name of output file; default = unnamed.nc"
  echo "    BOOTFILE  optional name of input file; default = $PISM_DATANAME"
  echo
  echo "consider setting optional environment variables (see script for meaning):"
  echo "    EXSTEP       spacing in years between -extra_files outputs; defaults to 100"
  echo "    PISM_DO      set to 'echo' if no run desired; defaults to empty"
  echo "    PISM_MPIDO   defaults to 'mpiexec -n'"
  echo "    PISM_PREFIX  set to path to pismr executable if desired; defaults to empty"
  echo "    PISM_EXEC    defaults to 'pismr'"
  echo "    REGRIDFILE   set to file name to regrid from; defaults to empty (no regrid)"
  echo "    REGRIDVARS   set to vars to regrid *if* REGRIDFILE set;"
  echo "                 defaults to 'bmelt,enthalpy,litho_temp,thk,tillwat'"
  echo
  echo "example usage 1:"
  echo "    $ ./spinup.sh 4 const 1000 20 sia"
  echo "  does spinup with 4 processors, constant-climate, 1000 year run, 20 km"
  echo "  grid, and non-sliding SIA stress balance, and which bootstraps and"
  echo "  outputs to default files"
  echo
  echo "example usage 2:"
  echo "    $ PISM_DO=echo ./spinup.sh 128 paleo 100.0 5 hybrid out.nc boot.nc &> foo.sh"
  echo "  creates a script foo.sh that will do spinup with 128 processors,"
  echo "  simulated paleo-climate, '-ys -100 -ye 0' run, 5 km grid, sliding with"
  echo "  SSA stress balance, outputs to out.nc, and which bootstraps from NetCDF"
  echo "  file boot.nc (which might come from previous PISM run, even at a"
  echo "  different grid spacing)"
  echo
  exit
fi

if [ -n "${SCRIPTNAME:+1}" ] ; then
  echo "[SCRIPTNAME=$SCRIPTNAME (already set)]"
  echo ""
else
  SCRIPTNAME="#(spinup.sh)"
fi

if [ $# -gt 7 ] ; then
  echo "$SCRIPTNAME WARNING: ignoring arguments after argument 7 ..."
fi

NN="$1" # first arg is number of processes
DURATION=$3
RUNSTARTEND="-ys -$DURATION -ye 0"

# set coupler from argument 2
if [ "$2" = "const" ]; then
  climname="constant-climate"
  INLIST="$PISM_DATANAME"
  COUPLER="-surface given -surface_given_file $PISM_DATANAME"
elif [ "$2" = "paleo" ]; then
  climname="paleo-climate"
  INLIST="$PISM_DATANAME $PISM_TEMPSERIES $PISM_SLSERIES"
  COUPLER=" -bed_def lc -atmosphere searise_greenland,delta_T,paleo_precip -surface pdd -atmosphere_paleo_precip_file $PISM_TEMPSERIES -atmosphere_delta_T_file $PISM_TEMPSERIES -ocean constant,delta_SL -ocean_delta_SL_file $PISM_SLSERIES"
else
  echo "invalid second argument; must be in $CLIMLIST"
  exit
fi

# decide on grid and skip from argument 4
COARSESKIP=10
FINESKIP=20
FINESTSKIP=50
VDIMS="-Lz 4000 -Lbz 2000 -skip -skip_max "
COARSEVGRID="-Mz 101 -Mbz 11 -z_spacing equal ${VDIMS} ${COARSESKIP}"
FINEVGRID="-Mz 201 -Mbz 21 -z_spacing equal ${VDIMS} ${FINESKIP}"
FINESTVGRID="-Mz 401 -Mbz 41 -z_spacing equal ${VDIMS} ${FINESTSKIP}"
if [ "$4" -eq "20" ]; then
  dx=20
  myMx=76
  myMy=141
  vgrid=$COARSEVGRID
elif [ "$4" -eq "10" ]; then
  dx=10
  myMx=151
  myMy=281
  vgrid=$FINEVGRID
elif [ "$4" -eq "5" ]; then
  # "native" resolution in data file, with 561 x 301 grid
  dx=5
  myMx=301
  myMy=561
  vgrid=$FINEVGRID
elif [ "$4" -eq "3" ]; then
  dx=3
  myMx=501
  myMy=934
  vgrid=$FINESTVGRID
elif [ "$4" -eq "2" ]; then
  dx=2
  myMx=750
  myMy=1400
  vgrid=$FINESTVGRID
else
  echo "invalid fourth argument: must be in $GRIDLIST"
  exit
fi

# set stress balance from argument 5
PHYS="-ocean_kill $PISM_DATANAME -sia_e 3.0"
if [ "$5" = "sia" ]; then
  PHYS=$PHYS
elif [ "$5" = "hybrid" ]; then
  TILLPHI="-topg_to_phi 15.0,40.0,-300.0,700.0"
  PHYS="${PHYS} -ssa_sliding ${TILLPHI} -pseudo_plastic -pseudo_plastic_q 0.23 -till_effective_fraction_overburden 0.01 -tauc_slippery_grounding_lines"
else
  echo "invalid fifth argument; must be in $DYNALIST"
  exit
fi

# set output filename from argument 6
if [ -z "$6" ]; then
  OUTNAME=unnamed.nc
else
  OUTNAME=$6
fi

# set bootstrapping input filename from argument 6
if [ -z "$7" ]; then
  INNAME=$PISM_DATANAME
else
  INNAME=$7
fi

# now we have read options ... we know enough to report to user ...
echo
echo "# ======================================================================="
echo "# PISM std Greenland spinup:"
echo "#    $NN processors, $DURATION a run, $dx km grid, $climname, $5 dynamics"
echo "# ======================================================================="

echo "$SCRIPTNAME              NN = $NN"

# set MPIDO if using different MPI execution command, for example:
#  $ export PISM_MPIDO="aprun -n "
if [ -n "${PISM_MPIDO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME      PISM_MPIDO = $PISM_MPIDO  (already set)"
else
  PISM_MPIDO="mpiexec -n "
  echo "$SCRIPTNAME      PISM_MPIDO = $PISM_MPIDO"
fi

# check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var DO is already set
  echo "$SCRIPTNAME         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO="" 
fi

# prefix to pism (not to executables)
if [ -n "${PISM_PREFIX:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX  (already set)"
else
  PISM_PREFIX=""    # just a guess
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX"
fi

# set PISM_EXEC if using different executables, for example:
#  $ export PISM_EXEC="pismr -cold"
if [ -n "${PISM_EXEC:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC  (already set)"
else
  PISM_EXEC="pismr"
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC"
fi

# set EXSTEP to default if not set
if [ -n "${EXSTEP:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME          EXSTEP = $EXSTEP  (already set)"
else
  EXSTEP="100"
  echo "$SCRIPTNAME          EXSTEP = $EXSTEP"
fi

# if REGRIDFILE set then form regridcommand
if [ -n "${REGRIDFILE:+1}" ] ; then  # check if env var is already set
  if [ -n "${REGRIDVARS:+1}" ] ; then  # check if env var is already set
    echo "$SCRIPTNAME      REGRIDVARS = $REGRIDVARS  (already set)"
  else
    REGRIDVARS='litho_temp,thk,enthalpy,tillwat,bmelt'
  fi
  # note: other vars which are "state":  Href, dbdt, shelfbtemp, shelfbmassflux
  echo "$SCRIPTNAME      REGRIDFILE = $REGRIDFILE; we will regrid these vars:"
  echo "$SCRIPTNAME                     $REGRIDVARS"
  regridcommand="-regrid_file $REGRIDFILE -regrid_vars $REGRIDVARS"
  if [ "$2" = "paleo" ]; then
    regridcommand="$regridcommand -regrid_bed_special"
  fi
else
  regridcommand=""
fi

# actually check for input files
for INPUT in $INLIST $INNAME; do
  if [ -e "$INPUT" ] ; then  # check if file exist
    echo "$SCRIPTNAME           input   $INPUT (found)"
  else
    echo "$SCRIPTNAME           input   $INPUT (MISSING!!)"
    echo
    echo "$SCRIPTNAME           please run ./preprocess.sh, exiting"
    echo
    exit
  fi
done

# show remaining setup options:
PISM="${PISM_PREFIX}${PISM_EXEC}"
echo "$SCRIPTNAME      executable = '$PISM'"
echo "$SCRIPTNAME         coupler = '$COUPLER'"
echo "$SCRIPTNAME        dynamics = '$PHYS'"

# set up diagnostics
TSNAME=ts_$OUTNAME
TSTIMES=-$DURATION:yearly:0
EXNAME=ex_$OUTNAME
EXTIMES=-$DURATION:$EXSTEP:0
# check_stationarity.py can be applied to $EXNAME
EXVARS="diffusivity,temppabase,tempicethk_basal,bmelt,tillwat,csurf,mask,thk,topg,usurf"
if [ "$5" = "hybrid" ]; then
  EXVARS="${EXVARS},hardav,cbase,tauc"
fi
DIAGNOSTICS="-ts_file $TSNAME -ts_times $TSTIMES -extra_file $EXNAME -extra_times $EXTIMES -extra_vars $EXVARS"

# construct command
cmd="$PISM_MPIDO $NN $PISM -boot_file $INNAME -Mx $myMx -My $myMy $vgrid $RUNSTARTEND $regridcommand $COUPLER $PHYS $DIAGNOSTICS -o $OUTNAME"
echo
$PISM_DO $cmd
