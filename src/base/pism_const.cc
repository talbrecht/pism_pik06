// Copyright (C) 2007--2008 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <petsc.h>
#include <petscfix.h>
#include "nc_util.hh"
#include "pism_const.hh"

PetscErrorCode getFlowLawNumber(PetscInt &flowLawNum, const PetscInt defaultFLN) {
    PetscErrorCode ierr;
    PetscTruth     flowlawSet = PETSC_FALSE, useGK = PETSC_FALSE, useGKAGE = PETSC_FALSE;

    flowLawNum = defaultFLN;

    ierr = PetscOptionsGetInt(PETSC_NULL, "-law", &flowLawNum, &flowlawSet); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-gk_age", &useGKAGE); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-gk", &useGK); CHKERRQ(ierr);
    if ((useGKAGE == PETSC_TRUE) || (useGK == PETSC_TRUE)) {
      flowLawNum = 4;
    }
    return 0;
}


PetscErrorCode userChoosesIceType(MPI_Comm com, IceType* &ice) {
    PetscErrorCode ierr;
    ierr = userChoosesIceType(com, ice, 0); CHKERRQ(ierr); // use Paterson-Budd by default
    return 0;
}


PetscErrorCode userChoosesIceType(MPI_Comm com, IceType* &ice, const PetscInt defaultFLN) {
    PetscErrorCode ierr;
    PetscInt       myflowLawNum;
    
    ierr = getFlowLawNumber(myflowLawNum, defaultFLN); CHKERRQ(ierr);

    ierr = verbPrintf(3,com, 
        "  [using flow law %d (where 0=Paterson-Budd,1=cold P-B,2=warm P-B,"
        "3=Hooke,4=Goldsby-Kohlstedt)]\n",
        myflowLawNum); CHKERRQ(ierr);
    
    switch (myflowLawNum) {
      case 0: // Paterson-Budd
        ice = new ThermoGlenIce;  
        break;
      case 1: // cold part of P-B
        ice = new ThermoGlenArrIce;  
        break;
      case 2: // warm part of P-B
        ice = new ThermoGlenArrIceWarm;  
        break;
      case 3: // Hooke
        ice = new ThermoGlenIceHooke;
        break;
      case 4: // Goldsby Kohlstedt
        ice = new HybridIce;  
        break;
      case 5: // Goldsby Kohlstedt stripped down
        ice = new HybridIceStripped;  
        break;
      default:
	{
	  PetscPrintf(com, "\n PISM ERROR: flow law number must be 0,1,2,3,4\n");
	  PetscEnd();
	}
    }
    return 0;
}



// Verbosity level version of PetscPrintf.  We print according to whether 
// (thresh <= verbosityLevel), in which case print, or 
// (thresh > verbosityLevel) in which case no print.
//
//   level  option        meaning
//   -----  ------        -------
//   0      -verbose 0    never print to std out AT ALL!
//
//   1      -verbose 1    less verbose than default: thresh must be 1 to print
//
//   2     [-verbose 2]   default
//
//   3      -verbose      somewhat verbose
//         [-verbose 3]   
//   4      -vverbose     fairly verbose
//         [-verbose 4]
//   5      -vvverbose    very verbose: if level this high then (thresh <= level) 
//         [-verbose 5]       always, so print everything
//
// note: 1 <= thresh <= 5  enforced in verbPrintf() below

static PetscInt verbosityLevel;

PetscErrorCode setVerbosityLevel(PetscInt level) {
  if ((level < 0) || (level > 5)) {
     SETERRQ(1,"verbosity level invalid");
  }
  verbosityLevel = level;
  return 0;  
}

PetscInt getVerbosityLevel() {
  return verbosityLevel;
}


PetscErrorCode verbosityLevelFromOptions() {
  // verbosity options: more info to standard out
  PetscErrorCode ierr;
  PetscInt     myverbosityLevel;
  PetscTruth   verbose, verbosityLevelSet;
  PetscInt     DEFAULT_VERBOSITY_LEVEL = 2;
  
  ierr = setVerbosityLevel(DEFAULT_VERBOSITY_LEVEL);  
  ierr = PetscOptionsGetInt(PETSC_NULL, "-verbose", &myverbosityLevel, &verbosityLevelSet); CHKERRQ(ierr);
  if (verbosityLevelSet == PETSC_TRUE) {
    ierr = setVerbosityLevel(myverbosityLevel);
  } else {
    ierr = PetscOptionsHasName(PETSC_NULL, "-verbose", &verbose); CHKERRQ(ierr);
    if (verbose == PETSC_TRUE)   ierr = setVerbosityLevel(3);
  }
  ierr = PetscOptionsHasName(PETSC_NULL, "-vverbose", &verbose); CHKERRQ(ierr);
  if (verbose == PETSC_TRUE)   ierr = setVerbosityLevel(4);
  ierr = PetscOptionsHasName(PETSC_NULL, "-vvverbose", &verbose); CHKERRQ(ierr);
  if (verbose == PETSC_TRUE)   ierr = setVerbosityLevel(5);
  return 0;
}


PetscErrorCode verbPrintf(const int thresh, 
                          MPI_Comm comm,const char format[],...)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  size_t         len;
  char           *buffer,*sub1,*sub2;
  const char     *nformat;
  PetscReal      value;

  extern FILE *petsc_history;

  if ((thresh < 1) || (thresh > 5)) { SETERRQ(1,"invalid threshold in verbPrintf()"); }

  PetscFunctionBegin;
  if (!comm) comm = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  if (!rank && ((verbosityLevel >= thresh) || petsc_history) ) {
    va_list Argp;
    va_start(Argp,format);

    ierr = PetscStrstr(format,"%A",&sub1);CHKERRQ(ierr);
    if (sub1) {
      ierr = PetscStrstr(format,"%",&sub2);CHKERRQ(ierr);
      if (sub1 != sub2) SETERRQ(PETSC_ERR_ARG_WRONG,"%%A format must be first in format string");
      ierr    = PetscStrlen(format,&len);CHKERRQ(ierr);
      ierr    = PetscMalloc((len+16)*sizeof(char),&buffer);CHKERRQ(ierr);
      ierr    = PetscStrcpy(buffer,format);CHKERRQ(ierr);
      ierr    = PetscStrstr(buffer,"%",&sub2);CHKERRQ(ierr);
      sub2[0] = 0;
      value   = (double)va_arg(Argp,double);
      if (PetscAbsReal(value) < 1.e-12) {
        ierr    = PetscStrcat(buffer,"< 1.e-12");CHKERRQ(ierr);
      } else {
        ierr    = PetscStrcat(buffer,"%g");CHKERRQ(ierr);
        va_end(Argp);
        va_start(Argp,format);
      }
      ierr    = PetscStrcat(buffer,sub1+2);CHKERRQ(ierr);
      nformat = buffer;
    } else {
      nformat = format;
    }
    if (verbosityLevel >= thresh) {
      ierr = PetscVFPrintf(PETSC_STDOUT,nformat,Argp);CHKERRQ(ierr);
    }
    if (petsc_history) { // always print to history
      ierr = PetscVFPrintf(petsc_history,nformat,Argp);CHKERRQ(ierr);
    }
    va_end(Argp);
    if (sub1) {ierr = PetscFree(buffer);CHKERRQ(ierr);}
  }
  PetscFunctionReturn(0);
}

//! Comparison function from the glibc manual
//! http://www.gnu.org/software/libc/manual/index.html
int compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  
  return (*da > *db) - (*da < *db);
}
