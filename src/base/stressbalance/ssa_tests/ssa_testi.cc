// Copyright (C) 2010--2011 Ed Bueler, Constantine Khroulev, and David Maxwell
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

static char help[] =
  "\nSSAFD_TEST\n"
  "  Testing program for the finite element implementation of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof. Uses verification test I. Also may be used in a PISM\n"
  "  software (regression) test.\n\n";

#include "pism_const.hh"
#include "iceModelVec.hh"
#include "flowlaws.hh" // IceFlowLaw
#include "materials.hh" // IceBasalResistancePlasticLaw
#include "PISMIO.hh"
#include "NCVariable.hh"
#include "SSAFEM.hh"
#include "SSAFD.hh"
#include "exactTestsIJ.h"
#include "SSATestCase.hh"

class SSATestCaseI: public SSATestCase
{
public:
  SSATestCaseI( MPI_Comm com, PetscMPIInt rank, 
                 PetscMPIInt size, NCConfigVariable &c ): 
                 SSATestCase(com,rank,size,c)
  { };
  
protected:
  virtual PetscErrorCode initializeGrid(PetscInt Mx,PetscInt My);

  virtual PetscErrorCode initializeSSAModel();

  virtual PetscErrorCode initializeSSACoefficients();

  virtual PetscErrorCode exactSolution(PetscInt i, PetscInt j, 
    PetscReal x, PetscReal y, PetscReal *u, PetscReal *v );

};

const PetscScalar m_schoof = 10; // (pure number)
const PetscScalar L_schoof = 40e3; // meters
const PetscScalar aspect_schoof = 0.05; // (pure)
const PetscScalar H0_schoof = aspect_schoof * L_schoof; 
                                       // = 2000 m THICKNESS
const PetscScalar B_schoof = 3.7e8; // Pa s^{1/3}; hardness 
                                     // given on p. 239 of Schoof; why so big?
const PetscScalar p_schoof = 4.0/3.0; // = 1 + 1/n


PetscErrorCode SSATestCaseI::initializeGrid(PetscInt Mx,PetscInt My)
{
  PetscReal Ly = 3*L_schoof;  // 300.0 km half-width (L=40.0km in Schoof's choice of variables)
  PetscReal Lx = PetscMax(60.0e3, ((Mx - 1) / 2) * (2.0 * Ly / (My - 1)) );
  init_shallow_grid(grid,Lx,Ly,Mx,My,NONE);
  return 0;
}


PetscErrorCode SSATestCaseI::initializeSSAModel()
{
  basal = new IceBasalResistancePlasticLaw(
         config.get("plastic_regularization") / secpera,
         config.get_flag("do_pseudo_plastic_till"),
         config.get("pseudo_plastic_q"),
         config.get("pseudo_plastic_uthreshold") / secpera);

  CustomGlenIce *glenIce = new CustomGlenIce(grid.com, "", config);
  glenIce->setHardness(B_schoof);
  ice = glenIce;

  enthalpyconverter = new EnthalpyConverter(config);
  return 0;
}

PetscErrorCode SSATestCaseI::initializeSSACoefficients()
{
  PetscErrorCode ierr;
  PetscScalar    **ph, **pbed;
  
  ierr = mask.set(MASK_DRAGGING_SHEET); CHKERRQ(ierr);
  ierr = thickness.set(H0_schoof); CHKERRQ(ierr);

  // ssa->strength_extension->set_min_thickness(2*H0_schoof);
  
  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  config.set_flag("compute_surf_grad_inward_ssa", true);
  // config.set("epsilon_ssa", 0.0);  // don't use this lower bound

  PetscScalar **ptauc;

  ierr = tauc.get_array(ptauc); CHKERRQ(ierr);

  PetscScalar standard_gravity = config.get("standard_gravity");

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar y = grid.y[j];
      const PetscScalar theta = atan(0.001);   /* a slope of 1/1000, a la Siple streams */
      const PetscScalar f = ice->rho * standard_gravity * H0_schoof * tan(theta);
      ptauc[i][j] = f * pow(PetscAbs(y / L_schoof), m_schoof);
    }
  }
  ierr = tauc.end_access(); CHKERRQ(ierr);
  ierr = tauc.beginGhostComm(); CHKERRQ(ierr);
  ierr = tauc.endGhostComm(); CHKERRQ(ierr);


  ierr = vel_bc.begin_access(); CHKERRQ(ierr);
  ierr = mask.begin_access(); CHKERRQ(ierr);
  ierr = surface.get_array(ph); CHKERRQ(ierr);    
  ierr = bed.get_array(pbed); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar junk, myu, myv;
      const PetscScalar myx = grid.x[i], myy=grid.y[j];
      // eval exact solution; will only use exact vels if at edge
      exactI(m_schoof, myx, myy, &(pbed[i][j]), &junk, &myu, &myv); 
      ph[i][j] = pbed[i][j] + H0_schoof;
      
      bool edge = ( (j == 0) || (j == grid.My - 1) ) || ( (i==0) || (i==grid.Mx-1) );
      if (edge) {
        mask(i,j) = MASK_SHEET;
        vel_bc(i,j).u = myu;
        vel_bc(i,j).v = myv;
      }
    }
  } 
  ierr = vel_bc.end_access(); CHKERRQ(ierr);
  ierr = mask.end_access(); CHKERRQ(ierr);
  ierr = surface.end_access(); CHKERRQ(ierr);    
  ierr = bed.end_access(); CHKERRQ(ierr);

  // communicate what we have set
  ierr = surface.beginGhostComm(); CHKERRQ(ierr);
  ierr = surface.endGhostComm(); CHKERRQ(ierr);
  ierr = bed.beginGhostComm(); CHKERRQ(ierr);
  ierr = bed.endGhostComm(); CHKERRQ(ierr);
  ierr = mask.beginGhostComm(); CHKERRQ(ierr);
  ierr = mask.endGhostComm(); CHKERRQ(ierr);
  ierr = vel_bc.beginGhostComm(); CHKERRQ(ierr);
  ierr = vel_bc.endGhostComm(); CHKERRQ(ierr);

  ierr = ssa->set_boundary_conditions(mask, vel_bc); CHKERRQ(ierr); 

  return 0;
}


PetscErrorCode SSATestCaseI::exactSolution(PetscInt /*i*/, PetscInt /*j*/, 
                                           PetscReal x, PetscReal y,
                                           PetscReal *u, PetscReal *v)
{
  PetscReal junk1, junk2;
  exactI(m_schoof, x,y, &junk1, &junk2,u,v); 
  return 0;
}


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;  // won't be used except for rank,size
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);
  
  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {  
    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    ierr = setVerbosityLevel(5); CHKERRQ(ierr);

    PetscTruth usage_set, help_set;
    ierr = PetscOptionsHasName(PETSC_NULL, "-usage", &usage_set); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-help", &help_set); CHKERRQ(ierr);
    if ((usage_set==PETSC_TRUE) || (help_set==PETSC_TRUE)) {
      PetscPrintf(com,
                  "\n"
                  "usage of SSA_TESTJ:\n"
                  "  run ssafe_test -Mx <number> -My <number> -ssa <fd|fem>\n"
                  "\n");
    }

    // Parameters that can be overridden by command line options
    PetscInt Mx=11;
    PetscInt My=61;
    string output_file = "ssa_test_i.nc";
    string driver = "fem";

    ierr = PetscOptionsBegin(com, "", "SSAFD_TEST options", ""); CHKERRQ(ierr);
    {
      bool flag;
      ierr = PISMOptionsInt("-Mx", "Number of grid points in the X direction", 
                                                      Mx, flag); CHKERRQ(ierr);
      ierr = PISMOptionsInt("-My", "Number of grid points in the Y direction", 
                                                      My, flag); CHKERRQ(ierr);
      ierr = PISMOptionsString("-ssa_method", "Algorithm for computing the SSA solution",
                                                  driver, flag); CHKERRQ(ierr);                                                      
      ierr = PISMOptionsString("-o", "Set the output file name", 
                                              output_file, flag); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    // Determine the kind of solver to use.
    SSAFactory ssafactory;
    if(driver.compare("fem") == 0) ssafactory = SSAFEMFactory;
    else if(driver.compare("fd") == 0) ssafactory = SSAFDFactory;
    else SETERRQ(1,"SSA algorithm argument should be one of -ssa fe or -ssa fem");

    SSATestCaseI testcase(com,rank,size,config);
    ierr = testcase.init(Mx,My,ssafactory); CHKERRQ(ierr);
    ierr = testcase.run(); CHKERRQ(ierr);
    ierr = testcase.report(); CHKERRQ(ierr);
    ierr = testcase.write(output_file); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}