// Copyright (C) 2004--2014 Constantine Khroulev, Ed Bueler, Jed Brown, Torsten Albrecht
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include "SSA.hh"
#include "Mask.hh"
#include "basal_resistance.hh"
#include "PISMVars.hh"
#include "pism_options.hh"
#include "flowlaw_factory.hh"
#include "PIO.hh"
#include "enthalpyConverter.hh"

SSA::SSA(IceGrid &g, EnthalpyConverter &e, const PISMConfig &c)
  : ShallowStressBalance(g, e, c)
{
  mask = NULL;
  thickness = NULL;
  tauc = NULL;
  surface = NULL;
  bed = NULL;
  enthalpy = NULL;
  driving_stress_x = NULL;
  driving_stress_y = NULL;
  gl_mask = NULL;

  strength_extension = new SSAStrengthExtension(config);
  allocate();
}

SSA::~SSA() { 
  if (deallocate() != 0) {
    PetscPrintf(grid.com, "FATAL ERROR: SSA de-allocation failed.\n");
    PISMEnd();
  }
  delete strength_extension;
}


//! \brief Initialize a generic regular-grid SSA solver.
PetscErrorCode SSA::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = ShallowStressBalance::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"* Initializing the SSA stress balance...\n"); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com,
                    "  [using the %s flow law]\n", flow_law->name().c_str()); CHKERRQ(ierr);
  
  if (config.get_flag("sub_groundingline")) {
    gl_mask = dynamic_cast<IceModelVec2S*>(vars.get("gl_mask"));
    if (gl_mask == NULL)
      SETERRQ(grid.com, 1, "subgrid_grounding_line_position is not available");
  }

  mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  if (mask == NULL)
    SETERRQ(grid.com, 1, "mask is not available");

  thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (thickness == NULL)
    SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  tauc = dynamic_cast<IceModelVec2S*>(vars.get("tauc"));
  if (tauc == NULL)
    SETERRQ(grid.com, 1, "tauc is not available");

  surface = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  driving_stress_x = dynamic_cast<IceModelVec2S*>(vars.get("ssa_driving_stress_x"));
  driving_stress_y = dynamic_cast<IceModelVec2S*>(vars.get("ssa_driving_stress_y"));
  if( (driving_stress_x==NULL) || (driving_stress_y==NULL) ) {
    if(surface == NULL) {
      SETERRQ(grid.com, 1,
              "neither surface_altitude nor the pair ssa_driving_stress_x/y is available");
    }
 }

  bed = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (bed == NULL)
    SETERRQ(grid.com, 1, "bedrock_altitude is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(vars.get("enthalpy"));
  if (enthalpy == NULL)
    SETERRQ(grid.com, 1, "enthalpy is not available");
  
  // Check if PISM is being initialized from an output file from a previous run
  // and read the initial guess (unless asked not to).
  bool i_set;
  std::string filename;
  ierr = PISMOptionsString("-i", "PISM input file",
                           filename, i_set); CHKERRQ(ierr);

  if (i_set) {
    bool dont_read_initial_guess, u_ssa_found, v_ssa_found;
    unsigned int start;
    PIO nc(grid, "guess_mode");

    ierr = PISMOptionsIsSet("-dontreadSSAvels", dont_read_initial_guess); CHKERRQ(ierr);

    ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
    ierr = nc.inq_var("u_ssa", u_ssa_found); CHKERRQ(ierr); 
    ierr = nc.inq_var("v_ssa", v_ssa_found); CHKERRQ(ierr); 
    ierr = nc.inq_nrecords(start); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr); 
    start -= 1;

    if (u_ssa_found && v_ssa_found &&
        (! dont_read_initial_guess)) {
      ierr = verbPrintf(3,grid.com,"Reading u_ssa and v_ssa...\n"); CHKERRQ(ierr);

      ierr = m_velocity.read(filename, start); CHKERRQ(ierr);
    }

  } else {
    ierr = m_velocity.set(0.0); CHKERRQ(ierr); // default initial guess
  }

  if (config.get_flag("ssa_dirichlet_bc")) {
    bc_locations = dynamic_cast<IceModelVec2Int*>(vars.get("bcflag"));
    if (bc_locations == NULL) SETERRQ(grid.com, 1, "bc_locations is not available");

    m_vel_bc = dynamic_cast<IceModelVec2V*>(vars.get("vel_ssa_bc"));
    if (m_vel_bc == NULL) SETERRQ(grid.com, 1, "vel_ssa_bc is not available");
  }

  return 0;
}

//! \brief Allocate objects which any SSA solver would use.
PetscErrorCode SSA::allocate() {
  PetscErrorCode ierr;

  ierr = taud.create(grid, "taud", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = taud.set_attrs("diagnostic",
                        "X-component of the driving shear stress at the base of ice",
                        "Pa", "", 0); CHKERRQ(ierr);
  ierr = taud.set_attrs("diagnostic",
                        "Y-component of the driving shear stress at the base of ice",
                        "Pa", "", 1); CHKERRQ(ierr);

  ierr = m_velocity_old.create(grid, "velocity_old", WITH_GHOSTS); CHKERRQ(ierr);
  ierr = m_velocity_old.set_attrs("internal",
                                  "old SSA velocity field; used for re-trying with a different epsilon",
                                  "m s-1", ""); CHKERRQ(ierr);

  // override velocity metadata
  std::vector<std::string> long_names;
  long_names.push_back("SSA model ice velocity in the X direction");
  long_names.push_back("SSA model ice velocity in the Y direction");
  ierr = m_velocity.rename("_ssa",long_names,""); CHKERRQ(ierr);

  int dof=2, stencil_width=1;
  ierr = grid.get_dm(dof, stencil_width, SSADA); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(SSADA, &SSAX); CHKERRQ(ierr);

  {
    IceFlowLawFactory ice_factory(grid.com, "ssa_", config, &EC);
    ice_factory.removeType(ICE_GOLDSBY_KOHLSTEDT);

    ierr = ice_factory.setType(config.get_string("ssa_flow_law")); CHKERRQ(ierr);

    ierr = ice_factory.setFromOptions(); CHKERRQ(ierr);
    ierr = ice_factory.create(&flow_law); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode SSA::deallocate() {
  PetscErrorCode ierr;

  if (SSAX != PETSC_NULL) {
    ierr = VecDestroy(&SSAX); CHKERRQ(ierr);
  }

  if (flow_law != NULL) {
    delete flow_law;
    flow_law = NULL;
  }

  return 0;
}


//! \brief Update the SSA solution.
PetscErrorCode SSA::update(bool fast, IceModelVec2S &melange_back_pressure) {
  PetscErrorCode ierr;

  if (fast)
    return 0;

  (void) melange_back_pressure;

  ierr = solve();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISM ERROR: SSA solver failed.\n");
    return ierr;
  }

  ierr = compute_basal_frictional_heating(m_velocity, *tauc, *mask,
                                          basal_frictional_heating); CHKERRQ(ierr);
  
  return 0;
}

//! \brief Compute the gravitational driving stress.
/*!
Computes the gravitational driving stress at the base of the ice:
\f[ \tau_d = - \rho g H \nabla h \f]

If configuration parameter `surface_gradient_method` = `eta` then the surface
gradient \f$\nabla h\f$ is computed by the gradient of the transformed variable
\f$\eta= H^{(2n+2)/n}\f$ (frequently, \f$\eta= H^{8/3}\f$). The idea is that
this quantity is more regular at ice sheet margins, and so we get a better
surface gradient. When the thickness at a grid point is very small (below \c
minThickEtaTransform in the procedure), the formula is slightly modified to
give a lower driving stress. The transformation is not used in floating ice.
 */
PetscErrorCode SSA::compute_driving_stress(IceModelVec2V &result) {
  PetscErrorCode ierr;

  IceModelVec2S &thk = *thickness; // to improve readability (below)

  const double n = flow_law->exponent(), // frequently n = 3
    etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
    invpow  = 1.0 / etapow,  // = 3/8
    dinvpow = (- n - 2.0) / (2.0 * n + 2.0); // = -5/8
  const double minThickEtaTransform = 5.0; // m
  const double dx=grid.dx, dy=grid.dy;

  bool cfbc = config.get_flag("calving_front_stress_boundary_condition");
  bool compute_surf_grad_inward_ssa = config.get_flag("compute_surf_grad_inward_ssa");
  bool use_eta = (config.get_string("surface_gradient_method") == "eta");

  ierr =   surface->begin_access();    CHKERRQ(ierr);
  ierr =       bed->begin_access();  CHKERRQ(ierr);
  ierr =      mask->begin_access();  CHKERRQ(ierr);
  ierr =        thk.begin_access();  CHKERRQ(ierr);

  MaskQuery m(*mask);

  ierr = result.begin_access(); CHKERRQ(ierr);

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const double pressure = EC.getPressureFromDepth(thk(i,j)); // FIXME issue #15
      if (pressure <= 0.0) {
        result(i,j).u = 0.0;
        result(i,j).v = 0.0;
      } else {
        double h_x = 0.0, h_y = 0.0;
        // FIXME: we need to handle grid periodicity correctly.
        if (m.grounded(i,j) && (use_eta == true)) {
                // in grounded case, differentiate eta = H^{8/3} by chain rule
          if (thk(i,j) > 0.0) {
            const double myH = (thk(i,j) < minThickEtaTransform ?
                                     minThickEtaTransform : thk(i,j));
            const double eta = pow(myH, etapow), factor = invpow * pow(eta, dinvpow);
            h_x = factor * (pow(thk(i+1,j),etapow) - pow(thk(i-1,j),etapow)) / (2*dx);
            h_y = factor * (pow(thk(i,j+1),etapow) - pow(thk(i,j-1),etapow)) / (2*dy);
          }
          // now add bed slope to get actual h_x,h_y
          // FIXME: there is no reason to assume user's bed is periodized
          h_x += bed->diff_x(i,j);
          h_y += bed->diff_y(i,j);
        } else {  // floating or eta transformation is not used
          if (compute_surf_grad_inward_ssa) {
            // Special case for verification tests.
            h_x = surface->diff_x_p(i,j);
            h_y = surface->diff_y_p(i,j);
          } else {              // general case

            // To compute the x-derivative we use
            // * away from the grounding line -- 2nd order centered difference
            //
            // * at the grounded cell near the grounding line -- 1st order
            //   one-sided difference using the grounded neighbor
            //
            // * at the floating cell near the grounding line -- 1st order
            //   one-sided difference using the floating neighbor
            //
            // All three cases can be combined by writing h_x as the weighted
            // average of one-sided differences, with weights of 0 if a finite
            // difference is not used and 1 if it is.
            //
            // The y derivative is handled the same way.

            // x-derivative
            {
              double west = 1, east = 1;
              if ((m.grounded(i,j) && m.floating_ice(i+1,j)) || (m.floating_ice(i,j) && m.grounded(i+1,j)) ||
                  (m.floating_ice(i,j) && m.ice_free_ocean(i+1,j)))
                east = 0;
              if ((m.grounded(i,j) && m.floating_ice(i-1,j)) || (m.floating_ice(i,j) && m.grounded(i-1,j)) ||
                  (m.floating_ice(i,j) && m.ice_free_ocean(i-1,j)))
                west = 0;

              // This driving stress computation has to match the calving front
              // stress boundary condition in SSAFD::assemble_rhs().
              if (cfbc) {
                if (m.icy(i,j) && m.ice_free(i+1,j))
                  east = 0;
                if (m.icy(i,j) && m.ice_free(i-1,j))
                  west = 0;
              }

              if (east + west > 0)
                h_x = 1.0 / (west + east) * (west * surface->diff_x_stagE(i-1,j) +
                                             east * surface->diff_x_stagE(i,j));
              else
                h_x = 0.0;
            }

            // y-derivative
            {
              double south = 1, north = 1;
              if ((m.grounded(i,j) && m.floating_ice(i,j+1)) || (m.floating_ice(i,j) && m.grounded(i,j+1)) ||
                  (m.floating_ice(i,j) && m.ice_free_ocean(i,j+1)))
                north = 0;
              if ((m.grounded(i,j) && m.floating_ice(i,j-1)) || (m.floating_ice(i,j) && m.grounded(i,j-1)) ||
                  (m.floating_ice(i,j) && m.ice_free_ocean(i,j-1)))
                south = 0;

              // This driving stress computation has to match the calving front
              // stress boundary condition in SSAFD::assemble_rhs().
              if (cfbc) {
                if (m.icy(i,j) && m.ice_free(i,j+1))
                  north = 0;
                if (m.icy(i,j) && m.ice_free(i,j-1))
                  south = 0;
              }

              if (north + south > 0)
                h_y = 1.0 / (south + north) * (south * surface->diff_y_stagN(i,j-1) +
                                               north * surface->diff_y_stagN(i,j));
              else
                h_y = 0.0;
            }

          } // end of "general case"

        } // end of "floating or eta transformation is not used"

        result(i,j).u = - pressure * h_x;
        result(i,j).v = - pressure * h_y;
      } // end of "(pressure > 0)"
    } // inner loop (j)
  } // outer loop (i)

  ierr =        thk.end_access(); CHKERRQ(ierr);
  ierr =       bed->end_access(); CHKERRQ(ierr);
  ierr =   surface->end_access(); CHKERRQ(ierr);
  ierr =      mask->end_access(); CHKERRQ(ierr);
  ierr =     result.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode SSA::stdout_report(std::string &result) {
  result = stdout_ssa;
  return 0;
}


//! \brief Set the initial guess of the SSA velocity.
PetscErrorCode SSA::set_initial_guess(IceModelVec2V &guess) {
  PetscErrorCode ierr;
  ierr = m_velocity.copy_from(guess); CHKERRQ(ierr);
  return 0;
}


void SSA::add_vars_to_output(std::string /*keyword*/, std::set<std::string> &result) {
  result.insert("vel_ssa");
}


PetscErrorCode SSA::define_variables(std::set<std::string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "vel_ssa")) {
    ierr = m_velocity.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode SSA::write_variables(std::set<std::string> vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, "vel_ssa")) {
    ierr = m_velocity.write(nc); CHKERRQ(ierr);
  }

  return 0;
}

void SSA::get_diagnostics(std::map<std::string, PISMDiagnostic*> &dict,
                          std::map<std::string, PISMTSDiagnostic*> &/*ts_dict*/) {
    dict["taud"] = new SSA_taud(this, grid, *variables);
    dict["taud_mag"] = new SSA_taud_mag(this, grid, *variables);
    dict["taub"] = new SSA_taub(this, grid, *variables);
    dict["beta"] = new SSA_beta(this, grid, *variables);
}

SSA_taud::SSA_taud(SSA *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<SSA>(m, g, my_vars) {

  dof = 2;
  vars.resize(dof,  NCSpatialVariable(g.get_unit_system()));
  // set metadata:
  vars[0].init_2d("taud_x", grid);
  vars[1].init_2d("taud_y", grid);

  set_attrs("X-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  set_attrs("Y-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 1);

  for (int k = 0; k < dof; ++k)
    vars[k].set_string("comment",
                       "this is the driving stress used by the SSA solver");
}

PetscErrorCode SSA_taud::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2V *result = new IceModelVec2V;
  ierr = result->create(grid, "result", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];

  ierr = model->compute_driving_stress(*result); CHKERRQ(ierr);

  output = result;
  return 0;
}

SSA_taud_mag::SSA_taud_mag(SSA *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<SSA>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("taud_mag", grid);

  set_attrs("magnitude of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  vars[0].set_string("comment",
                     "this is the magnitude of the driving stress used by the SSA solver");
}

PetscErrorCode SSA_taud_mag::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  // Allocate memory:
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "taud_mag", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  IceModelVec* tmp;
  SSA_taud diag(model, grid, variables);

  ierr = diag.compute(tmp);

  IceModelVec2V *taud = dynamic_cast<IceModelVec2V*>(tmp);
  if (taud == NULL)
    SETERRQ(grid.com, 1, "expected an IceModelVec2V, but dynamic_cast failed");

  ierr = taud->magnitude(*result); CHKERRQ(ierr);

  delete tmp;

  output = result;
  return 0;
}

SSA_taub::SSA_taub(SSA *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<SSA>(m, g, my_vars) {
  dof = 2;
  vars.resize(dof, NCSpatialVariable(g.get_unit_system()));
  // set metadata:
  vars[0].init_2d("taub_x", grid);
  vars[1].init_2d("taub_y", grid);

  set_attrs("X-component of the shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  set_attrs("Y-component of the shear stress at the base of ice", "",
            "Pa", "Pa", 1);

  for (int k = 0; k < dof; ++k)
    vars[k].set_string("comment",
                       "this field is purely diagnostic (not used by the model)");
}

PetscErrorCode SSA_taub::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2V *result = new IceModelVec2V;
  ierr = result->create(grid, "result", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];

  IceModelVec2V *velocity;
  ierr = model->get_2D_advective_velocity(velocity); CHKERRQ(ierr);

  IceModelVec2V &vel = *velocity;
  IceModelVec2S &tauc = *model->tauc;
  IceModelVec2Int &mask = *model->mask;

  MaskQuery m(mask);

  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = tauc.begin_access(); CHKERRQ(ierr);
  ierr = vel.begin_access(); CHKERRQ(ierr);
  ierr = mask.begin_access(); CHKERRQ(ierr);
  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (m.grounded_ice(i,j)) {
        double beta = model->basal_sliding_law->drag(tauc(i,j), vel(i,j).u, vel(i,j).v);
        (*result)(i,j).u = - beta * vel(i,j).u;
        (*result)(i,j).v = - beta * vel(i,j).v;
      } else {
        (*result)(i,j).u = 0.0;
        (*result)(i,j).v = 0.0;
      }
    }
  }
  ierr = mask.end_access(); CHKERRQ(ierr);
  ierr = vel.end_access(); CHKERRQ(ierr);
  ierr = tauc.end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);

  output = result;

  return 0;
}

SSA_beta::SSA_beta(SSA *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<SSA>(m, g, my_vars) {
  // set metadata:
  vars[0].init_2d("beta", grid);

  set_attrs("basal drag coefficient", "", "Pa s / m", "Pa s / m", 0);
}

PetscErrorCode SSA_beta::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  // Allocate memory:
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "beta", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  IceModelVec2S *tauc = model->tauc;
  IceModelVec2V *velocity;
  ierr = model->get_2D_advective_velocity(velocity); CHKERRQ(ierr);

  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = tauc->begin_access(); CHKERRQ(ierr);
  ierr = velocity->begin_access(); CHKERRQ(ierr);
  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      (*result)(i,j) = model->basal_sliding_law->drag((*tauc)(i,j),
                                                      (*velocity)(i,j).u, (*velocity)(i,j).v);
    }
  }
  ierr = velocity->end_access(); CHKERRQ(ierr);
  ierr = tauc->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}
