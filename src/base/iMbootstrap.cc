// Copyright (C) 2004-2014 Jed Brown, Nathan Shemonski, Ed Bueler and
// Constantine Khroulev
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

#include "iceModel.hh"
#include "PIO.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "enthalpyConverter.hh"
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include <cmath>                // for erf() in method 1 in putTempAtDepth()
#include <assert.h>

//! Read file and use heuristics to initialize PISM from typical 2d data available through remote sensing.
/*! 
This procedure is called by the base class when option `-boot_file` is used.

See chapter 4 of the User's Manual.  We read only 2D information from the bootstrap file.
 */
PetscErrorCode IceModel::bootstrapFromFile(std::string filename) {
  PetscErrorCode  ierr;

  // Bootstrap 2D fields:
  ierr = bootstrap_2d(filename); CHKERRQ(ierr);

  // Regrid 2D fields:
  ierr = regrid(2); CHKERRQ(ierr);

  // Check the consistency of geometry fields:
  ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr); 

  ierr = verbPrintf(2, grid.com,
                    "getting surface B.C. from couplers...\n"); CHKERRQ(ierr);

  // Update couplers (because heuristics in bootstrap_3d() might need boundary
  // conditions provided by couplers):
  ierr = init_step_couplers(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "bootstrapping 3D variables...\n"); CHKERRQ(ierr);

  // Fill 3D fields using heuristics:
  ierr = bootstrap_3d(); CHKERRQ(ierr);

  ierr = regrid(3); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, "done reading %s; bootstrapping done\n",filename.c_str()); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::bootstrap_2d(std::string filename) {
  PetscErrorCode ierr;

  PIO nc(grid, "guess_mode");
  ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, 
                    "bootstrapping by PISM default method from file %s\n", filename.c_str()); CHKERRQ(ierr);

  // report on resulting computational box, rescale grid, actually create a
  // local interpolation context
  ierr = verbPrintf(2, grid.com, 
                    "  rescaling computational box for ice from -boot_file file and\n"
                    "    user options to dimensions:\n"
                    "    [%6.2f km, %6.2f km] x [%6.2f km, %6.2f km] x [0 m, %6.2f m]\n",
                    (grid.x0 - grid.Lx)/1000.0,
                    (grid.x0 + grid.Lx)/1000.0,
                    (grid.y0 - grid.Ly)/1000.0,
                    (grid.y0 + grid.Ly)/1000.0,
                    grid.Lz); 
  CHKERRQ(ierr);

  std::string usurf_name;
  bool hExists=false, maskExists=false, usurf_found_by_std_name;
  ierr = nc.inq_var("usurf", "surface_altitude", hExists, usurf_name, usurf_found_by_std_name); CHKERRQ(ierr);
  ierr = nc.inq_var("mask", maskExists); CHKERRQ(ierr);

  std::string lon_name, lat_name;
  bool lonExists=false, latExists=false, lon_found_by_std_name, lat_found_by_std_name;
  ierr = nc.inq_var("lon", "longitude", lonExists, lon_name, lon_found_by_std_name); CHKERRQ(ierr);
  ierr = nc.inq_var("lat", "latitude",  latExists, lat_name, lat_found_by_std_name); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  // now work through all the 2d variables, regridding if present and otherwise
  // setting to default values appropriately

  if (maskExists) {
    ierr = verbPrintf(2, grid.com, 
                      "  WARNING: 'mask' found; IGNORING IT!\n"); CHKERRQ(ierr);
  }

  if (hExists) {
    ierr = verbPrintf(2, grid.com, 
                      "  WARNING: surface elevation 'usurf' found; IGNORING IT!\n");
                      CHKERRQ(ierr);
  }

  ierr = verbPrintf(2, grid.com, 
                    "  reading 2D model state variables by regridding ...\n"); CHKERRQ(ierr);

  ierr = vLongitude.regrid(filename, OPTIONAL); CHKERRQ(ierr);
  if (!lonExists) {
    vLongitude.metadata().set_string("missing_at_bootstrap","true");
  }

  ierr =  vLatitude.regrid(filename, OPTIONAL); CHKERRQ(ierr);
  if (!latExists) {
    vLatitude.metadata().set_string("missing_at_bootstrap","true");
  }

  ierr = bed_topography.regrid(filename, OPTIONAL, config.get("bootstrapping_bed_value_no_var")); CHKERRQ(ierr);
  ierr = basal_melt_rate.regrid(filename, OPTIONAL, config.get("bootstrapping_bmelt_value_no_var")); CHKERRQ(ierr);
  ierr = geothermal_flux.regrid(filename, OPTIONAL, config.get("bootstrapping_geothermal_flux_value_no_var")); CHKERRQ(ierr);
  ierr = bed_uplift_rate.regrid(filename, OPTIONAL, config.get("bootstrapping_uplift_value_no_var")); CHKERRQ(ierr);


  if (config.get_flag("drainageBasins")) {
    ierr =    vBasinMask.regrid(filename, OPTIONAL, config.get("bootstrapping_basinmask_value_no_var")); CHKERRQ(ierr);
  }

  ierr = ice_thickness.regrid(filename, OPTIONAL, config.get("bootstrapping_H_value_no_var")); CHKERRQ(ierr);
  // check the range of the ice thickness
  {
    double thk_min = 0.0, thk_max = 0.0;
    ierr = ice_thickness.range(thk_min, thk_max); CHKERRQ(ierr);

    if (thk_max >= grid.Lz + 1e-6) {
      PetscPrintf(grid.com,
                  "PISM ERROR: Maximum ice thickness (%f meters)\n"
                  "            exceeds the height of the computational domain (%f meters).\n"
                  "            Stopping...\n", thk_max, grid.Lz);
      PISMEnd();
    }
  }

  if (config.get_flag("part_grid")) {
    // Read the Href field from an input file. This field is
    // grid-dependent, so interpolating it from one grid to a
    // different one does not make sense in general.
    // (IceModel::Href_cleanup() will take care of the side effects of
    // such interpolation, though.)
    //
    // On the other hand, we need to read it in to be able to re-start
    // from a PISM output file using the -boot_file option.
    ierr = vHref.regrid(filename, OPTIONAL, 0.0); CHKERRQ(ierr);
  }

  if (config.get_string("calving_methods").find("eigen_calving") != std::string::npos) {
    ierr = strain_rates.set(0.0); CHKERRQ(ierr);
  }

  if (config.get_flag("ssa_dirichlet_bc")) {
    // Do not use Dirichlet B.C. anywhere if bcflag is not present.
    ierr = vBCMask.regrid(filename, OPTIONAL, 0.0); CHKERRQ(ierr);
    // In the absence of u_ssa_bc and v_ssa_bc in the file the only B.C. that
    // makes sense is the zero Dirichlet B.C.
    ierr = vBCvel.regrid(filename, OPTIONAL,  0.0); CHKERRQ(ierr);
  }

  // check if Lz is valid
  double thk_min, thk_max;
  ierr = ice_thickness.range(thk_min, thk_max); CHKERRQ(ierr);

  if (thk_max > grid.Lz) {
    PetscPrintf(grid.com,
                "PISM ERROR: Max. ice thickness (%3.3f m) exceeds the height of the computational domain (%3.3f m).\n"
                "            Exiting...\n", thk_max, grid.Lz);
    PISMEnd();
  }

  return 0;
}

PetscErrorCode IceModel::bootstrap_3d() {
  PetscErrorCode ierr;

  // set the initial age of the ice if appropriate
  if (config.get_flag("do_age")) {
    ierr = verbPrintf(2, grid.com, 
      "  setting initial age to %.4f years\n", config.get("initial_age_of_ice_years"));
      CHKERRQ(ierr);
      tau3.set(config.get("initial_age_of_ice_years", "years", "seconds"));
  }
  
  ierr = verbPrintf(2, grid.com, "  filling ice temperatures using surface temps (and %s)\n",
                    (config.get_string("bootstrapping_temperature_heuristic") == "quartic_guess"
                     ? "quartic guess sans smb" : "mass balance for velocity estimate") );
     CHKERRQ(ierr);
  ierr = putTempAtDepth(); CHKERRQ(ierr);

  if (config.get_flag("do_cold_ice_methods") == false) {
    ierr = verbPrintf(2, grid.com,
                      "  ice enthalpy set from temperature, as cold ice (zero liquid fraction)\n");
    CHKERRQ(ierr);
  }

  return 0;
}


//! Create a temperature field within the ice from provided ice thickness, surface temperature, surface mass balance, and geothermal flux.
/*!
In bootstrapping we need to determine initial values for the temperature within
the ice (and the bedrock).  There are various data available at bootstrapping,
but not the 3D temperature field needed as initial values for the temperature.  Here 
we take a "guess" based on an assumption of steady state and a simple model of
the vertical velocity in the column.  The rule is certainly heuristic but it
seems to work well anyway.

The result is *not* the temperature field which is in steady state with the ice
dynamics.  Spinup is most-definitely needed in many applications.  Such spinup
usually starts from the temperature field computed by this procedure and then
runs for a long time (e.g. \f$10^4\f$ to \f$10^6\f$ years), possibly with fixed
geometry, to get closer to thermomechanically-coupled equilibrium.

Consider a horizontal grid point.  Suppose the surface temperature
\f$T_s\f$, surface mass balance \f$m\f$, and geothermal flux \f$g\f$ are given at that location.
Within the column denote the temperature by \f$T(z)\f$ at height \f$z\f$ above
the base of the ice.  Suppose the column of ice has height \f$H\f$, the ice
thickness.

There are two alternative bootstrap methods determined by the configuration parameter
`config.get("bootstrapping_temperature_heuristic"))`. Allowed values are `"smb"` and `"quartic_guess"`.

1. If the `smb` method is chosen, which is the default, and if \f$m>0\f$,
then the method sets the ice
temperature to the solution of the steady problem [\ref Paterson]
  \f[\rho_i c w \frac{\partial T}{\partial z} = k_i \frac{\partial^2 T}{\partial z^2} \qquad \text{with boundary conditions} \qquad T(H) = T_s \quad \text{and} \quad \frac{\partial T}{\partial z}(0) = - \frac{g}{k_i}, \f]
where the vertical velocity is linear between the surface value \f$w=-m\f$ and
a velocity of zero at the base:
  \f[w(z) = - m z / H.\f]
(Note that because \f$m>0\f$, this vertical velocity is downward.)
This is a two-point boundary value problem for a linear ODE.  In fact, if
\f$K = k_i / (\rho_i c)\f$ then we can write the ODE as
  \f[K T'' + \frac{m z}{H} T' = 0.\f]
Then let
  \f[C_0 = \frac{g \sqrt{\pi H K}}{k_i \sqrt{2 m}}, \qquad \gamma_0 = \sqrt{\frac{mH}{2K}}.\f]
(Note \f$\gamma_0\f$ is, up to a constant, the square root of the Peclet number
[\ref Paterson]; compare [\ref vanderWeletal2013].)  The solution to the
two-point boundary value problem is then
  \f[T(z) = T_s + C_0 \left(\operatorname{erf}(\gamma_0) - \operatorname{erf}\left(\gamma_0 \frac{z}{H}\right)\right).\f]
If `usesmb` is true and \f$m \le 0\f$, then the velocity in the column, relative
to the base, is taken to be zero.  Thus the solution is
  \f[ T(z) = \frac{g}{k_i} \left( H - z \right) + T_s, \f]
a straight line whose slope is determined by the geothermal flux and whose value
at the ice surface is the surface temperature, \f$T(H) = T_s\f$.
2. If the `quartic_guess` method is chosen, the "quartic guess" formula which was in older
versions of PISM is used.  Namely, within the ice we set
\f[T(z) = T_s + \alpha (H-z)^2 + \beta (H-z)^4\f]
where \f$\alpha,\beta\f$ are chosen so that
\f[\frac{\partial T}{\partial z}\Big|_{z=0} = - \frac{g}{k_i} \qquad \text{and} \qquad \frac{\partial T}{\partial z}\Big|_{z=H/4} = - \frac{g}{2 k_i}.\f]
The purpose of the second condition is that when ice is advecting downward then
the temperature gradient is much larger in roughly the bottom quarter of the
ice column.  However, without the surface mass balance, much less the solution
of the stress balance equations, we cannot estimate the vertical velocity, so
we make such a rough guess.

In either case the temperature within the ice is not allowed to exceed the
pressure-melting temperature.

We set \f$T(z)=T_s\f$ above the top of the ice.

This method determines \f$T(0)\f$, the ice temperature at the ice base.  This
temperature is used by PISMBedThermalUnit::bootstrap() to determine a
bootstrap temperature profile in the bedrock.
*/
PetscErrorCode IceModel::putTempAtDepth() {
  PetscErrorCode  ierr;

  double *T = new double[grid.Mz];
  const bool do_cold = config.get_flag("do_cold_ice_methods"),
             usesmb  = config.get_string("bootstrapping_temperature_heuristic") == "smb";
  const double
    ice_k = config.get("ice_thermal_conductivity"),
    melting_point_temp = config.get("water_melting_point_temperature"),
    ice_density = config.get("ice_density"),
    beta_CC_grad = config.get("beta_CC") * ice_density * config.get("standard_gravity"),
    KK = ice_k / (ice_density * config.get("ice_specific_heat_capacity"));

  assert(surface != NULL);

  // Note that surface->update was called in bootstrapFromFile()
  {
    ierr = surface->ice_surface_temperature(ice_surface_temp); CHKERRQ(ierr);
    if (usesmb == true) {
      ierr = surface->ice_surface_mass_flux(climatic_mass_balance); CHKERRQ(ierr);
      // convert from [kg m-2 s-1] to [m / s]
      ierr = climatic_mass_balance.scale(1.0 / config.get("ice_density")); CHKERRQ(ierr);
    }
  }

  IceModelVec3 *result;
  if (do_cold)
    result = &T3;
  else
    result = &Enth3;

  ierr = ice_surface_temp.begin_access();      CHKERRQ(ierr);
  ierr = climatic_mass_balance.begin_access(); CHKERRQ(ierr);
  ierr = ice_thickness.begin_access();         CHKERRQ(ierr);
  ierr = geothermal_flux.begin_access();       CHKERRQ(ierr);
  ierr = result->begin_access();               CHKERRQ(ierr);

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const double HH = ice_thickness(i,j),
                        Ts = ice_surface_temp(i,j),
                        gg = geothermal_flux(i,j);
      const unsigned int ks = grid.kBelowHeight(HH);

      // within ice
      if (usesmb == true) { // method 1:  includes surface mass balance in estimate
        const double mm = climatic_mass_balance(i,j);
        if (mm <= 0.0) { // negative or zero surface mass balance case: linear
          for (unsigned int k = 0; k < ks; k++) {
            const double z = grid.zlevels[k],
              Tpmp = melting_point_temp - beta_CC_grad * (HH - z);
            T[k] = gg / ice_k * (HH - z) + Ts;
            T[k] = PetscMin(Tpmp,T[k]);
          }
        } else { // positive surface mass balance case
          const double C0 = (gg * sqrt(M_PI * HH * KK)) / (ice_k * sqrt(2.0 * mm)),
            gamma0 = sqrt(mm * HH / (2.0 * KK));

          for (unsigned int k = 0; k < ks; k++) {
            const double z = grid.zlevels[k],
              Tpmp = melting_point_temp - beta_CC_grad * (HH - z);
            T[k] = Ts + C0 * ( erf(gamma0) - erf(gamma0 * z / HH) );
            T[k] = PetscMin(Tpmp,T[k]);
          }
        }
      } else { // method 2:  does not use smb
        const double beta = (4.0/21.0) * (gg / (2.0 * ice_k * HH * HH * HH)),
                          alpha = (gg / (2.0 * HH * ice_k)) - 2.0 * HH * HH * beta;
        for (unsigned int k = 0; k < ks; k++) {
          const double depth = HH - grid.zlevels[k],
                            Tpmp = melting_point_temp - beta_CC_grad * depth,
                            d2 = depth * depth;
          T[k] = PetscMin(Tpmp, Ts + alpha * d2 + beta * d2 * d2);
        }
      }

      // above ice
      for (unsigned int k = ks; k < grid.Mz; k++)
        T[k] = ice_surface_temp(i,j);

      // convert to enthalpy if that's what we are calculating
      if (!do_cold) {
        for (unsigned int k = 0; k < grid.Mz; ++k) {
          const double depth = HH - grid.zlevels[k];
          const double pressure = EC->getPressureFromDepth(depth);
          // reuse T to store enthalpy; assume that the ice is cold
          ierr = EC->getEnthPermissive(T[k], 0.0, pressure, T[k]); CHKERRQ(ierr);
        }
      }

      ierr = result->setInternalColumn(i,j,T); CHKERRQ(ierr);

    }
  }

  ierr = ice_thickness.end_access();         CHKERRQ(ierr);
  ierr = geothermal_flux.end_access();       CHKERRQ(ierr);
  ierr = result->end_access();               CHKERRQ(ierr);
  ierr = ice_surface_temp.end_access();      CHKERRQ(ierr);
  ierr = climatic_mass_balance.end_access(); CHKERRQ(ierr);

  delete [] T;

  ierr = result->update_ghosts(); CHKERRQ(ierr);

  if (do_cold) {
    ierr = compute_enthalpy_cold(T3, Enth3); CHKERRQ(ierr);
  }

  return 0;
}

