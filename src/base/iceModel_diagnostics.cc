// Copyright (C) 2010, 2011, 2012, 2013, 2014 Constantine Khroulev
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

#include "pism_options.hh"
#include "iceModel_diagnostics.hh"
#include "PISMDiagnostic.hh"
#include "Mask.hh"

#include "PISMBedDef.hh"
#include "PISMYieldStress.hh"
#include "PISMHydrology.hh"
#include "PISMStressBalance.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "enthalpyConverter.hh"
#include "ShallowStressBalance.hh"
#include "SSB_Modifier.hh"
#include "bedrockThermalUnit.hh"
#include "IceModelVec3Custom.hh"

PetscErrorCode IceModel::init_diagnostics() {
  bool print_list_and_stop = false;

  PetscErrorCode ierr = PISMOptionsIsSet("-list_diagnostics",
                                         "List available diagnostic quantities and stop",
                                         print_list_and_stop); CHKERRQ(ierr);

  // Add IceModel diagnostics:
  diagnostics["cts"]              = new IceModel_cts(this, grid, variables);
  diagnostics["enthalpybase"]     = new IceModel_enthalpybase(this, grid, variables);
  diagnostics["enthalpysurf"]     = new IceModel_enthalpysurf(this, grid, variables);
  diagnostics["hardav"]           = new IceModel_hardav(this, grid, variables);
  diagnostics["liqfrac"]          = new IceModel_liqfrac(this, grid, variables);
  diagnostics["proc_ice_area"]    = new IceModel_proc_ice_area(this, grid, variables);
  diagnostics["rank"]             = new IceModel_rank(this, grid, variables);
  diagnostics["temp"]             = new IceModel_temp(this, grid, variables);
  diagnostics["temp_pa"]          = new IceModel_temp_pa(this, grid, variables);
  diagnostics["tempbase"]         = new IceModel_tempbase(this, grid, variables);
  diagnostics["tempicethk"]       = new IceModel_tempicethk(this, grid, variables);
  diagnostics["tempicethk_basal"] = new IceModel_tempicethk_basal(this, grid, variables);
  diagnostics["temppabase"]       = new IceModel_temppabase(this, grid, variables);
  diagnostics["tempsurf"]         = new IceModel_tempsurf(this, grid, variables);
  diagnostics["dHdt"]             = new IceModel_dHdt(this, grid, variables);

  if (flux_divergence.was_created() ||
      print_list_and_stop) {
    diagnostics["flux_divergence"] = new IceModel_flux_divergence(this, grid, variables);
  }

  if (climatic_mass_balance_cumulative.was_created() ||
      print_list_and_stop) {
    diagnostics["climatic_mass_balance_cumulative"] = new IceModel_climatic_mass_balance_cumulative(this, grid, variables);
  }

  if (nonneg_flux_2D_cumulative.was_created() ||
      print_list_and_stop) {
    diagnostics["nonneg_flux_cumulative"] = new IceModel_nonneg_flux_2D_cumulative(this, grid, variables);
  }

  if (grounded_basal_flux_2D_cumulative.was_created() ||
      print_list_and_stop) {
    diagnostics["grounded_basal_flux_cumulative"] = new IceModel_grounded_basal_flux_2D_cumulative(this, grid, variables);
  }

  if (floating_basal_flux_2D_cumulative.was_created() ||
      print_list_and_stop) {
    diagnostics["floating_basal_flux_cumulative"] = new IceModel_floating_basal_flux_2D_cumulative(this, grid, variables);
  }

  if (discharge_flux_2D_cumulative.was_created() ||
      print_list_and_stop) {
    diagnostics["discharge_flux_cumulative"] = new IceModel_discharge_flux_2D_cumulative(this, grid, variables);
  }

#if (PISM_USE_PROJ4==1)
  if (mapping.has_attribute("proj4")) {
    std::string proj4 = mapping.get_string("proj4");
    diagnostics["lat_bounds"] = new IceModel_lat_lon_bounds(this, grid, variables, "lat", proj4);
    diagnostics["lon_bounds"] = new IceModel_lat_lon_bounds(this, grid, variables, "lon", proj4);
  }
#elif (PISM_USE_PROJ4==0)
  // do nothing
#else  // PISM_USE_PROJ4 is not set
#error "PISM build system error: PISM_USE_PROJ4 is not set."
#endif

  ts_diagnostics["ivol"]          = new IceModel_ivol(this, grid, variables);
  ts_diagnostics["slvol"]         = new IceModel_slvol(this, grid, variables);
  ts_diagnostics["divoldt"]       = new IceModel_divoldt(this, grid, variables);
  ts_diagnostics["iarea"]         = new IceModel_iarea(this, grid, variables);
  ts_diagnostics["imass"]         = new IceModel_imass(this, grid, variables);
  ts_diagnostics["dimassdt"]      = new IceModel_dimassdt(this, grid, variables);
  ts_diagnostics["ivoltemp"]      = new IceModel_ivoltemp(this, grid, variables);
  ts_diagnostics["ivolcold"]      = new IceModel_ivolcold(this, grid, variables);
  ts_diagnostics["ivolg"]         = new IceModel_ivolg(this, grid, variables);
  ts_diagnostics["ivolf"]         = new IceModel_ivolf(this, grid, variables);
  ts_diagnostics["iareatemp"]     = new IceModel_iareatemp(this, grid, variables);
  ts_diagnostics["iareacold"]     = new IceModel_iareacold(this, grid, variables);
  ts_diagnostics["iareag"]        = new IceModel_iareag(this, grid, variables);
  ts_diagnostics["iareaf"]        = new IceModel_iareaf(this, grid, variables);
  ts_diagnostics["dt"]            = new IceModel_dt(this, grid, variables);
  ts_diagnostics["max_diffusivity"] = new IceModel_max_diffusivity(this, grid, variables);
  ts_diagnostics["ienthalpy"]     = new IceModel_ienthalpy(this, grid, variables);
  ts_diagnostics["max_hor_vel"]   = new IceModel_max_hor_vel(this, grid, variables);

  ts_diagnostics["surface_ice_flux"]   = new IceModel_surface_flux(this, grid, variables);
  ts_diagnostics["surface_ice_flux_cumulative"]   = new IceModel_surface_flux_cumulative(this, grid, variables);
  ts_diagnostics["grounded_basal_ice_flux"]     = new IceModel_grounded_basal_flux(this, grid, variables);
  ts_diagnostics["grounded_basal_ice_flux_cumulative"]     = new IceModel_grounded_basal_flux_cumulative(this, grid, variables);
  ts_diagnostics["sub_shelf_ice_flux"] = new IceModel_sub_shelf_flux(this, grid, variables);
  ts_diagnostics["sub_shelf_ice_flux_cumulative"] = new IceModel_sub_shelf_flux_cumulative(this, grid, variables);
  ts_diagnostics["nonneg_rule_flux"]   = new IceModel_nonneg_flux(this, grid, variables);
  ts_diagnostics["nonneg_rule_flux_cumulative"]   = new IceModel_nonneg_flux_cumulative(this, grid, variables);
  ts_diagnostics["discharge_flux"]    = new IceModel_discharge_flux(this, grid, variables);
  ts_diagnostics["discharge_flux_cumulative"]    = new IceModel_discharge_flux_cumulative(this, grid, variables);
  ts_diagnostics["H_to_Href_flux"] = new IceModel_H_to_Href_flux(this, grid, variables);
  ts_diagnostics["Href_to_H_flux"] = new IceModel_Href_to_H_flux(this, grid, variables);
  ts_diagnostics["sum_divQ_flux"]  = new IceModel_sum_divQ_flux(this, grid, variables);

  // Get diagnostics supported by the stress balance object:
  if (stress_balance != NULL)
    stress_balance->get_diagnostics(diagnostics, ts_diagnostics);

  // Get diagnostics supported by the surface model:
  if (surface != NULL)
    surface->get_diagnostics(diagnostics, ts_diagnostics);

  // Get diagnostics supported by the ocean model:
  if (ocean != NULL)
    ocean->get_diagnostics(diagnostics, ts_diagnostics);

  // Get diagnostics supported by the bed deformation model:
  if (beddef != NULL)
    beddef->get_diagnostics(diagnostics, ts_diagnostics);

  if (basal_yield_stress_model != NULL)
    basal_yield_stress_model->get_diagnostics(diagnostics, ts_diagnostics);

  if (subglacial_hydrology != NULL)
    subglacial_hydrology->get_diagnostics(diagnostics, ts_diagnostics);

  if (print_list_and_stop) {
    ierr = list_diagnostics(); CHKERRQ(ierr);

    PISMEnd();
  }

  return 0;
}

PetscErrorCode IceModel::list_diagnostics() {

  PetscPrintf(grid.com, "\n");

  // quantities with dedicated storage
  {
    std::set<std::string> list = variables.keys();

    if (beddef != NULL)
      beddef->add_vars_to_output("big", list);

    if (btu != NULL)
      btu->add_vars_to_output("big", list);

    if (basal_yield_stress_model != NULL)
      basal_yield_stress_model->add_vars_to_output("big", list);

    if (subglacial_hydrology != NULL)
      subglacial_hydrology->add_vars_to_output("big", list);

    if (stress_balance != NULL)
      stress_balance->add_vars_to_output("big", list);

    if (ocean != NULL)
      ocean->add_vars_to_output("big", list);

    if (surface != NULL)
      surface->add_vars_to_output("big", list);

    for (unsigned int d = 3; d > 1; --d) {

      PetscPrintf(grid.com,
                  "======== Available %dD quantities with dedicated storage ========\n",
                  d);

      std::set<std::string>::iterator j = list.begin();
      while(j != list.end()) {
        IceModelVec *v = variables.get(*j);

        if (v != NULL && v->get_ndims() == d) {
          NCSpatialVariable var = v->metadata();

          std::string
            name                = var.get_name(),
            units               = var.get_string("units"),
            glaciological_units = var.get_string("glaciological_units"),
            long_name           = var.get_string("long_name");

          if (glaciological_units.empty() == false)
            units = glaciological_units;

            PetscPrintf(grid.com,
                        "   Name: %s [%s]\n"
                        "       - %s\n\n", name.c_str(), units.c_str(), long_name.c_str());
        }

        ++j;
      }
    }

  } // done with quantities with dedicated storage

  // 2D and 3D diagnostics
  for (unsigned int d = 3; d > 1; --d) {

    PetscPrintf(grid.com,
                "======== Available %dD diagnostic quantities ========\n",
                d);

    std::map<std::string, PISMDiagnostic*>::iterator j = diagnostics.begin();
    while (j != diagnostics.end()) {
      PISMDiagnostic *diag = j->second;

      std::string name           = j->first,
        units               = diag->get_metadata().get_string("units"),
        glaciological_units = diag->get_metadata().get_string("glaciological_units");

      if (glaciological_units.empty() == false)
        units = glaciological_units;

      if (diag->get_metadata().get_n_spatial_dimensions() == d) {

        PetscPrintf(grid.com, "   Name: %s [%s]\n", name.c_str(), units.c_str());

        for (int k = 0; k < diag->get_nvars(); ++k) {
          NCSpatialVariable var = diag->get_metadata(k);

          std::string long_name = var.get_string("long_name");

          PetscPrintf(grid.com, "      -  %s\n", long_name.c_str());
        }

        PetscPrintf(grid.com, "\n");

      }

      ++j;
    }
  }

  // scalar time-series
  PetscPrintf(grid.com, "======== Available time-series ========\n");

  std::map<std::string, PISMTSDiagnostic*>::iterator j = ts_diagnostics.begin();
  while (j != ts_diagnostics.end()) {
    PISMTSDiagnostic *diag = j->second;

    std::string name = j->first,
      long_name = diag->get_string("long_name"),
      units = diag->get_string("units"),
      glaciological_units = diag->get_string("glaciological_units");

    if (glaciological_units.empty() == false)
      units = glaciological_units;

    PetscPrintf(grid.com,
                "   Name: %s [%s]\n"
                "      -  %s\n\n",
                name.c_str(), units.c_str(), long_name.c_str());

    ++j;
  }

  return 0;
}


IceModel_hardav::IceModel_hardav(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("hardav", grid);

  const double power = 1.0 / grid.config.get("Glen_exponent");
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);

  set_attrs("vertical average of ice hardness", "",
            unitstr, unitstr, 0);

  vars[0].set_double("valid_min", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

//! \brief Computes vertically-averaged ice hardness.
PetscErrorCode IceModel_hardav::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  const double fillval = grid.config.get("fill_value");
  double *Eij; // columns of enthalpy values

  IceFlowLaw *flow_law = model->stress_balance->get_stressbalance()->get_flow_law();
  if (flow_law == NULL) {
    flow_law = model->stress_balance->get_ssb_modifier()->get_flow_law();
    if (flow_law == NULL) {
      PetscPrintf(grid.com, "ERROR: Can't compute vertically-averaged hardness: no flow law is used.\n");
      PISMEnd();
    }
  }

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "hardav", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  MaskQuery mask(model->vMask);

  ierr = model->vMask.begin_access(); CHKERRQ(ierr);
  ierr = model->Enth3.begin_access(); CHKERRQ(ierr);
  ierr = model->ice_thickness.begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = model->Enth3.getInternalColumn(i,j,&Eij); CHKERRQ(ierr);
      const double H = model->ice_thickness(i,j);
      if (mask.icy(i, j)) {
        (*result)(i,j) = flow_law->averaged_hardness(H, grid.kBelowHeight(H),
                                                     &grid.zlevels[0], Eij);
      } else { // put negative value below valid range
        (*result)(i,j) = fillval;
      }
    }
  }
  ierr = model->Enth3.end_access(); CHKERRQ(ierr);
  ierr = model->ice_thickness.end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = model->vMask.end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}


IceModel_rank::IceModel_rank(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("rank", grid);

  set_attrs("processor rank", "", "1", "", 0);
  vars[0].set_time_independent(true);
}

PetscErrorCode IceModel_rank::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "rank", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  ierr = result->begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i)
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j)
      (*result)(i,j) = grid.rank;
  ierr = result->end_access();

  output = result;
  return 0;
}


IceModel_cts::IceModel_cts(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_3d("cts", grid, g.zlevels);

  set_attrs("cts = E/E_s(p), so cold-temperate transition surface is at cts = 1", "",
            "", "", 0);
}

PetscErrorCode IceModel_cts::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  // update vertical levels (in case the grid was extended
  vars[0].set_levels(grid.zlevels);

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "cts", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  ierr = model->setCTSFromEnthalpy(*result); CHKERRQ(ierr);

  output = result;
  return 0;
}

IceModel_proc_ice_area::IceModel_proc_ice_area(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("proc_ice_area", grid);

  set_attrs("number of cells containing ice in a processor's domain", "",
            "", "", 0);
  vars[0].set_time_independent(true);
}

PetscErrorCode IceModel_proc_ice_area::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *thickness;

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  IceModelVec2Int *ice_mask;
  ice_mask = dynamic_cast<IceModelVec2Int*>(variables.get("mask"));
  if (ice_mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "proc_ice_area", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  int ice_filled_cells = 0;

  MaskQuery mask(*ice_mask);

  ierr = ice_mask->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i)
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j)
      if (mask.icy(i, j)) {
        ice_filled_cells += 1;
      }
  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = ice_mask->end_access(); CHKERRQ(ierr);

  ierr = result->begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      (*result)(i,j) = ice_filled_cells;
    }
  }
  ierr = result->end_access();

  output = result;
  return 0;
}


IceModel_temp::IceModel_temp(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_3d("temp", grid, g.zlevels);

  set_attrs("ice temperature", "land_ice_temperature", "K", "K", 0);
  vars[0].set_double("valid_min", 0);
}

PetscErrorCode IceModel_temp::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  // update vertical levels (in case the grid was extended
  vars[0].set_levels(grid.zlevels);

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "temp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  IceModelVec2S *thickness;
  IceModelVec3 *enthalpy;

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(variables.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(grid.com, 1, "enthalpy is not available");

  double *Tij, *Enthij; // columns of these values
  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = result->getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = enthalpy->getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (unsigned int k=0; k <grid.Mz; ++k) {
        const double depth = (*thickness)(i,j) - grid.zlevels[k];
        ierr = model->EC->getAbsTemp(Enthij[k],
                                     model->EC->getPressureFromDepth(depth),
                                     Tij[k]);
        if (ierr) {
          PetscPrintf(grid.com,
                      "\n\nEnthalpyConverter.getAbsTemp() error at i=%d,j=%d,k=%d\n\n",
                      i,j,k);
        }
        CHKERRQ(ierr);
      }
    }
  }
  ierr = enthalpy->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}


IceModel_temp_pa::IceModel_temp_pa(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_3d("temp_pa", grid, g.zlevels);

  set_attrs("pressure-adjusted ice temperature (degrees above pressure-melting point)", "",
            "deg_C", "deg_C", 0);
  vars[0].set_double("valid_max", 0);
}

PetscErrorCode IceModel_temp_pa::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  bool cold_mode = grid.config.get_flag("do_cold_ice_methods");
  double melting_point_temp = grid.config.get("water_melting_point_temperature");

  // update vertical levels (in case the grid was extended
  vars[0].set_levels(grid.zlevels);

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "temp_pa", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  IceModelVec2S *thickness;
  IceModelVec3 *enthalpy;

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(variables.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(grid.com, 1, "enthalpy is not available");

  double *Tij, *Enthij; // columns of these values
  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = result->getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = enthalpy->getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (unsigned int k=0; k < grid.Mz; ++k) {
        const double depth = (*thickness)(i,j) - grid.zlevels[k],
          p = model->EC->getPressureFromDepth(depth);
        ierr = model->EC->getPATemp(Enthij[k], p, Tij[k]);
        if (ierr) {
          PetscPrintf(grid.com,
                      "\n\nEnthalpyConverter.getAbsTemp() error at i=%d,j=%d,k=%d\n\n",
                      i,j,k);
        }
        CHKERRQ(ierr);

        if (cold_mode) { // if ice is temperate then its pressure-adjusted temp
          // is 273.15
          if ( model->EC->isTemperate(Enthij[k],p) && ((*thickness)(i,j) > 0)) {
            Tij[k] = melting_point_temp;
          }
        }

      }
    }
  }
  ierr = enthalpy->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);

  ierr = result->shift(-melting_point_temp); CHKERRQ(ierr);

  output = result;
  return 0;
}

IceModel_temppabase::IceModel_temppabase(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("temppabase", grid);

  set_attrs("pressure-adjusted ice temperature at the base of ice", "",
            "Celsius", "Celsius", 0);
}

PetscErrorCode IceModel_temppabase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  bool cold_mode = grid.config.get_flag("do_cold_ice_methods");
  double melting_point_temp = grid.config.get("water_melting_point_temperature");

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "temp_pa_base", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  IceModelVec2S *thickness;
  IceModelVec3 *enthalpy;

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(variables.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(grid.com, 1, "enthalpy is not available");

  double *Enthij; // columns of these values
  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = enthalpy->getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);

      const double depth = (*thickness)(i,j),
        p = model->EC->getPressureFromDepth(depth);
      ierr = model->EC->getPATemp(Enthij[0], p,
                                  (*result)(i,j));
      if (ierr) {
        PetscPrintf(grid.com,
                    "\n\nEnthalpyConverter.getAbsTemp() error at i=%d,j=%d\n\n",
                    i,j);
      }
      CHKERRQ(ierr);

      if (cold_mode) { // if ice is temperate then its pressure-adjusted temp
        // is 273.15
        if ( model->EC->isTemperate(Enthij[0],p) && ((*thickness)(i,j) > 0)) {
          (*result)(i,j) = melting_point_temp;
        }
      }
    }
  }
  ierr = enthalpy->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);

  ierr = result->shift(-melting_point_temp); CHKERRQ(ierr);

  output = result;
  return 0;
}

IceModel_enthalpysurf::IceModel_enthalpysurf(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("enthalpysurf", grid);

  set_attrs("ice enthalpy at 1m below the ice surface", "",
            "J kg-1", "J kg-1", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

PetscErrorCode IceModel_enthalpysurf::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "enthalpysurf", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  double fill_value = grid.config.get("fill_value");

  // compute levels corresponding to 1 m below the ice surface:

  ierr = model->ice_thickness.begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      (*result)(i,j) = PetscMax(model->ice_thickness(i,j) - 1.0, 0.0);
    }
  }
  ierr = result->end_access(); CHKERRQ(ierr);

  ierr = model->Enth3.getSurfaceValues(*result, *result); CHKERRQ(ierr);  // z=0 slice

  ierr = result->begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (model->ice_thickness(i,j) <= 1.0)
        (*result)(i,j) = fill_value;
    }
  }
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = model->ice_thickness.end_access(); CHKERRQ(ierr);


  output = result;
  return 0;
}

IceModel_enthalpybase::IceModel_enthalpybase(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("enthalpybase", grid);

  set_attrs("ice enthalpy at the base of ice", "",
            "J kg-1", "J kg-1", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

PetscErrorCode IceModel_enthalpybase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "enthalpybase", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  ierr = model->Enth3.getHorSlice(*result, 0.0); CHKERRQ(ierr);  // z=0 slice

  ierr = result->mask_by(model->ice_thickness, grid.config.get("fill_value")); CHKERRQ(ierr);

  output = result;
  return 0;
}


IceModel_tempbase::IceModel_tempbase(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("tempbase", grid);

  set_attrs("ice temperature at the base of ice", "",
            "K", "K", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

PetscErrorCode IceModel_tempbase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result, *thickness;

  IceModel_enthalpybase enth(model, grid, variables);

  ierr = enth.compute(output); CHKERRQ(ierr);
  result = dynamic_cast<IceModelVec2S*>(output);
  if (result == NULL) SETERRQ(grid.com, 1, "dynamic_cast failure");

  // result contains basal enthalpy; note that it is allocated by
  // enth.compute().

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  MaskQuery mask(model->vMask);

  ierr = model->vMask.begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);

  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      double depth = (*thickness)(i,j),
        pressure = model->EC->getPressureFromDepth(depth);
      if (mask.icy(i, j)) {
        ierr = model->EC->getAbsTemp((*result)(i,j),
                                     pressure,
                                     (*result)(i,j));
        CHKERRQ(ierr);
      } else {
        (*result)(i,j) = grid.config.get("fill_value");
      }
    }
  }

  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = model->vMask.end_access(); CHKERRQ(ierr);

  result->metadata() = vars[0];
  output = result;
  return 0;
}

IceModel_tempsurf::IceModel_tempsurf(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("tempsurf", grid);

  set_attrs("ice temperature at 1m below the ice surface", "",
            "K", "K", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

PetscErrorCode IceModel_tempsurf::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result, *thickness;

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  IceModel_enthalpysurf enth(model, grid, variables);

  ierr = enth.compute(output); CHKERRQ(ierr);
  result = dynamic_cast<IceModelVec2S*>(output);
  if (result == NULL) SETERRQ(grid.com, 1, "dynamic_cast failure");

  // result contains surface enthalpy; note that it is allocated by
  // enth.compute().

  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);

  double depth = 1.0,
    pressure = model->EC->getPressureFromDepth(depth);
  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if ((*thickness)(i,j) > 1) {
        ierr = model->EC->getAbsTemp((*result)(i,j),
                                     pressure,
                                     (*result)(i,j)); CHKERRQ(ierr);
      } else {
        (*result)(i,j) = grid.config.get("fill_value");
      }
    }
  }

  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);

  result->metadata() = vars[0];
  output = result;
  return 0;
}


IceModel_liqfrac::IceModel_liqfrac(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_3d("liqfrac", grid, g.zlevels);

  set_attrs("liquid water fraction in ice (between 0 and 1)", "",
            "1", "1", 0);
  vars[0].set_double("valid_min", 0);
  vars[0].set_double("valid_max", 1);
}

PetscErrorCode IceModel_liqfrac::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  // update vertical levels (in case the grid was extended
  vars[0].set_levels(grid.zlevels);

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "liqfrac", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  bool cold_mode = grid.config.get_flag("do_cold_ice_methods");

  if (cold_mode) {
    ierr = result->set(0.0); CHKERRQ(ierr);
  } else {
    ierr = model->compute_liquid_water_fraction(model->Enth3, *result); CHKERRQ(ierr);
  }

  output = result;
  return 0;
}

IceModel_tempicethk::IceModel_tempicethk(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("tempicethk", grid);

  set_attrs("temperate ice thickness (total column content)", "",
            "m", "m", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

PetscErrorCode IceModel_tempicethk::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "tempicethk", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  double *Enth;

  MaskQuery mask(model->vMask);

  ierr = model->vMask.begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = model->Enth3.begin_access(); CHKERRQ(ierr);
  ierr = model->ice_thickness.begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (mask.icy(i, j)) {
        ierr = model->Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
        double temperate_ice_thickness = 0.0;
        double ice_thickness = model->ice_thickness(i,j);
        const unsigned int ks = grid.kBelowHeight(ice_thickness);

        for (unsigned int k=0; k<ks; ++k) { // FIXME issue #15
          double pressure = model->EC->getPressureFromDepth(ice_thickness - grid.zlevels[k]);

          if (model->EC->isTemperate(Enth[k], pressure)) {
            temperate_ice_thickness += grid.zlevels[k+1] - grid.zlevels[k];
          }
        }

        double pressure = model->EC->getPressureFromDepth(ice_thickness - grid.zlevels[ks]);
        if (model->EC->isTemperate(Enth[ks], pressure)) {
          temperate_ice_thickness += ice_thickness - grid.zlevels[ks];
        }

        (*result)(i,j) = temperate_ice_thickness;
      } else {
        // ice-free
        (*result)(i,j) = grid.config.get("fill_value");
      }
    }
  }
  ierr = model->Enth3.end_access(); CHKERRQ(ierr);
  ierr = model->ice_thickness.end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = model->vMask.end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}

IceModel_tempicethk_basal::IceModel_tempicethk_basal(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("tempicethk_basal", grid);

  set_attrs("thickness of the basal layer of temperate ice", "",
            "m", "m", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

/*!
 * Uses linear interpolation to go beyond vertical grid resolution.
 */
PetscErrorCode IceModel_tempicethk_basal::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "tempicethk_basal", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  double *Enth;
  EnthalpyConverter *EC = model->EC;

  MaskQuery mask(model->vMask);

  const double fill_value = grid.config.get("fill_value");

  ierr = model->vMask.begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = model->ice_thickness.begin_access(); CHKERRQ(ierr);
  ierr = model->Enth3.begin_access(); CHKERRQ(ierr);
  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      double thk = model->ice_thickness(i,j);

      // if we have no ice, go on to the next grid point (this cell will be
      // marked as "missing" later)
      if (mask.ice_free(i, j)) {
        (*result)(i,j) = fill_value;
        continue;
      }

      ierr = model->Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
      double pressure;
      unsigned int ks = grid.kBelowHeight(thk),
        k = 0;

      while (k <= ks) {         // FIXME issue #15
        pressure = EC->getPressureFromDepth(thk - grid.zlevels[k]);

        if (EC->isTemperate(Enth[k],pressure))
          k++;
        else
          break;
      }
      // after this loop 'pressure' is equal to the pressure at the first level
      // that is cold

      // no temperate ice at all; go to the next grid point
      if (k == 0) {
        (*result)(i,j) = 0.0;
        continue;
      }

      // the whole column is temperate (except, possibly, some ice between
      // zlevels[ks] and the total thickness; we ignore it)
      if (k == ks + 1) {
        (*result)(i,j) = grid.zlevels[ks];
        continue;
      }

      double
        pressure_0 = EC->getPressureFromDepth(thk - grid.zlevels[k-1]),
        dz         = grid.zlevels[k] - grid.zlevels[k-1],
        slope1     = (Enth[k] - Enth[k-1]) / dz,
        slope2     = (EC->getEnthalpyCTS(pressure) - EC->getEnthalpyCTS(pressure_0)) / dz;

      if (slope1 != slope2) {
        (*result)(i,j) = grid.zlevels[k-1] +
          (EC->getEnthalpyCTS(pressure_0) - Enth[k-1]) / (slope1 - slope2);

        // check if the resulting thickness is valid:
        (*result)(i,j) = PetscMax((*result)(i,j), grid.zlevels[k-1]);
        (*result)(i,j) = PetscMin((*result)(i,j), grid.zlevels[k]);
      } else {
        SETERRQ4(grid.com, 1, "This should never happen: (i=%d, j=%d, k=%d, ks=%d)\n",
                 i, j, k, ks);
      }
    }
  }
  ierr = model->Enth3.end_access(); CHKERRQ(ierr);
  ierr = model->ice_thickness.end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = model->vMask.end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}

IceModel_flux_divergence::IceModel_flux_divergence(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("flux_divergence", grid);

  set_attrs("flux divergence", "", "m s-1", "m year-1", 0);
}

PetscErrorCode IceModel_flux_divergence::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "flux_divergence", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  ierr = result->copy_from(model->flux_divergence); CHKERRQ(ierr);

  output = result;
  return 0;
}

IceModel_climatic_mass_balance_cumulative::IceModel_climatic_mass_balance_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("climatic_mass_balance_cumulative", grid);

  set_attrs("cumulative ice-equivalent climatic mass balance", "",
            "kg m-2", "kg m-2", 0);
}

PetscErrorCode IceModel_climatic_mass_balance_cumulative::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "climatic_mass_balance_cumulative", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  ierr = result->copy_from(model->climatic_mass_balance_cumulative); CHKERRQ(ierr);

  output = result;
  return 0;
}

IceModel_ivol::IceModel_ivol(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "ivol", time_dimension_name);

  ts->set_units("m3", "");
  ts->set_dimension_units(time_units, "");

  ts->set_attr("long_name", "total ice volume");
  ts->set_attr("valid_min", 0.0);
}

PetscErrorCode IceModel_ivol::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->compute_ice_volume(value); CHKERRQ(ierr);

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_slvol::IceModel_slvol(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "slvol", time_dimension_name);

  ts->set_units("m", "");
  ts->set_dimension_units(time_units, "");

  ts->set_attr("long_name", "total sea-level relevant ice IN SEA-LEVEL EQUIVALENT");
  ts->set_attr("valid_min", 0.0);
}

PetscErrorCode IceModel_slvol::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->compute_sealevel_volume(value); CHKERRQ(ierr);

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_divoldt::IceModel_divoldt(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "divoldt", time_dimension_name);

  ts->set_units("m3 s-1", "");
  ts->set_dimension_units(time_units, "");
  ts->rate_of_change = true;

  ts->set_attr("long_name", "total ice volume rate of change");
}

PetscErrorCode IceModel_divoldt::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->compute_ice_volume(value); CHKERRQ(ierr);

  // note that "value" below *should* be the ice volume
  ts->append(value, a, b);

  return 0;
}


IceModel_iarea::IceModel_iarea(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "iarea", time_dimension_name);

  ts->set_units("m2", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "total ice area");
  ts->set_attr("valid_min", 0.0);
}

PetscErrorCode IceModel_iarea::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->compute_ice_area(value); CHKERRQ(ierr);

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_imass::IceModel_imass(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "imass", time_dimension_name);

  ts->set_units("kg", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "total ice mass");
  ts->set_attr("valid_min", 0.0);
}

PetscErrorCode IceModel_imass::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->compute_ice_volume(value); CHKERRQ(ierr);

  ierr = ts->append(value * grid.config.get("ice_density"), a, b); CHKERRQ(ierr);

  return 0;
}


IceModel_dimassdt::IceModel_dimassdt(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "dimassdt", time_dimension_name);

  ts->set_units("kg s-1", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "total ice mass rate of change");

  ts->rate_of_change = true;
}

PetscErrorCode IceModel_dimassdt::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->compute_ice_volume(value); CHKERRQ(ierr);

  ierr = ts->append(value * grid.config.get("ice_density"), a, b); CHKERRQ(ierr);

  return 0;
}


IceModel_ivoltemp::IceModel_ivoltemp(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "ivoltemp", time_dimension_name);

  ts->set_units("m3", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "total volume of temperate ice");
  ts->set_attr("valid_min", 0.0);
}

PetscErrorCode IceModel_ivoltemp::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->compute_ice_volume_temperate(value); CHKERRQ(ierr);

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}


IceModel_ivolcold::IceModel_ivolcold(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "ivolcold", time_dimension_name);

  ts->set_units("m3", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "total volume of cold ice");
  ts->set_attr("valid_min", 0.0);
}

PetscErrorCode IceModel_ivolcold::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->compute_ice_volume_cold(value); CHKERRQ(ierr);

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_iareatemp::IceModel_iareatemp(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "iareatemp", time_dimension_name);

  ts->set_units("m2", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "ice-covered area where basal ice is temperate");
  ts->set_attr("valid_min", 0.0);
}

PetscErrorCode IceModel_iareatemp::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->compute_ice_area_temperate(value); CHKERRQ(ierr);

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_iareacold::IceModel_iareacold(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "iareacold", time_dimension_name);

  ts->set_units("m2", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "ice-covered area where basal ice is cold");
  ts->set_attr("valid_min", 0.0);
}

PetscErrorCode IceModel_iareacold::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->compute_ice_area_cold(value); CHKERRQ(ierr);

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_ienthalpy::IceModel_ienthalpy(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "ienthalpy", time_dimension_name);

  ts->set_units("J", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "total ice enthalpy");
  ts->set_attr("valid_min", 0.0);
}

PetscErrorCode IceModel_ienthalpy::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->compute_ice_enthalpy(value); CHKERRQ(ierr);

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_iareag::IceModel_iareag(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "iareag", time_dimension_name);

  ts->set_units("m2", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "total grounded ice area");
}

PetscErrorCode IceModel_iareag::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->compute_ice_area_grounded(value); CHKERRQ(ierr);

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_iareaf::IceModel_iareaf(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "iareaf", time_dimension_name);

  ts->set_units("m2", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "total floating ice area");
}

PetscErrorCode IceModel_iareaf::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->compute_ice_area_floating(value); CHKERRQ(ierr);

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_dt::IceModel_dt(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "dt", time_dimension_name);

  ts->set_units("second", "year");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "mass continuity time step");
}

PetscErrorCode IceModel_dt::update(double a, double b) {
  PetscErrorCode ierr;

  ierr = ts->append(model->dt, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_max_diffusivity::IceModel_max_diffusivity(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "max_diffusivity", time_dimension_name);

  ts->set_units("m2 s-1", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "maximum diffusivity");
}

PetscErrorCode IceModel_max_diffusivity::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  ierr = model->stress_balance->get_max_diffusivity(value); CHKERRQ(ierr);

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_surface_flux::IceModel_surface_flux(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "surface_ice_flux", time_dimension_name);

  ts->set_units("kg s-1", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "total over ice domain of top surface ice mass flux");
  ts->rate_of_change = true;
}

PetscErrorCode IceModel_surface_flux::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  value = model->surface_ice_flux_cumulative;

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_surface_flux_cumulative::IceModel_surface_flux_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "surface_ice_flux_cumulative", time_dimension_name);

  ts->set_units("kg", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "cumulative total over ice domain of top surface ice mass flux");
}

PetscErrorCode IceModel_surface_flux_cumulative::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  value = model->surface_ice_flux_cumulative;

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_grounded_basal_flux::IceModel_grounded_basal_flux(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "grounded_basal_ice_flux", time_dimension_name);

  ts->set_units("kg s-1", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "total over grounded ice domain of basal mass flux");
  ts->rate_of_change = true;
}

PetscErrorCode IceModel_grounded_basal_flux::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  value = model->grounded_basal_ice_flux_cumulative;

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_grounded_basal_flux_cumulative::IceModel_grounded_basal_flux_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "grounded_basal_ice_flux_cumulative", time_dimension_name);

  ts->set_units("kg", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "cumulative total grounded basal mass flux");
}

PetscErrorCode IceModel_grounded_basal_flux_cumulative::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  value = model->grounded_basal_ice_flux_cumulative;

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_sub_shelf_flux::IceModel_sub_shelf_flux(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "sub_shelf_ice_flux", time_dimension_name);

  ts->set_units("kg s-1", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "total sub-shelf ice flux");
  ts->rate_of_change = true;
}

PetscErrorCode IceModel_sub_shelf_flux::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  value = model->sub_shelf_ice_flux_cumulative;

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_sub_shelf_flux_cumulative::IceModel_sub_shelf_flux_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "sub_shelf_ice_flux_cumulative", time_dimension_name);

  ts->set_units("kg", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "cumulative total sub-shelf ice flux");
}

PetscErrorCode IceModel_sub_shelf_flux_cumulative::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  value = model->sub_shelf_ice_flux_cumulative;

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_nonneg_flux::IceModel_nonneg_flux(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "nonneg_flux", time_dimension_name);

  ts->set_units("kg s-1", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "'numerical' ice flux resulting from enforcing the 'thk >= 0' rule");
  ts->rate_of_change = true;
}

PetscErrorCode IceModel_nonneg_flux::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  value = model->nonneg_rule_flux_cumulative;

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_nonneg_flux_cumulative::IceModel_nonneg_flux_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "nonneg_flux_cumulative", time_dimension_name);

  ts->set_units("kg", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "cumulative 'numerical' ice flux resulting from enforcing the 'thk >= 0' rule");
}

PetscErrorCode IceModel_nonneg_flux_cumulative::update(double a, double b) {
  PetscErrorCode ierr;
  double value;

  value = model->nonneg_rule_flux_cumulative;

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_discharge_flux::IceModel_discharge_flux(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "discharge_flux", time_dimension_name);

  ts->set_units("kg s-1", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "discharge (calving & icebergs) flux");
  ts->rate_of_change = true;
}

PetscErrorCode IceModel_discharge_flux::update(double a, double b) {
  PetscErrorCode ierr;
  double value = model->discharge_flux_cumulative;

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_discharge_flux_cumulative::IceModel_discharge_flux_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "discharge_flux_cumulative", time_dimension_name);

  ts->set_units("kg", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "cumulative discharge (calving etc.) flux");
}

PetscErrorCode IceModel_discharge_flux_cumulative::update(double a, double b) {
  PetscErrorCode ierr;
  double value = model->discharge_flux_cumulative;

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_dHdt::IceModel_dHdt(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("dHdt", grid);

  set_attrs("ice thickness rate of change", "tendency_of_land_ice_thickness",
            "m s-1", "m year-1", 0);

  vars[0].set_double("valid_min",  grid.convert(-1e6, "m/year", "m/s"));
  vars[0].set_double("valid_max",  grid.convert( 1e6, "m/year", "m/s"));
  vars[0].set_double("_FillValue", grid.config.get("fill_value", "m/year", "m/s"));
  vars[0].set_string("cell_methods", "time: mean");

  last_ice_thickness.create(grid, "last_ice_thickness", WITHOUT_GHOSTS);
  last_ice_thickness.set_attrs("internal",
                               "ice thickness at the time of the last report of dHdt",
                               "m", "land_ice_thickness");

  last_report_time = GSL_NAN;
}

PetscErrorCode IceModel_dHdt::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "dHdt", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  if (gsl_isnan(last_report_time)) {
    ierr = result->set(grid.convert(2e6, "m/year", "m/s")); CHKERRQ(ierr);
  } else {

    ierr = result->begin_access(); CHKERRQ(ierr);
    ierr = last_ice_thickness.begin_access(); CHKERRQ(ierr);
    ierr = model->ice_thickness.begin_access(); CHKERRQ(ierr);

    double dt = grid.time->current() - last_report_time;
    for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
      for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
        (*result)(i, j) = (model->ice_thickness(i, j) - last_ice_thickness(i, j)) / dt;
      }
    }

    ierr = model->ice_thickness.end_access(); CHKERRQ(ierr);
    ierr = last_ice_thickness.end_access(); CHKERRQ(ierr);
    ierr = result->end_access(); CHKERRQ(ierr);

  }

  // Save the ice thickness and the corresponding time:
  ierr = this->update_cumulative(); CHKERRQ(ierr);

  output = result;
  return 0;
}

PetscErrorCode IceModel_dHdt::update_cumulative() {
  PetscErrorCode ierr;
  ierr = model->ice_thickness.begin_access(); CHKERRQ(ierr);
  ierr = last_ice_thickness.begin_access(); CHKERRQ(ierr);

  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      last_ice_thickness(i, j) = model->ice_thickness(i, j);
    }
  }

  ierr = last_ice_thickness.end_access(); CHKERRQ(ierr);
  ierr = model->ice_thickness.end_access(); CHKERRQ(ierr);

  last_report_time = grid.time->current();

  return 0;
}

IceModel_ivolg::IceModel_ivolg(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "ivolg", time_dimension_name);

  ts->set_units("m3", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "total grounded ice volume");
}

PetscErrorCode IceModel_ivolg::update(double a, double b) {
  PetscErrorCode ierr;
  double volume=0.0, value;

  MaskQuery mask(model->vMask);

  ierr = model->ice_thickness.begin_access(); CHKERRQ(ierr);
  ierr = model->vMask.begin_access(); CHKERRQ(ierr);
  ierr = model->cell_area.begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (mask.grounded_ice(i,j))
        volume += model->cell_area(i,j) * model->ice_thickness(i,j);;
    }
  }
  ierr = model->cell_area.end_access(); CHKERRQ(ierr);
  ierr = model->vMask.end_access(); CHKERRQ(ierr);
  ierr = model->ice_thickness.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&volume, &value, grid.com); CHKERRQ(ierr);

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_ivolf::IceModel_ivolf(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "ivolf", time_dimension_name);

  ts->set_units("m3", "");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "total floating ice volume");
}

PetscErrorCode IceModel_ivolf::update(double a, double b) {
  PetscErrorCode ierr;
  double volume=0.0, value;

  MaskQuery mask(model->vMask);

  ierr = model->ice_thickness.begin_access(); CHKERRQ(ierr);
  ierr = model->vMask.begin_access(); CHKERRQ(ierr);
  ierr = model->cell_area.begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (mask.floating_ice(i,j))
        volume += model->cell_area(i,j) * model->ice_thickness(i,j);;
    }
  }
  ierr = model->cell_area.end_access(); CHKERRQ(ierr);
  ierr = model->vMask.end_access(); CHKERRQ(ierr);
  ierr = model->ice_thickness.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&volume, &value, grid.com); CHKERRQ(ierr);

  ierr = ts->append(value, a, b); CHKERRQ(ierr);

  return 0;
}

//! \brief Reports the maximum horizontal absolute velocity component over the grid.
/*!
 * This is the value used by the adaptive time-stepping code in the CFL condition
 * for horizontal advection (i.e. in energy and mass conservation time steps).
 *
 * This is not the maximum horizontal speed, but rather the maximum of components.
 *
 * Note that this picks up the value computed during the time-step taken at a
 * reporting time. (It is not the "average over the reporting interval computed using
 * differencing in time", as other rate-of-change diagnostics.)
 */
IceModel_max_hor_vel::IceModel_max_hor_vel(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "max_hor_vel", time_dimension_name);

  ts->set_units("m/second", "m/year");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "maximum abs component of horizontal ice velocity over grid in last time step during time-series reporting interval");
}

PetscErrorCode IceModel_max_hor_vel::update(double a, double b) {
  PetscErrorCode ierr;

  double gmaxu = model->gmaxu, gmaxv = model->gmaxv;

  ierr = ts->append(gmaxu > gmaxv ? gmaxu : gmaxv, a, b); CHKERRQ(ierr);

  return 0;
}

IceModel_H_to_Href_flux::IceModel_H_to_Href_flux(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {
  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "H_to_Href_flux", time_dimension_name);

  ts->set_units("kg s-1", "kg s-1");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "mass flux from thk to Href");
  ts->rate_of_change = true;
}

PetscErrorCode IceModel_H_to_Href_flux::update(double a, double b) {
  PetscErrorCode ierr;

  ierr = ts->append(model->H_to_Href_flux_cumulative, a, b); CHKERRQ(ierr);

  return 0;
}


IceModel_Href_to_H_flux::IceModel_Href_to_H_flux(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {
  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "Href_to_H_flux", time_dimension_name);

  ts->set_units("kg s-1", "kg s-1");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "mass flux from Href to thk");
  ts->rate_of_change = true;
}

PetscErrorCode IceModel_Href_to_H_flux::update(double a, double b) {
  PetscErrorCode ierr;

  ierr = ts->append(model->Href_to_H_flux_cumulative, a, b); CHKERRQ(ierr);

  return 0;
}



IceModel_sum_divQ_flux::IceModel_sum_divQ_flux(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMTSDiag<IceModel>(m, g, my_vars) {
  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "sum_divQ_flux", time_dimension_name);

  ts->set_units("kg s-1", "kg s-1");
  ts->set_dimension_units(time_units, "");
  ts->set_attr("long_name", "sum(divQ)");
  ts->rate_of_change = true;
}

PetscErrorCode IceModel_sum_divQ_flux::update(double a, double b) {
  PetscErrorCode ierr;

  ierr = ts->append(model->sum_divQ_SIA_cumulative + model->sum_divQ_SSA_cumulative,
                    a, b); CHKERRQ(ierr);

  return 0;
}


IceModel_nonneg_flux_2D_cumulative::IceModel_nonneg_flux_2D_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("nonneg_flux_cumulative", grid);

  set_attrs("cumulative non-negative rule (thk >= 0) flux (positive means ice gain)",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
}

PetscErrorCode IceModel_nonneg_flux_2D_cumulative::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "nonneg_flux_cumulative", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  ierr = result->copy_from(model->nonneg_flux_2D_cumulative); CHKERRQ(ierr);

  output = result;
  return 0;
}

IceModel_grounded_basal_flux_2D_cumulative::IceModel_grounded_basal_flux_2D_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("grounded_basal_flux_cumulative", grid);

  set_attrs("cumulative grounded basal flux (positive means ice gain)",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
}

PetscErrorCode IceModel_grounded_basal_flux_2D_cumulative::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "grounded_basal_flux_cumulative", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  ierr = result->copy_from(model->grounded_basal_flux_2D_cumulative); CHKERRQ(ierr);

  output = result;
  return 0;
}

IceModel_floating_basal_flux_2D_cumulative::IceModel_floating_basal_flux_2D_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("floating_basal_flux_cumulative", grid);

  set_attrs("cumulative floating basal flux (positive means ice gain)",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
}

PetscErrorCode IceModel_floating_basal_flux_2D_cumulative::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "floating_basal_flux_cumulative", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  ierr = result->copy_from(model->floating_basal_flux_2D_cumulative); CHKERRQ(ierr);

  output = result;
  return 0;
}


IceModel_discharge_flux_2D_cumulative::IceModel_discharge_flux_2D_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("discharge_flux_cumulative", grid);

  set_attrs("cumulative ice discharge (calving) flux (negative means ice loss)",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
}

PetscErrorCode IceModel_discharge_flux_2D_cumulative::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "discharge_flux_cumulative", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  ierr = result->copy_from(model->discharge_flux_2D_cumulative); CHKERRQ(ierr);

  output = result;
  return 0;
}

#if (PISM_USE_PROJ4==1)
IceModel_lat_lon_bounds::IceModel_lat_lon_bounds(IceModel *m, IceGrid &g, PISMVars &my_vars,
                                                 std::string var_name, std::string proj_string)
  :PISMDiag<IceModel>(m, g, my_vars) {
  assert(var_name == "lat" || var_name == "lon");
  m_var_name = var_name;

  // set metadata:
  std::vector<double> levels(4);
  for (int k = 0; k < 4; ++k)
    levels[k] = k;

  vars[0].init_3d(m_var_name + "_bounds", grid, levels);
  vars[0].get_z().set_name("grid_corners");
  vars[0].get_z().clear_all_strings();
  vars[0].get_z().clear_all_doubles();
  vars[0].set_time_independent(true);

  if (m_var_name == "lon") {
    set_attrs("longitude bounds", "", "degree_east", "degree_east", 0);
    vars[0].set_double("valid_min", -180);
    vars[0].set_double("valid_max", 180);
  } else {
    set_attrs("latitude bounds", "", "degree_north", "degree_north", 0);
    vars[0].set_double("valid_min", -90);
    vars[0].set_double("valid_max", 90);
  }

  lonlat = pj_init_plus("+proj=latlong +datum=WGS84 +ellps=WGS84");
  if (lonlat == NULL) {
    PetscPrintf(grid.com, "PISM ERROR: projection initialization failed "
                "('+proj=latlong +datum=WGS84 +ellps=WGS84').\n");
    PISMEnd();
  }

  pism = pj_init_plus(proj_string.c_str());
  if (pism == NULL) {
    PetscPrintf(grid.com, "PISM ERROR: proj.4 string '%s' is invalid. Exiting...\n",
                proj_string.c_str());
    PISMEnd();
  }

}

IceModel_lat_lon_bounds::~IceModel_lat_lon_bounds() {
  pj_free(pism);
  pj_free(lonlat);
}

PetscErrorCode IceModel_lat_lon_bounds::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3Custom *result = new IceModelVec3Custom;
  std::map<std::string,std::string> attrs;
  std::vector<double> indices(4);

  ierr = result->create(grid, m_var_name + "_bounds", "grid_corners",
                        indices, attrs); CHKERRQ(ierr);
  result->metadata() = vars[0];

  double dx2 = 0.5 * grid.dx, dy2 = 0.5 * grid.dy;
  double x_offsets[] = {-dx2, dx2, dx2, -dx2};
  double y_offsets[] = {-dy2, -dy2, dy2, dy2};

  bool latitude = true;
  if (m_var_name == "lon")
    latitude = false;

  ierr =  result->begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      double x0 = grid.x[i], y0 = grid.y[j];
      double *values;

      ierr = result->getInternalColumn(i,j,&values); CHKERRQ(ierr);

      for (int k = 0; k < 4; ++k) {
        double
          x = x0 + x_offsets[k],
          y = y0 + y_offsets[k];

        // compute lon,lat coordinates:
        pj_transform(pism, lonlat, 1, 1, &x, &y, NULL);

        // NB! proj.4 converts x,y pairs into lon,lat pairs in *radians*.

        if (latitude)
          values[k] = y * RAD_TO_DEG;
        else
          values[k] = x * RAD_TO_DEG;
      }
    }
  }
  ierr =  result->end_access(); CHKERRQ(ierr);

  output = result;

  return 0;
}
#elif (PISM_USE_PROJ4==0)
  // do nothing
#else  // PISM_USE_PROJ4 is not set
#error "PISM build system error: PISM_USE_PROJ4 is not set."
#endif
