// Copyright (C) 2012-2015 Ricarda Winkelmann, Ronja Reese and Torsten Albrecht
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

#include "POoceanboxmodel.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
//#include "pism_options.hh" //05
#include "PISMConfig.hh" //06




/*!
 The aim of these routines is to simulate the ocean circulation underneath the ice shelves and compute the basal melt/refreezing rates according to the ocean box model described in olbers_hellmer10.
*/




POoceanboxmodel::POBMConstants::POBMConstants(const PISMConfig &config) {
 
  PetscErrorCode ierr;

  numberOfBasins = 20;

  continental_shelf_depth = -800;

  T_dummy = -1.5; //FIXME why these values? 
  S_dummy = 34.5; //

  earth_grav = config.get("standard_gravity");
  rhoi = config.get("ice_density");
  rhow = config.get("sea_water_density");
  rho_star  = 1033;         // kg/m^3
  nu        = rhoi / rho_star;       // no unit

  latentHeat = config.get("water_latent_heat_fusion");
  c_p_ocean = 3974.0;       // J/(K*kg), specific heat capacity of ocean mixed layer
  lambda    = latentHeat / c_p_ocean;   // °C, NOTE K vs °C

  a         = -0.057;       // °C/psu
  b         = 0.0832;       // °C
  c         = 7.64e-4;      // °C/dbar

  alpha     = 7.5e-5;       // 1/°C, NOTE K vs °C
  beta      = 7.7e-4;       // 1/psu

  gamma_T = 1e-6;
  value_C = 5e6;

  // other ice shelves
  gamma_T_o   = 1.0e-4; //config.get("gamma_T"); //1e-4;     
  // m/s, thermal exchange velocity for Beckmann-Goose parameterization
  meltFactor = 0.002;     // FIXME!!!! (Wrong value!) FIXME config
  meltSalinity = 35.0;
  b2 = 0.0939;

 }


const int POoceanboxmodel::box_unidentified = -99;     // This should never show up in the .nc-files.
const int POoceanboxmodel::box_neighboring = -1; // This should never show up in the .nc-files.
const int POoceanboxmodel::box_noshelf = 0;
const int POoceanboxmodel::box_GL = 1;  // ocean box covering the grounding line region
const int POoceanboxmodel::box_IF = 2;  // ocean box covering the rest of the ice shelf
const int POoceanboxmodel::box_other = 3;  // ice_shelf but there is no GL_box in the corresponding basin

const int POoceanboxmodel::maskfloating = MASK_FLOATING;  
const int POoceanboxmodel::maskocean = MASK_ICE_FREE_OCEAN;  
const int POoceanboxmodel::maskgrounded = MASK_GROUNDED;  

const int POoceanboxmodel::imask_inner = 2;
const int POoceanboxmodel::imask_outer = 0;
const int POoceanboxmodel::imask_exclude = 1;
const int POoceanboxmodel::imask_unidentified = -1;



POoceanboxmodel::POoceanboxmodel(IceGrid &g, const PISMConfig &conf)
  : PGivenClimate<POModifier,PISMOceanModel>(g, conf, NULL)
{
  PetscErrorCode ierr = allocate_POoceanboxmodel(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();
}

POoceanboxmodel::~POoceanboxmodel() {
  // empty
}

PetscErrorCode POoceanboxmodel::allocate_POoceanboxmodel() {
  PetscErrorCode ierr;
  option_prefix   = "-ocean_oceanboxmodel";

  // will be de-allocated by the parent's destructor
  theta_ocean    = new IceModelVec2T;
  salinity_ocean = new IceModelVec2T;

  m_fields["theta_ocean"]     = theta_ocean;
  m_fields["salinity_ocean"]  = salinity_ocean;

  ierr = process_options(); CHKERRQ(ierr);

  std::map<std::string, std::string> standard_names;
  ierr = set_vec_parameters(standard_names); CHKERRQ(ierr);

  ierr = theta_ocean->create(grid, "theta_ocean", false); CHKERRQ(ierr);
  ierr = theta_ocean->set_attrs("climate_forcing",
                                "absolute potential temperature of the adjacent ocean",
                                "Kelvin", ""); CHKERRQ(ierr);

  ierr = salinity_ocean->create(grid, "salinity_ocean", false); CHKERRQ(ierr);
  ierr = salinity_ocean->set_attrs("climate_forcing",
                                   "salinity of the adjacent ocean",
                                   "g/kg", ""); CHKERRQ(ierr);

  ierr = shelfbtemp.create(grid, "shelfbtemp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = shelfbtemp.set_attrs("climate_forcing",
                              "absolute temperature at ice shelf base",
                              "Kelvin", ""); CHKERRQ(ierr);

  ierr = shelfbmassflux.create(grid, "shelfbmassflux", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = shelfbmassflux.set_attrs("climate_forcing",
                                  "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                                  "kg m-2 s-1", ""); CHKERRQ(ierr);
  ierr = shelfbmassflux.set_glaciological_units("kg m-2 year-1"); CHKERRQ(ierr);
  shelfbmassflux.write_in_glaciological_units = true;


  //////////////////////////////////////////////////////////////////////////


  // mask to identify the ocean boxes
  ierr = BOXMODELmask.create(grid, "BOXMODELmask", WITH_GHOSTS); CHKERRQ(ierr);
  ierr = BOXMODELmask.set_attrs("model_state", "mask displaying ocean box model grid","", ""); CHKERRQ(ierr);

  // mask to identify the grounded ice rises
  ierr = ICERISESmask.create(grid, "ICERISESmask", WITH_GHOSTS); CHKERRQ(ierr);
  ierr = ICERISESmask.set_attrs("model_state", "mask displaying ice rises","", ""); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-exclude_icerises", exicerises_set); CHKERRQ(ierr);

  // mask displaying continental shelf - region where mean salinity and ocean temperature is calculated
  ierr = OCEANMEANmask.create(grid, "OCEANMEANmask", WITH_GHOSTS); CHKERRQ(ierr);
  ierr = OCEANMEANmask.set_attrs("model_state", "mask displaying ocean region for parameter input","", ""); CHKERRQ(ierr);


  // salinity
  ierr = Soc.create(grid, "Soc", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = Soc.set_attrs("model_state", "ocean salinity field","", "ocean salinity field"); CHKERRQ(ierr);  //NOTE unit=psu

  ierr = Soc_base.create(grid, "Soc_base", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = Soc_base.set_attrs("model_state", "ocean base salinity field","", "ocean base salinity field"); CHKERRQ(ierr);  //NOTE unit=psu


  // temperature
  ierr = Toc.create(grid, "Toc", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = Toc.set_attrs("model_state", "ocean temperature field","K", "ocean temperature field"); CHKERRQ(ierr);

  ierr = Toc_base.create(grid, "Toc_base", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = Toc_base.set_attrs("model_state", "ocean base temperature","K", "ocean base temperature"); CHKERRQ(ierr);

  ierr = Toc_inCelsius.create(grid, "Toc_inCelsius", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = Toc_inCelsius.set_attrs("model_state", "ocean box model temperature field","degree C", "ocean box model temperature field"); CHKERRQ(ierr);

  ierr = T_star.create(grid, "T_star", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = T_star.set_attrs("model_state", "T_star field","degree C", "T_star field"); CHKERRQ(ierr);

  ierr = Toc_anomaly.create(grid, "Toc_anomaly", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = Toc_anomaly.set_attrs("model_state", "ocean temperature anomaly","K", "ocean temperature anomaly"); CHKERRQ(ierr);


  // overturning rate
  ierr = overturning.create(grid, "overturning", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = overturning.set_attrs("model_state", "cavity overturning","m^3 s-1", "cavity overturning"); CHKERRQ(ierr); // no CF standard_name ??

  // heat flux
  ierr = heatflux.create(grid, "ocean heat flux", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = heatflux.set_attrs("climate_state", "ocean heat flux", "W/m^2", ""); CHKERRQ(ierr);

  // basal melt rate
  ierr = basalmeltrate_shelf.create(grid, "basal melt rate from ocean box model", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.set_attrs("climate_state", "basal melt rate from ocean box model", "m/s", ""); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  basalmeltrate_shelf.write_in_glaciological_units = true;

 ///////// forcing  /////////////////////////////////////////////////////////////////////////////////


  // option for scalar forcing of ocean temperature
  ierr = PISMOptionsIsSet("-ocean_obm_deltaT", ocean_oceanboxmodel_deltaT_set); CHKERRQ(ierr);

  if (ocean_oceanboxmodel_deltaT_set) {
    bool delta_T_set;
    std::string delta_T_file;

    ierr = PISMOptionsString("-ocean_obm_deltaT",
                             "Specifies the ocean temperature offsets file to use with -ocean_obm_deltaT",
                             delta_T_file, delta_T_set); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com, 
                      "  reading delta_T data from forcing file %s for -ocean_obm_deltaT actions ...\n",
                      delta_T_file.c_str());  CHKERRQ(ierr);

    delta_T = new Timeseries(&grid, "delta_T",grid.config.get_string("time_dimension_name"));
    ierr = delta_T->set_units("Kelvin", ""); CHKERRQ(ierr);
    ierr = delta_T->set_dimension_units(grid.time->units_string(), ""); CHKERRQ(ierr);
    ierr = delta_T->set_attr("long_name", "ocean temperature offsets"); CHKERRQ(ierr);
    //ierr = delta_T->read(delta_T_file, grid.time->use_reference_date()); CHKERRQ(ierr);

    PIO nc(grid.com, "netcdf3", grid.get_unit_system());
    ierr = nc.open(delta_T_file, PISM_NOWRITE); CHKERRQ(ierr);
    {
      ierr = delta_T->read(nc, grid.time); CHKERRQ(ierr);
    }
    ierr = nc.close(); CHKERRQ(ierr);

    bool delta_T_factor_set;
    delta_T_factor=1.0;

    ierr = PISMOptionsReal("-ocean_obm_factor","ocean_obm_factor set",delta_T_factor, delta_T_factor_set); CHKERRQ(ierr);

  }
  return 0;
}


PetscErrorCode POoceanboxmodel::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the ocean box model (based on Olbers & Hellmer (2010)...\n"); CHKERRQ(ierr);

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (ice_thickness == NULL) {SETERRQ(grid.com, 1, "ERROR: ice thickness is not available");}

  topg = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (topg == NULL) SETERRQ(grid.com, 1, "ERROR: bedrock topography is not available");

  mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  if (!mask) { SETERRQ(grid.com, 1, "ERROR: mask is not available"); }

  basins = dynamic_cast<IceModelVec2S*>(vars.get("drainage_basins")); //if option drainageBasins set
  if (!basins) { SETERRQ(grid.com, 1, "ERROR: drainage basins is not available"); }

  bool omeans_set;
  ierr = PISMOptionsIsSet("-ocean_means", omeans_set); CHKERRQ(ierr);

  //FIXME: not necessarry when -ocean_means set
  ierr = theta_ocean->init(filename, bc_period, bc_reference_time); CHKERRQ(ierr);
  ierr = salinity_ocean->init(filename, bc_period, bc_reference_time); CHKERRQ(ierr);

  // read time-independent data right away:
  if (theta_ocean->get_n_records() == 1 && salinity_ocean->get_n_records() == 1) {
    ierr = update(grid.time->current(), 0); CHKERRQ(ierr); // dt is irrelevant
  }

  POBMConstants cc(config);
  ierr = initBasinsOptions(cc); CHKERRQ(ierr);


  return 0;
}

void POoceanboxmodel::add_vars_to_output(std::string keyword, std::set<std::string> &result) {
  PGivenClimate<POModifier,PISMOceanModel>::add_vars_to_output(keyword, result);

  if (keyword != "none" && keyword != "small") {
    result.insert("shelfbtemp");
    result.insert("shelfbmassflux");
  }
}


PetscErrorCode POoceanboxmodel::define_variables(std::set<std::string> vars, const PIO &nc,
                                           PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  ierr = PGivenClimate<POModifier,PISMOceanModel>::define_variables(vars, nc, nctype); CHKERRQ(ierr);

  if (set_contains(vars, "shelfbtemp")) {
    ierr = shelfbtemp.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    ierr = shelfbmassflux.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode POoceanboxmodel::initBasinsOptions(const POBMConstants &cc) {
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"0b : set number of Basins\n"); CHKERRQ(ierr);

  //set number of basins per option
  bool number_of_basins_set;
  numberOfBasins = cc.numberOfBasins;
  ierr = PISMOptionsInt("-number_of_basins","Number of Drainage Basins",numberOfBasins, number_of_basins_set); CHKERRQ(ierr);


  Toc_base_vec.resize(numberOfBasins);
  Soc_base_vec.resize(numberOfBasins);
  gamma_T_star_vec.resize(numberOfBasins);
  C_vec.resize(numberOfBasins);

  counter.resize(numberOfBasins);
  counter_GLbox.resize(numberOfBasins);
  counter_CFbox.resize(numberOfBasins);
  //k_basins.resize(numberOfBasins);;

  mean_salinity_GLbox_vector.resize(numberOfBasins);
  mean_meltrate_GLbox_vector.resize(numberOfBasins);
  mean_overturning_GLbox_vector.resize(numberOfBasins);


  //set gamma_T and value_C per option
  bool gamma_T_set;
  gamma_T = cc.gamma_T;
  ierr = PISMOptionsReal("-gamma_T","-gamma_T",gamma_T, gamma_T_set); CHKERRQ(ierr);

  bool value_C_set;
  value_C = cc.value_C;
  ierr = PISMOptionsReal("-value_C","-value_C",value_C, value_C_set); CHKERRQ(ierr);

  ///////////////////////////////////////////////////////////////////////////////////
  // data have been calculated previously for the 18 Rignot basins
  //const double Toc_base_schmidtko[18] = {0.0,271.56415203,271.63356482,271.42074233,271.46720524,272.30038929,271.52821139,271.5440751,271.58201494,272.90159695,273.61058862,274.19203524,274.32083917,272.55938554,271.35349906,271.39337366,271.49926019,271.49473924}; //Schmidtko
  //const double Soc_base_schmidtko[18] = {0.0,34.49308909,34.50472554,34.70187911,34.65306507,34.7137078,34.74817136,34.89206844,34.78056731,34.60528314,34.72521443,34.86210624,34.83836297,34.73392011,34.83617088,34.82137147,34.69477334,34.48145265}; //Schmidtko

  // data have been calculated previously for the 20 Zwally basins
  const double Toc_base_schmidtko[20] = {0.0,271.39431005,271.49081157,271.49922596,271.56714804,271.63507013,271.42228667,271.46720524,272.42253843,271.53779093,271.84942002,271.31676801,271.56846696,272.79372542,273.61694268,274.19168456,274.31958227,273.38372579,271.91951514,271.35349906}; //Schmidtko

  const double Soc_base_schmidtko[20] = {0.0,34.82193374,34.69721226,34.47641407,34.48950162,34.50258917,34.70101507,34.65306507,34.73295029,34.74859586,34.8368573,34.9529016,34.79486795,34.58380953,34.7260615,34.86198383,34.8374212 ,34.70418016,34.75598208,34.83617088}; //Schmidtko


  //const double Toc_base_woa[18] = {0.0,272.28351693,272.10101401,271.65965597,271.50766979,273.02732277,272.12473624,271.79505722,271.93548261,273.37866926,272.98126135,273.73564726,273.95971315,273.02383769,272.56732024,271.75152607,271.93962932,272.46601985}; //World Ocean Atlas
  //const double Soc_base_woa[18] = {0.0,34.48230812,34.46499742,34.51939747,34.40695979,34.62947893,34.59932424,34.77118004,34.71666183,34.54603987,34.44824601,34.62416923,34.57034648,34.60459029,34.68356516,34.67159838,34.62308218,34.49961882}; //World Ocean Atlas

  const double Toc_base_woa[20] = {272.99816667,271.27814004,272.1840257,272.04435251,272.20415662,272.36396072,271.48763831,271.99695864,272.06504052,272.27114732,272.66657018,271.18920729,271.74067699,273.01811291,272.15295572,273.08542047,272.74584469,273.14263356,272.58496563,272.45217911}; //World Ocean Atlas
  const double Soc_base_woa[20] = {34.6810522,34.78161073,34.67151084,34.66538478,34.67127468,34.67716458,34.75327377,34.69213327,34.72086382,34.70670158,34.71210592,34.80229468,34.76588022,34.69745763,34.7090778,34.68690903,34.66379606,34.64572337,34.6574402,34.65813983}; //World Ocean Atlas


  bool ocean_mean_set;
  std::string data_input;
  ierr = PISMOptionsString("-ocean_means", "Input data name",
                             data_input, ocean_mean_set); CHKERRQ(ierr);

  /////////////////////////////////////////////////////////////////////////////////////

  for(int k=0;k<numberOfBasins;k++) {
      if (ocean_mean_set){
        if (data_input=="schmidtko"){
          Toc_base_vec[k] = Toc_base_schmidtko[k] - 273.15;
          Soc_base_vec[k] = Soc_base_schmidtko[k];}
        else if (data_input=="woa"){
          Toc_base_vec[k] = Toc_base_woa[k] - 273.15;
          Soc_base_vec[k] = Soc_base_woa[k];}
        else{
          Toc_base_vec[k] = cc.T_dummy; //dummy, FIXME why these values? 
          Soc_base_vec[k] = cc.S_dummy; //dummy
        }
      }
      else{
        Toc_base_vec[k] = cc.T_dummy; //dummy, FIXME why these values? 
        Soc_base_vec[k] = cc.S_dummy; //dummy
      }

      gamma_T_star_vec[k]= gamma_T; 
      C_vec[k]           = value_C;
  }

  ierr = verbPrintf(5, grid.com,"     Using %d drainage basins and default values: \n"  
                                "     gamma_T_star= %.2e, C = %.2e... \n"  
                                 , numberOfBasins, gamma_T, value_C); CHKERRQ(ierr);

  if (!ocean_mean_set) {
    ierr = verbPrintf(5, grid.com,"  calculate Soc and Toc from thetao and salinity... \n"); CHKERRQ(ierr);       
  
    // set continental shelf depth // FIXME -800 might be too high for Antacrtica
    continental_shelf_depth = cc.continental_shelf_depth;
    ierr = PISMOptionsReal("-continental_shelf_depth","-continental_shelf_depth", continental_shelf_depth, continental_shelf_depth_set); CHKERRQ(ierr);
    if (continental_shelf_depth_set) {
      ierr = verbPrintf(5, grid.com,
                        "  Depth of continental shelf for computation of temperature and salinity input\n"
                        "  is set for whole domain to continental_shelf_depth=%.0f meter\n", continental_shelf_depth); CHKERRQ(ierr);
    }
  }

  return 0;
}



PetscErrorCode POoceanboxmodel::update(double my_t, double my_dt) {

  // Make sure that sea water salinity and sea water potential
  // temperature fields are up to date:
  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  ierr = theta_ocean->average(m_t, m_dt); CHKERRQ(ierr);
  ierr = salinity_ocean->average(m_t, m_dt); CHKERRQ(ierr);

  if ((delta_T != NULL) && ocean_oceanboxmodel_deltaT_set) {
    temp_anomaly = (*delta_T)(my_t + 0.5*my_dt);
    ierr = verbPrintf(4, grid.com,"0a : set global temperature anomaly = %.3f\n",temp_anomaly); CHKERRQ(ierr);
  }

  bool omeans_set;
  ierr = PISMOptionsIsSet("-ocean_means", omeans_set); CHKERRQ(ierr);

  POBMConstants cc(config);
  ierr = initBasinsOptions(cc); CHKERRQ(ierr);

  ierr = roundBasins(); CHKERRQ(ierr); 
  if (omeans_set){
    ierr = verbPrintf(4, grid.com,"0c : reading mean salinity and temperatures\n"); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(4, grid.com,"0c : calculating mean salinity and temperatures\n"); CHKERRQ(ierr);
    ierr = identifyMASK(OCEANMEANmask,"ocean"); CHKERRQ(ierr);
    ierr = computeOCEANMEANS(); CHKERRQ(ierr);   
  }


  //geometry of ice shelves and temperatures
  ierr = verbPrintf(4, grid.com,"A  : calculating shelf_base_temperature\n"); CHKERRQ(ierr);
  if (exicerises_set) {
    ierr = identifyMASK(ICERISESmask,"icerises");}
  ierr = extentOfIceShelves(); CHKERRQ(ierr); 
  ierr = identifyBOXMODELmask(); CHKERRQ(ierr); 
  ierr = oceanTemperature(cc); CHKERRQ(ierr);
  ierr = Toc.copy_to(shelfbtemp); CHKERRQ(ierr);


  //basal melt rates underneath ice shelves
  ierr = verbPrintf(4, grid.com,"B  : calculating shelf_base_mass_flux\n"); CHKERRQ(ierr);
  ierr = basalMeltRateForGroundingLineBox(cc); CHKERRQ(ierr);
  ierr = basalMeltRateForIceFrontBox(cc); CHKERRQ(ierr); // TODO Diese Routinen woanders aufrufen (um Dopplung zu vermeiden)
  ierr = basalMeltRateForOtherShelves(cc); CHKERRQ(ierr);  //Assumes that mass flux is proportional to the shelf-base heat flux.
  //const double secpera=31556926.0;
  ierr = basalmeltrate_shelf.scale(cc.rhoi); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.copy_to(shelfbmassflux); CHKERRQ(ierr); //TODO Check if scaling with ice density


  return 0;
}



//! Round basin mask non integer values to an integral value of the next neigbor
PetscErrorCode POoceanboxmodel::roundBasins() {
  PetscErrorCode ierr;
  //FIXME: THIS routine should be applied once in init, and roundbasins should be stored as field

  ierr = PISMOptionsIsSet("-round_basins", roundbasins_set); CHKERRQ(ierr);

  ierr = basins->begin_access();   CHKERRQ(ierr);

  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {

      double id_fractional = (*basins)(i,j),
            id_fr_ne = (*basins)(i+1,j+1),
            id_fr_nw = (*basins)(i-1,j+1),
            id_fr_sw = (*basins)(i-1,j-1),
            id_fr_se = (*basins)(i+1,j-1);

      int  id_rounded = static_cast<int>(round(id_fractional)),
                id_ro_ne = static_cast<int>(round(id_fr_ne)),
                id_ro_nw = static_cast<int>(round(id_fr_nw)),
                id_ro_sw = static_cast<int>(round(id_fr_sw)),
                id_ro_se = static_cast<int>(round(id_fr_se)),
                id = -1;
  
      if (roundbasins_set){

        if( PetscAbs(id_fractional - static_cast<float>(id_rounded)) > 0.0){ //if id_fractional differs from integer value  
    
          if (id_fr_sw == static_cast<float>(id_ro_sw) && id_fr_sw != 0){ 
            id = id_ro_sw;
          } else if (id_fr_se == static_cast<float>(id_ro_se) && id_fr_se != 0){ 
            id = id_ro_se;
          } else if (id_fr_nw == static_cast<float>(id_ro_nw) && id_fr_nw != 0){ 
            id = id_ro_nw;
          } else if (id_fr_ne == static_cast<float>(id_ro_ne) && id_fr_ne != 0){ 
            id = id_ro_ne;
          } else { //if no neigbour has an integer id
            id = id_rounded;
          }
        } else { // if id_rounded == id_fractional
          id = id_rounded;
        }
      } else { // if -round_basins not set
        id = id_rounded;
      }
      (*basins)(i,j)=id;
    }
  }

  ierr = basins->end_access();   CHKERRQ(ierr); 

  return 0;
}


//! Identify 
//!   ocean: identify ocean up to continental shelf without detached submarine islands regions
//!   icerises: identify grounded regions without detached ice rises

PetscErrorCode POoceanboxmodel::identifyMASK(IceModelVec2S &inputmask, std::string masktype) {
  
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"0b1: in identifyMASK rountine\n"); CHKERRQ(ierr);

  int seed_x = (grid.Mx - 1)/2,
           seed_y = (grid.My - 1)/2;

  double linner_identified = 0.0,
              all_inner_identified = 1.0,
              previous_step_identified = 0.0;;


  ierr = inputmask.begin_access();   CHKERRQ(ierr);
  inputmask.set(imask_unidentified);
  if ((seed_x >= grid.xs) && (seed_x < grid.xs+grid.xm) && (seed_y >= grid.ys)&& (seed_y <= grid.ys+grid.ym)){
    inputmask(seed_x,seed_y)=imask_inner;
  }
  ierr = inputmask.end_access();   CHKERRQ(ierr);

  int iteration_round = 0;
  // find inner region first
  while(all_inner_identified > previous_step_identified){

    iteration_round+=1;
    previous_step_identified = all_inner_identified;

    ierr = inputmask.begin_access();   CHKERRQ(ierr);
    ierr = mask->begin_access();   CHKERRQ(ierr);
    ierr = topg->begin_access(); CHKERRQ(ierr); 

    for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {

        bool masktype_condition = false;
        if (masktype=="ocean"){
          masktype_condition = ((*mask)(i,j)!=maskocean || (*topg)(i,j) >= continental_shelf_depth);}
        else if (masktype=="icerises"){
          masktype_condition = ((*mask)(i,j)==maskgrounded);
        }

        if (masktype_condition && inputmask(i,j)==imask_unidentified &&
          (inputmask(i,j+1)==imask_inner || inputmask(i,j-1)==imask_inner || 
           inputmask(i+1,j)==imask_inner || inputmask(i-1,j)==imask_inner)){
           inputmask(i,j)=imask_inner;
           linner_identified+=1;
        }
        else if (masktype_condition == false){
          inputmask(i,j)=imask_outer;
        }

        //ierr = verbPrintf(4, grid.com,"!!! %d %d, %.0f \n",i,j,inputmask(i,j)); CHKERRQ(ierr);
      }
    }

    ierr = mask->end_access();   CHKERRQ(ierr);   
    ierr = topg->end_access(); CHKERRQ(ierr);  

    ierr = inputmask.end_access();   CHKERRQ(ierr);
    //ierr = inputmask.beginGhostComm(); CHKERRQ(ierr);
    //ierr = inputmask.endGhostComm(); CHKERRQ(ierr);
    ierr = inputmask.update_ghosts(); CHKERRQ(ierr);

    ierr = PISMGlobalSum(&linner_identified, &all_inner_identified, grid.com); CHKERRQ(ierr);
  
  }

  //set value for excluded areas (ice rises or submarine islands)
  ierr = inputmask.begin_access();   CHKERRQ(ierr);
  ierr = mask->begin_access();   CHKERRQ(ierr);

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {

      if (inputmask(i,j)==imask_unidentified){
        inputmask(i,j)=imask_exclude;
      }

      if (masktype=="ocean"){ //exclude ice covered parts
        if ((*mask)(i,j)!=maskocean && inputmask(i,j) == imask_inner){
          inputmask(i,j) = imask_outer; 
        }
      }
    }
  }

  ierr = inputmask.end_access();   CHKERRQ(ierr); 
  ierr = mask->end_access();   CHKERRQ(ierr);   


  return 0;
}



//! When ocean_given is set compute mean salinity and temperature in each basin.
PetscErrorCode POoceanboxmodel::computeOCEANMEANS() {
  // FIXME currently the mean is also calculated over submarine islands which are higher than continental_shelf_depth

  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"0b2: in computeOCEANMEANS routine \n"); CHKERRQ(ierr);

  double lm_count[numberOfBasins]; //count cells to take mean over for each basin
  double m_count[numberOfBasins];
  double lm_Sval[numberOfBasins]; //add salinity for each basin
  double lm_Tval[numberOfBasins]; //add temperature for each basin
  double m_Tval[numberOfBasins]; 
  double m_Sval[numberOfBasins]; 

  for(int k=0;k<numberOfBasins;k++){
    m_count[k]=0.0;
    lm_count[k]=0.0;
    lm_Sval[k]=0.0;
    lm_Tval[k]=0.0;
    m_Tval[k]=0.0; 
    m_Sval[k]=0.0;
  }

  ierr = OCEANMEANmask.begin_access();   CHKERRQ(ierr); 
  ierr = theta_ocean->begin_access();   CHKERRQ(ierr); 
  ierr = salinity_ocean->begin_access();   CHKERRQ(ierr); //salinitiy
  ierr = basins->begin_access();   CHKERRQ(ierr);

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (OCEANMEANmask(i,j) == imask_inner ){
        int shelf_id =(*basins)(i,j);
        lm_count[shelf_id]+=1;
        lm_Sval[shelf_id]+=(*salinity_ocean)(i,j);
        lm_Tval[shelf_id]+=(*theta_ocean)(i,j);
      } //if
    } //j
  } //i
  
  ierr = OCEANMEANmask.end_access();   CHKERRQ(ierr); 
  ierr = theta_ocean->end_access();   CHKERRQ(ierr);  
  ierr = salinity_ocean->end_access();   CHKERRQ(ierr);  //salinity
  ierr = basins->end_access();   CHKERRQ(ierr);

  
  for(int k=0;k<numberOfBasins;k++) {
    ierr = PISMGlobalSum(&lm_count[k], &m_count[k], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lm_Sval[k], &m_Sval[k], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lm_Tval[k], &m_Tval[k], grid.com); CHKERRQ(ierr);
    
    if(k>0 && m_count[k]==0){ //if basin is not dummy basin 0 or there are no ocean cells in this basin to take the mean over.
      ierr = verbPrintf(2, grid.com,"PISM_WARNING: basin %d contains no ocean mean cells, no mean salinity or temperature values are computed! \n ", k); CHKERRQ(ierr);   
    } else {
      m_Sval[k] = m_Sval[k] / m_count[k];
      m_Tval[k] = m_Tval[k] / m_count[k]; 
    
      Toc_base_vec[k]=m_Tval[k] - 273.15;
      Soc_base_vec[k]=m_Sval[k];
      ierr = verbPrintf(4, grid.com,"  %d: temp =%.3f, salinity=%.3f\n", k, Toc_base_vec[k], Soc_base_vec[k]); CHKERRQ(ierr);
    }
  }

  return 0;
}



//! Compute the extent of the ice shelves of each basin/region (i.e. counter) and find grounding line and calving front in BOXMODELmask

PetscErrorCode POoceanboxmodel::extentOfIceShelves() {
  
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"A1b: in extent of ice shelves rountine\n"); CHKERRQ(ierr);

  double lcounter_box_unidentified = 0; //count the total amount of unidentified shelf boxes.
  double lcounter[numberOfBasins];
  double lcounter_CFbox[numberOfBasins];
  double lcounter_GLbox[numberOfBasins];

  for (int k=0;k<numberOfBasins;k++){ 
    lcounter[k]=0.0; 
    lcounter_CFbox[k]=0.0;
    lcounter_GLbox[k]=0.0;
  }
  
  ierr = mask->begin_access();   CHKERRQ(ierr); 
  ierr = basins->begin_access();   CHKERRQ(ierr); 
  ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);
  if (exicerises_set) { ierr = ICERISESmask.begin_access(); CHKERRQ(ierr);}

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {

      if ((*mask)(i,j)==maskfloating){ //if this is a ice shelf cell
        int shelf_id =(*basins)(i,j);
        lcounter[shelf_id]++;

        bool neighbor_to_land;
        if (exicerises_set) {
          neighbor_to_land = (  ICERISESmask(i,j+1)==imask_inner || ICERISESmask(i,j-1)==imask_inner || 
                                ICERISESmask(i+1,j)==imask_inner || ICERISESmask(i-1,j)==imask_inner || 
                                ICERISESmask(i+1,j+1)==imask_inner || ICERISESmask(i+1,j-1)==imask_inner ||
                                ICERISESmask(i-1,j+1)==imask_inner || ICERISESmask(i-1,j-1)==imask_inner );
        } else {
          neighbor_to_land = (  (*mask)(i,j+1)<maskfloating || (*mask)(i,j-1)<maskfloating || 
                                (*mask)(i+1,j)<maskfloating || (*mask)(i-1,j)<maskfloating || 
                                (*mask)(i+1,j+1)<maskfloating || (*mask)(i+1,j-1)<maskfloating ||
                                (*mask)(i-1,j+1)<maskfloating || (*mask)(i-1,j-1)<maskfloating );
        }


        if (neighbor_to_land ){
             // i.e. there is a grounded neighboring cell (which is not ice rise!)
            BOXMODELmask(i,j) = box_GL;
            lcounter_GLbox[shelf_id]++;

        } else if ((*mask)(i,j+1)==maskocean || (*mask)(i,j-1)==maskocean || 
                   (*mask)(i+1,j)==maskocean || (*mask)(i-1,j)==maskocean || 
                   (*mask)(i+1,j+1)==maskocean || (*mask)(i+1,j-1)==maskocean ||
                   (*mask)(i-1,j+1)==maskocean || (*mask)(i-1,j-1)==maskocean ){ // i.e. there is an ocean neighboring cell 
            BOXMODELmask(i,j) = box_IF;
            lcounter_CFbox[shelf_id]++;
        } else { // i.e., all other floating boxes
            BOXMODELmask(i,j) = box_unidentified;
            lcounter_box_unidentified++;
        }

      }else{ // i.e., not floating
        BOXMODELmask(i,j) = box_noshelf;
      }
    }
  }

  if (exicerises_set) { 
    ierr = ICERISESmask.end_access();   CHKERRQ(ierr);}
  ierr = mask->end_access();   CHKERRQ(ierr); 
  ierr = basins->end_access();   CHKERRQ(ierr);
  ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);

  //ierr = BOXMODELmask.beginGhostComm(); CHKERRQ(ierr);
  //ierr = BOXMODELmask.endGhostComm(); CHKERRQ(ierr);
  ierr = BOXMODELmask.update_ghosts(); CHKERRQ(ierr);
  

  ierr = PISMGlobalSum(&lcounter_box_unidentified, &counter_box_unidentified, grid.com); CHKERRQ(ierr);
  for(int k=0;k<numberOfBasins;k++) {
    ierr = PISMGlobalSum(&lcounter[k], &counter[k], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lcounter_CFbox[k], &counter_CFbox[k], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lcounter_GLbox[k], &counter_GLbox[k], grid.com); CHKERRQ(ierr);
   
    ierr = verbPrintf(5, grid.com,"  %d: cnt[k]=%.0f, cnt_CFbox=%.0f, cnt_GLbox=%.0f\n", k, counter[k], counter_CFbox[k], counter_GLbox[k]); CHKERRQ(ierr);
  }
 
  return 0;
}




//! Compute the boxmodelmask, calculate the extent of each box in each region
PetscErrorCode POoceanboxmodel::identifyBOXMODELmask() {

  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"A1c: in identify boxmodel mask rountine\n"); CHKERRQ(ierr);
  
  double lcounter_box_unidentified = counter_box_unidentified+1.0; 
  
  while((counter_box_unidentified > 0.0)&&(lcounter_box_unidentified != counter_box_unidentified)){

      lcounter_box_unidentified = counter_box_unidentified;
      ierr = verbPrintf(5, grid.com,"     cnt_box_unidentified=%.0f\n", lcounter_box_unidentified); CHKERRQ(ierr);

      for(int l=0;l<3;l++) { //FIXME size depends on how often this routine is called
        ierr = extendIFBox(); CHKERRQ(ierr);}
      ierr = extendGLBox(); CHKERRQ(ierr); 

  }

  // FIXME How to handle shelf-lakes (shelf which lies in the sheet), at the moment: GL_Box, better: Beckmann-Goose or no melting at all?
  if(counter_box_unidentified>0.0){ //if there are still floating cells which were not labels before, i.e. a shelf connected to an island (in the ex_icerises case)
    
    ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);
    ierr = basins->begin_access();   CHKERRQ(ierr);
    double lcounter_GLbox[numberOfBasins];
    double all_counter_GLbox[numberOfBasins];
    for (int k=0;k<numberOfBasins;k++){ lcounter_GLbox[k]=0.0; all_counter_GLbox[k]=0.0;}; 
    
    for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (BOXMODELmask(i,j)==box_unidentified ){ // i.e. this is an unidentified shelf cell with a neighbor that is in the IFbox
            ierr = verbPrintf(5, grid.com,"   left over i=%d, j=%d \n", i,j); CHKERRQ(ierr);
            BOXMODELmask(i,j)=box_GL;
            int shelf_id =(*basins)(i,j);
            lcounter_GLbox[shelf_id]++;
        }
      }
    }
    
    for(int k=0;k<numberOfBasins;k++) { 
      ierr = PISMGlobalSum(&lcounter_GLbox[k], &all_counter_GLbox[k], grid.com); CHKERRQ(ierr);
      counter_GLbox[k] += all_counter_GLbox[k];
    } 
    ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);
    ierr = basins->end_access();   CHKERRQ(ierr);
  }

  for(int k=0;k<numberOfBasins;k++){ ierr = verbPrintf(5, grid.com,"  %d: cnt[i] = %.0f, cnt_CFbox = %.0f, cnt_GLbox = %.0f, ratio_CF_box = %.3f, ratio_GL_box = %.3f\n", k, counter[k], counter_CFbox[k], counter_GLbox[k], counter_CFbox[k]/counter[k], counter_GLbox[k]/counter[k]); CHKERRQ(ierr);}
  return 0;
}





//! extent the grouningline box with neighboring unidentified shelf cells
PetscErrorCode POoceanboxmodel::extendGLBox() {
  PetscErrorCode ierr;

  double lcounter_box_unidentified = 0.0,
              all_counter_box_unidentified = 0.0;
  double lcounter_GLbox[numberOfBasins],
              all_counter_GLbox[numberOfBasins];

  for (int k=0;k<numberOfBasins;k++){ lcounter_GLbox[k]=0.0; all_counter_GLbox[k]=0.0; }; 

  ierr = mask->begin_access();   CHKERRQ(ierr);
  ierr = basins->begin_access();   CHKERRQ(ierr); 
  ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (BOXMODELmask(i,j)==box_unidentified && 
          (BOXMODELmask(i,j+1)==box_GL || BOXMODELmask(i,j-1)==box_GL || 
           BOXMODELmask(i+1,j)==box_GL || BOXMODELmask(i-1,j)==box_GL // || 
           //BOXMODELmask(i+1,j+1)==box_GL || BOXMODELmask(i+1,j-1)==box_GL || 
           //BOXMODELmask(i-1,j+1)==box_GL || BOXMODELmask(i-1,j-1)==box_GL
          )){ // i.e. this is an unidentified shelf cell with a neighbor that is in the GLbox
          BOXMODELmask(i,j)=box_neighboring; 
      }
    }
  }

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (BOXMODELmask(i,j)==box_neighboring){ 
          BOXMODELmask(i,j)=box_GL;
          lcounter_box_unidentified++; 
          int shelf_id=(*basins)(i,j);
          lcounter_GLbox[shelf_id]++;
      }
    }
  }


  ierr = mask->end_access();   CHKERRQ(ierr); 
  ierr = basins->end_access();   CHKERRQ(ierr);
  ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);
  //ierr = BOXMODELmask.beginGhostComm(); CHKERRQ(ierr);
  //ierr = BOXMODELmask.endGhostComm(); CHKERRQ(ierr);
  ierr = BOXMODELmask.update_ghosts(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&lcounter_box_unidentified, &all_counter_box_unidentified, grid.com); CHKERRQ(ierr);
  counter_box_unidentified -= all_counter_box_unidentified;

  for(int k=0;k<numberOfBasins;k++) { 
    ierr = PISMGlobalSum(&lcounter_GLbox[k], &all_counter_GLbox[k], grid.com); CHKERRQ(ierr);
    counter_GLbox[k] += all_counter_GLbox[k];
  } 

  return 0;
}



//! extend the ice_front box with neighboring unidentified shelf cells.
PetscErrorCode POoceanboxmodel::extendIFBox() {
  
  PetscErrorCode ierr;

  double lcounter_box_unidentified = 0.0; 
  double all_counter_box_unidentified = 0.0;
  double lcounter_CFbox[numberOfBasins];
  double all_counter_CFbox[numberOfBasins];
  for (int k=0;k<numberOfBasins;k++){ lcounter_CFbox[k] = 0.0; all_counter_CFbox[k] = 0.0;}

  ierr = mask->begin_access();   CHKERRQ(ierr); 
  ierr = basins->begin_access();   CHKERRQ(ierr);
  ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);
  
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (BOXMODELmask(i,j)==box_unidentified && 
          (BOXMODELmask(i,j+1)==box_IF || BOXMODELmask(i,j-1)==box_IF || 
           BOXMODELmask(i+1,j)==box_IF || BOXMODELmask(i-1,j)==box_IF  // || BOXMODELmask(i+1,j+1)==box_IF || BOXMODELmask(i+1,j-1)==box_IF || //FIXME eight neigbors?BOXMODELmask(i-1,j+1)==box_IF || BOXMODELmask(i-1,j-1)==box_IF
          )){ // i.e. this is an unidentified shelf cell with a neighbor that is in the IFbox
          BOXMODELmask(i,j)=box_neighboring;  
      }
    }
  }

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (BOXMODELmask(i,j)==box_neighboring){ 
          BOXMODELmask(i,j)=box_IF; 
          lcounter_box_unidentified++; 
          int shelf_id=(*basins)(i,j);
          lcounter_CFbox[shelf_id]++;
      }
    }
  }

  ierr = mask->end_access();   CHKERRQ(ierr); 
  ierr = basins->end_access();   CHKERRQ(ierr);
  ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);
  //ierr = BOXMODELmask.beginGhostComm(); CHKERRQ(ierr);
  //ierr = BOXMODELmask.endGhostComm(); CHKERRQ(ierr);
  ierr = BOXMODELmask.update_ghosts(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&lcounter_box_unidentified, &all_counter_box_unidentified, grid.com); CHKERRQ(ierr);
  counter_box_unidentified -= all_counter_box_unidentified;
  
  for(int k=0;k<numberOfBasins;k++) { 
    ierr = PISMGlobalSum(&lcounter_CFbox[k], &all_counter_CFbox[k], grid.com); CHKERRQ(ierr);
    counter_CFbox[k] += all_counter_CFbox[k];
  } 

  return 0;
}




/*!
Compute ocean temperature outside of the ice shelf cavities.
*/


PetscErrorCode POoceanboxmodel::oceanTemperature(const POBMConstants &cc) { 

  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"A2 : in ocean temp rountine\n"); CHKERRQ(ierr);

  ierr = mask->begin_access();   CHKERRQ(ierr);
  ierr = basins->begin_access();   CHKERRQ(ierr);
  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);

  ierr = Soc_base.begin_access();   CHKERRQ(ierr);
  ierr = Toc_base.begin_access();   CHKERRQ(ierr);
  ierr = Toc_anomaly.begin_access();   CHKERRQ(ierr);
  ierr = Toc.begin_access();   CHKERRQ(ierr);

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // make sure all temperatures are zero at the beginning of each timestep
      Toc(i,j) = 273.15; // in K
      Toc_base(i,j) = 273.15;  // in K
      Toc_anomaly(i,j) = 0.0;  // in K or °C
      Soc_base(i,j) = 0.0; // in psu


      if ((*mask)(i,j)==maskfloating){
        int shelf_id = (*basins)(i,j);
        Toc_base(i,j) = 273.15 + Toc_base_vec[shelf_id];
        Soc_base(i,j) =  Soc_base_vec[shelf_id];

        //! salinity and temperature for grounding line box
        if ( Soc_base(i,j) == 0.0 || Toc_base_vec[shelf_id] == 0.0 ) {
          ierr = verbPrintf(2, grid.com,"PISM_ERROR: Missing Soc_base and Toc_base for %d, %d, basin %d \n   Aborting... \n", i, j, shelf_id); CHKERRQ(ierr);
          PISMEnd();
        }


        // Add temperature anomalies from given nc-file  // FIXME different nc-files for each basin!
        if ((delta_T != NULL) && ocean_oceanboxmodel_deltaT_set) {
          //Toc_anomaly(i,j) = delta_T_factor * (*delta_T)(m_t + 0.5*m_dt);
          Toc_anomaly(i,j) = delta_T_factor * temp_anomaly;

        } else {

          Toc_anomaly(i,j) = 0.0;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //prevent ocean temp from being below pressure melting temperature
        //const double shelfbaseelev = - (cc.rhoi / cc.rhow) * (*ice_thickness)(i,j);

        const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4;
        // MUST be in dbar  // NOTE 1dbar = 10000 Pa = 1e4 kg m-1 s-2,
        // FIXME need to include atmospheric pressure?
        const double T_pmt = cc.a*Soc_base(i,j) + cc.b - cc.c*pressure;
        //ierr = verbPrintf(5, grid.com,"!!!!! T_pmt=%f, Ta=%f, Tb=%f, Toc=%f, Ta2=%f, Toc2=%f at %d,%d,%d\n", T_pmt,   Toc_anomaly(i,j),   Toc_base(i,j)-273.15,   Toc_base(i,j)-273.15 + Toc_anomaly(i,j),    PetscMax( T_pmt+273.15-Toc_base(i,j),   Toc_anomaly(i,j)),    Toc_base(i,j)-273.15+PetscMax( T_pmt+273.15-Toc_base(i,j),Toc_anomaly(i,j))  ,i  ,j,  shelf_id); CHKERRQ(ierr);

        Toc_anomaly(i,j) = PetscMax( T_pmt + 273.15 - Toc_base(i,j) , Toc_anomaly(i,j)); 
        /////////////////////////////////////////////////////////////////////////////////////////////////////

        Toc(i,j) = Toc_base(i,j) + Toc_anomaly(i,j); // in K


      } // end if herefloating
    } // end j
  } // end i

  ierr = mask->end_access();   CHKERRQ(ierr);
  ierr = basins->end_access();   CHKERRQ(ierr);
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);

  ierr = Soc_base.end_access();   CHKERRQ(ierr);
  ierr = Toc_base.end_access();   CHKERRQ(ierr);
  ierr = Toc_anomaly.end_access();   CHKERRQ(ierr);
  ierr = Toc.end_access();   CHKERRQ(ierr);

  return 0;
}



// NOTE Mean Gl_box meltrate is needed for basalMeltRateForIceFrontBox(). Here, mean is taken over all shelves in one drainage basin!

//! Compute the basal melt / refreezing rates for each shelf cell bordering the grounding line box
PetscErrorCode POoceanboxmodel::basalMeltRateForGroundingLineBox(const POBMConstants &cc) {
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"B1 : in basal melt rate gl rountine\n"); CHKERRQ(ierr);

  double lcounter_edge_of_GLbox_vector[numberOfBasins], 
              lmean_salinity_GLbox_vector[numberOfBasins], 
              lmean_meltrate_GLbox_vector[numberOfBasins], 
              lmean_overturning_GLbox_vector[numberOfBasins];

  for (int k=0;k<numberOfBasins;k++){ 
    lcounter_edge_of_GLbox_vector[k]=0.0;
    lmean_salinity_GLbox_vector[k]=0.0;
    lmean_meltrate_GLbox_vector[k]=0.0;
    lmean_overturning_GLbox_vector[k]=0.0;
  }


  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);
  ierr = basins->begin_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);
  ierr = T_star.begin_access(); CHKERRQ(ierr);
  ierr = Toc_base.begin_access(); CHKERRQ(ierr);
  ierr = Toc_anomaly.begin_access(); CHKERRQ(ierr);
  ierr = Toc_inCelsius.begin_access(); CHKERRQ(ierr);
  ierr = Toc.begin_access(); CHKERRQ(ierr);
  ierr = Soc_base.begin_access(); CHKERRQ(ierr);
  ierr = Soc.begin_access(); CHKERRQ(ierr);
  ierr = overturning.begin_access(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.begin_access(); CHKERRQ(ierr);

  double countHelpterm=0,
              lcountHelpterm=0;


  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {

      int shelf_id = (*basins)(i,j);

      // Make sure everything is at default values at the beginning of each timestep
      T_star(i,j) = 0.0; // in °C
      Toc_inCelsius(i,j) = 0.0; // in °C
      Soc(i,j) = 0.0; // in psu

      basalmeltrate_shelf(i,j) = 0.0;
      overturning(i,j) = 0.0;


      if (BOXMODELmask(i,j) == box_GL && shelf_id > 0.0){

        const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4; // MUST be in dbar  // NOTE 1dbar = 10000 Pa = 1e4 kg m-1 s-2
        // FIXME need to include atmospheric pressure?
        T_star(i,j) = cc.a*Soc_base(i,j) + cc.b - cc.c*pressure - (Toc_base(i,j) - 273.15 + Toc_anomaly(i,j)); // in °C

        double gamma_T_star,C1,g1;
        gamma_T_star = gamma_T_star_vec[shelf_id];
        C1 = C_vec[shelf_id];
        g1 = (counter_GLbox[shelf_id] * grid.dx * grid.dy) * gamma_T_star / (C1*cc.rho_star); 

        //! temperature for grounding line box
        double helpterm1 = g1/(cc.beta*(Soc_base(i,j) / (cc.nu*cc.lambda)) - cc.alpha);                  // in 1 / (1/°C) = °C
        double helpterm2 = (g1*T_star(i,j)) / (cc.beta*(Soc_base(i,j) / (cc.nu*cc.lambda)) - cc.alpha); // in °C / (1/°C) = °C^2

        if ((0.25*PetscSqr(helpterm1) -helpterm2) < 0.0) {
          //ierr = verbPrintf(5, grid.com,"PISM_ERROR: Tb=%f, Ta=%f, Tst=%f, Sb=%f  at %d, %d\n\n",Toc_base(i,j),Toc_anomaly(i,j),T_star(i,j),Soc_base(i,j),i,j); CHKERRQ(ierr);
          //ierr = verbPrintf(5, grid.com,"PISM_ERROR: h1=%f, h2=%f, h1sq=%f at %d, %d\n\n",helpterm1,helpterm2,0.25*PetscSqr(helpterm1),i,j); CHKERRQ(ierr);
          //ierr = verbPrintf(5, grid.com,"PISM_ERROR: square-root is negative! %f at %d, %d\n...with 0.25*helpterm^2=%f,helpterm2=%f,g1=%f,(cc.beta*(Soc_base(i,j)/(cc.nu*cc.lambda))-cc.alpha)=%f,Tstar=%f\n Not aborting, but setting sum to 0... \n", 0.25*PetscSqr(helpterm1) -helpterm2, i, j, 0.25*PetscSqr(helpterm1),helpterm2,g1,(cc.beta*(Soc_base(i,j) / (cc.nu*cc.lambda)) - cc.alpha),T_star(i,j)); CHKERRQ(ierr);
          //PISMEnd();
          helpterm2=0.25*PetscSqr(helpterm1);
          //FIXME: This might be wrong!
          lcountHelpterm+=1;
        }

        // NOTE Careful, Toc_base(i,j) is in K, Toc_inCelsius(i,j) NEEDS to be in °C!
        Toc_inCelsius(i,j) = (Toc_base(i,j)-273.15+Toc_anomaly(i,j)) - ( -0.5*helpterm1 + sqrt(0.25*PetscSqr(helpterm1) -helpterm2) );

        //! salinity for grounding line box
        Soc(i,j) = Soc_base(i,j) - (Soc_base(i,j) / (cc.nu*cc.lambda)) * ((Toc_base(i,j)-273.15+Toc_anomaly(i,j)) - Toc_inCelsius(i,j));  // in psu

        //! basal melt rate for grounding line box
        basalmeltrate_shelf(i,j) = (-gamma_T_star/(cc.nu*cc.lambda)) * (cc.a*Soc(i,j) + cc.b - cc.c*pressure - Toc_inCelsius(i,j));  // in m/s

        //! overturning
        // NOTE Actually, there is of course no overturning-FIELD, it is only a scalar for each shelf.
        // Here, I compute overturning as   MEAN[C1*cc.rho_star* (cc.beta*(Soc_base(i,j)-Soc(i,j)) - cc.alpha*((Toc_base(i,j)-273.15+Toc_anomaly(i,j))-Toc_inCelsius(i,j)))]
        // while in fact it should be   C1*cc.rho_star* (cc.beta*(Soc_base-MEAN[Soc(i,j)]) - cc.alpha*((Toc_base-273.15+Toc_anomaly)-MEAN[Toc_inCelsius(i,j)]))
        // which is the SAME since Soc_base, Toc_base and Toc_anomaly are the same FOR ALL i,j CONSIDERED, so this is just nomenclature!
        overturning(i,j) = C1*cc.rho_star* (cc.beta*(Soc_base(i,j)-Soc(i,j)) - cc.alpha*((Toc_base(i,j)-273.15+Toc_anomaly(i,j))-Toc_inCelsius(i,j))); // in m^3/s

        if (BOXMODELmask(i-1,j)==box_IF || BOXMODELmask(i+1,j)==box_IF || BOXMODELmask(i,j-1)==box_IF || BOXMODELmask(i,j+1)==box_IF){ 
        // i.e., if this cell is from the GL box and one of the neighbours is from the CF box - It is important to only take the border of the grounding line box 
        // to the calving front box into account, because the following mean value will be used to compute the value for the calving front box. I.e., this helps avoiding discontinuities!
          lcounter_edge_of_GLbox_vector[shelf_id]++;
          lmean_salinity_GLbox_vector[shelf_id] += Soc(i,j);
          lmean_meltrate_GLbox_vector[shelf_id] += basalmeltrate_shelf(i,j);
          lmean_overturning_GLbox_vector[shelf_id] += overturning(i,j);

          //ierr = verbPrintf(4, grid.com,"B1 : in basal melt rate gl rountine test1: %d,%d, %d: %.0f \n",i,j,shelf_id,lcounter_edge_of_GLbox_vector[shelf_id]); CHKERRQ(ierr);

        } // no else-case necessary since all variables are set to zero at the beginning of this routine

      }else { // i.e., not GL_box
          basalmeltrate_shelf(i,j) = 0.0;
      }
    } // end j
  } // end i

  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = basins->end_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);
  ierr = T_star.end_access(); CHKERRQ(ierr);
  ierr = Toc_base.end_access(); CHKERRQ(ierr);
  ierr = Toc_anomaly.end_access(); CHKERRQ(ierr);
  ierr = Toc_inCelsius.end_access(); CHKERRQ(ierr);
  ierr = Toc.end_access(); CHKERRQ(ierr);
  ierr = Soc_base.end_access(); CHKERRQ(ierr);
  ierr = Soc.end_access(); CHKERRQ(ierr);
  ierr = overturning.end_access(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.end_access(); CHKERRQ(ierr);


  for(int k=0;k<numberOfBasins;k++) {
    double counter_edge_of_GLbox_vector=0.0;
    ierr = PISMGlobalSum(&lcounter_edge_of_GLbox_vector[k], &counter_edge_of_GLbox_vector, grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lmean_meltrate_GLbox_vector[k], &mean_meltrate_GLbox_vector[k], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lmean_salinity_GLbox_vector[k], &mean_salinity_GLbox_vector[k], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lmean_overturning_GLbox_vector[k], &mean_overturning_GLbox_vector[k], grid.com); CHKERRQ(ierr);


    if (counter_edge_of_GLbox_vector>0.0){
      mean_salinity_GLbox_vector[k] = mean_salinity_GLbox_vector[k]/counter_edge_of_GLbox_vector;
      mean_meltrate_GLbox_vector[k] = mean_meltrate_GLbox_vector[k]/counter_edge_of_GLbox_vector;
      mean_overturning_GLbox_vector[k] = mean_overturning_GLbox_vector[k]/counter_edge_of_GLbox_vector;
    } else { // This means that there is no [cell from the GLbox neighboring a cell from the CFbox], NOT necessarily that there is no GLbox!
      mean_salinity_GLbox_vector[k]=0.0; mean_meltrate_GLbox_vector[k]=0.0; mean_overturning_GLbox_vector[k]=0.0;
    }

    ierr = verbPrintf(5, grid.com,"  %d: cnt=%.0f, sal=%.3f, melt=%.3e, over=%.1e \n", k,counter_edge_of_GLbox_vector,mean_salinity_GLbox_vector[k],mean_meltrate_GLbox_vector[k],mean_overturning_GLbox_vector[k]) ; CHKERRQ(ierr); 
  }

    ierr = PISMGlobalSum(&lcountHelpterm, &countHelpterm, grid.com); CHKERRQ(ierr);
    if (countHelpterm > 0) {
      ierr = verbPrintf(2, grid.com,"B1!: PISM_WARNING: square-root has been negative in %.0f cases!\n",countHelpterm); CHKERRQ(ierr);
    }

  return 0;
}




//! Compute the basal melt / refreezing rates for each shelf cell bordering the ice front box
PetscErrorCode POoceanboxmodel::basalMeltRateForIceFrontBox(const POBMConstants &cc) {
  PetscErrorCode ierr;  // FIXME redo all verbprintfs!
  ierr = verbPrintf(4, grid.com,"B2 : in bm ice front rountine\n"); CHKERRQ(ierr);

  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);
  ierr = basins->begin_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);

  ierr = T_star.begin_access(); CHKERRQ(ierr);
  ierr = Toc_base.begin_access(); CHKERRQ(ierr);
  ierr = Toc_anomaly.begin_access(); CHKERRQ(ierr);
  ierr = Toc_inCelsius.begin_access(); CHKERRQ(ierr);
  ierr = Toc.begin_access(); CHKERRQ(ierr);
  ierr = Soc_base.begin_access(); CHKERRQ(ierr);
  ierr = Soc.begin_access(); CHKERRQ(ierr);
  ierr = overturning.begin_access(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.begin_access(); CHKERRQ(ierr);

  double countk4=0,
              lcountk4=0,
              countGl0=0,
              lcountGl0=0,
              countSqr=0,
              lcountSqr=0,
              countMean0=0,
              lcountMean0=0;


  //! The ice front box = BOX I
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {  // FIXME REPAIR
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {

      int shelf_id = (*basins)(i,j);

      if (BOXMODELmask(i,j)==box_IF && shelf_id > 0.0){

        const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4; // MUST be in dbar  // NOTE 1dbar = 10000 Pa = 1e4 kg m-1 s-2
        // FIXME need to include atmospheric pressure?
        T_star(i,j) = cc.a*Soc_base(i,j) + cc.b - cc.c*pressure - (Toc_base(i,j) - 273.15 + Toc_anomaly(i,j)); // in °C

        double  gamma_T_star,area_GLbox,area_CFbox,mean_salinity_in_GLbox,mean_meltrate_in_GLbox,mean_overturning_in_GLbox;

        gamma_T_star = gamma_T_star_vec[shelf_id]; 
        area_CFbox = (counter_CFbox[shelf_id] * grid.dx * grid.dy); 
        area_GLbox = (counter_GLbox[shelf_id] * grid.dx * grid.dy);  
        mean_salinity_in_GLbox = mean_salinity_GLbox_vector[shelf_id];
        mean_meltrate_in_GLbox = mean_meltrate_GLbox_vector[shelf_id];
        mean_overturning_in_GLbox = mean_overturning_GLbox_vector[shelf_id];


        if (area_GLbox==0) { // if there is no grounding line box in the current basin, set BOXMODELmask to 3 and compute basal melt rate by Beckmann-Gose
            //ierr = verbPrintf(5, grid.com,"!!!! GLBOX =0 , basin = %d at %d,%d, \n   ", shelf_id,i,j); CHKERRQ(ierr);
            BOXMODELmask(i,j) = box_other;
            lcountGl0+=1;

        } else {

            double k1,k2,k3,k4,k5,k6;
            k1 = (area_CFbox*gamma_T_star)/(cc.nu*cc.lambda);
            // in (m^2*m/s)/(°C) = m^3/(s*°C)
            k2 = 1/(mean_overturning_in_GLbox + area_CFbox*gamma_T_star);
            // in s/m^3
            k3 = (area_CFbox*gamma_T_star*T_star(i,j) - cc.nu*cc.lambda*area_GLbox*mean_meltrate_in_GLbox);
            // in m^2 * m/s * °C = m^3 * °C / s
            k4 = (-k1*k2*area_CFbox*gamma_T_star*cc.a + k1*cc.a);
            // in m^3/(s*°C) * s/m^3 * m^2 * m/s * °C/psu = m^3/(s*psu)
            k5 = (mean_overturning_in_GLbox + Soc_base(i,j)*k1*k2*area_CFbox*gamma_T_star*cc.a  - k1*Soc_base(i,j)*cc.a - k1*T_star(i,j) + k1*k2*k3);
            // in m^3/s
            k6 = (k1*Soc_base(i,j)*T_star(i,j) - k1*k2*Soc_base(i,j)*k3  - area_GLbox*mean_meltrate_in_GLbox*mean_salinity_in_GLbox);
            // in psu * m^3/s

            //ierr = verbPrintf(2, grid.com,"!!!! ACF=%.3e, AGL=%.3e, MS=%.3f, MM=%.3f, MO=%.3f, basin = %d at %d,%d, \n   ",  area_CFbox,area_GLbox,mean_salinity_in_GLbox,mean_meltrate_in_GLbox,mean_overturning_in_GLbox,shelf_id,i,j); CHKERRQ(ierr);
            //ierr = verbPrintf(2, grid.com,"!!!! k1=%.3f, k2=%.3f, k3=%.3f, k4=%.3f, k5=%.3f, k6=%.3f, basin = %d at %d,%d, \n   ",k1,k2,k3,k4,k5,k6,shelf_id,i,j); CHKERRQ(ierr);

            //! salinity for calving front box
            if (k4 == 0.0) {
              //ierr = verbPrintf(5, grid.com,"PISM_ERROR: Division by zero! k4=%f at %d, %d\n   Aborting... \n", k4, i, j); CHKERRQ(ierr);
              //ierr = verbPrintf(5, grid.com,"PISM_ERROR: Probably mean_overturning_in_GLbox = %f is zero, check if there is a grounding line box in basin %d , \n   ", mean_overturning_in_GLbox, shelf_id); CHKERRQ(ierr);
              // FIXME rewrite this warning? Do not stop but calculate melt rates according to Beckmann-Gose?          
              //PISMEnd();
              lcountk4+=1;
              BOXMODELmask(i,j) = box_other;
              continue;
            } 

            if ((0.25*PetscSqr(k5/k4) - (k6/k4)) < 0.0) {
              //ierr = verbPrintf(5, grid.com,"PISM_ERROR: Square-root is negative! %f at %d, %d\n...with 0.25*PetscSqr((k5/k4))=%f,(k6/k4)=%f,k5=%f,k6=%f,k4=%f\n   Aborting... \n", 0.25*PetscSqr(k5/k4) - (k6/k4), i, j,0.25*PetscSqr((k5/k4)), (k6/k4),k5,k6,k4) ; CHKERRQ(ierr);
              //PISMEnd();
              lcountSqr+=1;
              BOXMODELmask(i,j) = box_other;
              continue;
            }

            Soc(i,j) = Soc_base(i,j) - ( -0.5*(k5/k4) - sqrt(0.25*PetscSqr(k5/k4) - (k6/k4)) ); // in psu

            //! temperature for calving front box
            // NOTE Careful, Toc_base(i,j) is in K, Toc_inCelsius(i,j) NEEDS to be in °C!
            Toc_inCelsius(i,j) = (Toc_base(i,j)-273.15+Toc_anomaly(i,j)) - ( k2*area_CFbox*gamma_T_star*cc.a*(Soc_base(i,j)-Soc(i,j)) - k2*k3 );

            //! basal melt rate for calving front box
            basalmeltrate_shelf(i,j) = (-gamma_T_star/(cc.nu*cc.lambda)) * (cc.a*Soc(i,j) + cc.b - cc.c*pressure - Toc_inCelsius(i,j)); // in m/s

            if (mean_salinity_in_GLbox == 0.0 || mean_meltrate_in_GLbox == 0.0 || mean_overturning_in_GLbox == 0.0){
              // this must not occur since there must always be a GL_box neighbor 
              //ierr = verbPrintf(5, grid.com, "PISM_ERROR: DETECTION CFBOX: There is no neighbouring grounding line box for this calving front box at %d,%d! \nThis will lead to a zero k4 and in turn to NaN in Soc, Toc_inCelsius and basalmeltrate_shelf. After the next massContExplicitStep(), H will be NaN, too! This will cause ks in temperatureStep() to be NaN and lead to a Segmentation Violation! \nIn particular: basin_id=%d, BOXMODELmask=%f, H=%f, T_star=%f, \narea_GLbox=%e, area_CFbox=%e, mean_salinity_in_GLbox=%f, mean_meltrate_in_GLbox=%e, mean_overturning_in_GLbox=%e, \nk1=%e,k2=%e,k3=%e,k4=%e,k5=%e,k6=%e, \nToc_base=%f, Toc_anomaly=%f, Toc_inCelsius=%f, Toc=%f, Soc_base=%f, Soc=%f, basalmeltrate_shelf=%e \n   Aborting... \n", i,j, shelf_id, BOXMODELmask(i,j), (*ice_thickness)(i,j), T_star(i,j), area_GLbox,area_CFbox,mean_salinity_in_GLbox,mean_meltrate_in_GLbox,mean_overturning_in_GLbox,k1,k2,k3,k4,k5,k6, Toc_base(i,j), Toc_anomaly(i,j), Toc_inCelsius(i,j), Toc(i,j), Soc_base(i,j), Soc(i,j), basalmeltrate_shelf(i,j)); CHKERRQ(ierr);
              //PISMEnd();
              lcountMean0+=1;
              BOXMODELmask(i,j) = box_other;
              continue;
            
            }
        }
      } // NOTE NO else-case, since  basalMeltRateForGroundingLineBox() and basalMeltRateForOtherShelves() cover all other cases and we would overwrite those results here.
    } // end j
  } // end i

  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = basins->end_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);

  ierr = T_star.end_access(); CHKERRQ(ierr);
  ierr = Toc_base.end_access(); CHKERRQ(ierr);
  ierr = Toc_anomaly.end_access(); CHKERRQ(ierr);
  ierr = Toc_inCelsius.end_access(); CHKERRQ(ierr);
  ierr = Toc.end_access(); CHKERRQ(ierr);
  ierr = Soc_base.end_access(); CHKERRQ(ierr);
  ierr = Soc.end_access(); CHKERRQ(ierr);
  ierr = overturning.end_access(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&lcountk4, &countk4, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lcountGl0, &countGl0, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lcountSqr, &countSqr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lcountMean0, &countMean0, grid.com); CHKERRQ(ierr);
    
  if (countk4 > 0) {
    ierr = verbPrintf(2, grid.com,"B2!: PISM_WARNING: k4 is zero in %.0f case(s)!\n",countk4); CHKERRQ(ierr);
  }
  if (countGl0 > 0) {
    ierr = verbPrintf(2, grid.com,"B2!: PISM_WARNING: no grounding line box in basin in %.0f case(s)!\n",countGl0); CHKERRQ(ierr);
  }
  if (countSqr > 0) {
    ierr = verbPrintf(2, grid.com,"B2!: PISM_WARNING: square root is negative in %.0f case(s)!\n",countSqr); CHKERRQ(ierr);
  }
  if (countMean0 > 0) {
    ierr = verbPrintf(2, grid.com,"B2!: PISM_WARNING: mean of salinity, meltrate or overturning is zero in %.0f case(s)!\n",countMean0); CHKERRQ(ierr);
  }

  return 0;
}




//! Convert Toc_inCelsius from °C to K and write into Toc for the .nc-file; NOTE It is crucial, that Toc_inCelsius is in °C for the computation of the basal melt rate
//! Compute the melt rate for all other ice shelves.
PetscErrorCode POoceanboxmodel::basalMeltRateForOtherShelves(const POBMConstants &cc) {
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"B3 : in bm others rountine\n"); CHKERRQ(ierr);


  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);
  ierr = basins->begin_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);
  ierr = Toc_base.begin_access(); CHKERRQ(ierr);
  ierr = Toc_anomaly.begin_access(); CHKERRQ(ierr);
  ierr = Toc_inCelsius.begin_access(); CHKERRQ(ierr);
  ierr = Toc.begin_access(); CHKERRQ(ierr);
  ierr = overturning.begin_access(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.begin_access(); CHKERRQ(ierr); // NOTE meltrate has units:   J m-2 s-1 / (J kg-1 * kg m-3) = m s-1
  ierr = heatflux.begin_access(); CHKERRQ(ierr);

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {

      int shelf_id = (*basins)(i,j);

      if (shelf_id == 0) { // boundary of computational domain

        basalmeltrate_shelf(i,j) = 0.0;

      } else if (BOXMODELmask(i,j)==box_other ) {

        Toc(i,j) = Toc_base(i,j) + Toc_anomaly(i,j); // in K, NOTE: Toc_base is already in K, so no (+273.15)
        // default: compute the melt rate from the temperature field according to beckmann_goosse03 (see below)
        
        const double shelfbaseelev = - (cc.rhoi / cc.rhow) * (*ice_thickness)(i,j);

        //FIXME: for consistency reasons there should be constants a,b,c, gamma_T used
        double T_f = 273.15 + (cc.a*cc.meltSalinity + cc.b2 + cc.c*shelfbaseelev); // add 273.15 to get it in Kelvin... 35 is the salinity

        heatflux(i,j) = cc.meltFactor * cc.rhow * cc.c_p_ocean * cc.gamma_T_o * (Toc(i,j) - T_f);  // in W/m^2
        basalmeltrate_shelf(i,j) = heatflux(i,j) / (cc.latentHeat * cc.rhoi); // in m s-1

      } else if (shelf_id > 0.0) {

        Toc(i,j) = 273.15 + Toc_inCelsius(i,j) + Toc_anomaly(i,j); // in K  // FIXME I think Toc should not occur in any of the routines before!
      
      } else { // This must not happen

        ierr = verbPrintf(2, grid.com,"PISM_ERROR: [rank %d] at %d, %d  -- basins(i,j)=%d causes problems.\n   Aborting... \n",grid.rank, i, j, shelf_id); CHKERRQ(ierr);
        PISMEnd();
      }
    } // end j
  } // end i



  ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = basins->end_access(); CHKERRQ(ierr);
  ierr = Toc_base.end_access(); CHKERRQ(ierr);
  ierr = Toc_anomaly.end_access(); CHKERRQ(ierr);
  ierr = Toc_inCelsius.end_access(); CHKERRQ(ierr);
  ierr = Toc.end_access(); CHKERRQ(ierr);
  ierr = overturning.end_access(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.end_access(); CHKERRQ(ierr);
  ierr = heatflux.end_access(); CHKERRQ(ierr);

  return 0;
}



PetscErrorCode POoceanboxmodel::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr = shelfbtemp.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode POoceanboxmodel::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr = shelfbmassflux.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode POoceanboxmodel::sea_level_elevation(double &result) {
  result = sea_level;
  return 0;
}

PetscErrorCode POoceanboxmodel::melange_back_pressure_fraction(IceModelVec2S &result) {
  PetscErrorCode ierr = result.set(0.0); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode POoceanboxmodel::write_variables(std::set<std::string> vars, const PIO& nc) {
  PetscErrorCode ierr;

  ierr = PGivenClimate<POModifier,PISMOceanModel>::write_variables(vars, nc); CHKERRQ(ierr);

  if (set_contains(vars, "shelfbtemp")) {
    ierr = shelfbtemp.write(nc); CHKERRQ(ierr);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    ierr = shelfbmassflux.write(nc); CHKERRQ(ierr);
  }

  
  if (set_contains(vars, "BOXMODELmask")) {         
    ierr = BOXMODELmask.write(nc); CHKERRQ(ierr);  
  }

  if (set_contains(vars, "OCEANMEANmask")) {        
    ierr = OCEANMEANmask.write(nc); CHKERRQ(ierr);  
  }

  if (exicerises_set) {
    if (set_contains(vars, "ICERISESmask")) {       
      ierr = ICERISESmask.write(nc); CHKERRQ(ierr);  
    }
  }

  if (set_contains(vars, "Soc")) {                  ierr = Soc.write(nc); CHKERRQ(ierr); }
  if (set_contains(vars, "Soc_base")) {             ierr = Soc_base.write(nc); CHKERRQ(ierr); }
  if (set_contains(vars, "Toc")) {                  ierr = Toc.write(nc); CHKERRQ(ierr); }
  if (set_contains(vars, "Toc_base")) {             ierr = Toc_base.write(nc); CHKERRQ(ierr);  }
  if (set_contains(vars, "Toc_inCelsius")) {        ierr = Toc_inCelsius.write(nc); CHKERRQ(ierr);  }
  if (set_contains(vars, "T_star")) {               ierr = T_star.write(nc); CHKERRQ(ierr);  }
  if (set_contains(vars, "Toc_anomaly")) {          ierr = Toc_anomaly.write(nc); CHKERRQ(ierr); }
  if (set_contains(vars, "overturning")) {          ierr = overturning.write(nc); CHKERRQ(ierr); }
  if (set_contains(vars, "heatflux")) {             ierr = heatflux.write(nc); CHKERRQ(ierr);  }
  if (set_contains(vars, "basalmeltrate_shelf")) {  ierr = basalmeltrate_shelf.write(nc); CHKERRQ(ierr); }

  return 0;
}

