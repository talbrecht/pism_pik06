// Copyright (C) 2008-2015 PISM Authors
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

// Implementation of the atmosphere model using constant-in-time precipitation
// and a cosine yearly cycle for near-surface air temperatures.

// This includes the SeaRISE Greenland parameterization.

#include "PATemperaturePIK.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include "PISMTime.hh"
#include <assert.h>
#include "PISMConfig.hh"

///// PATemperaturePIK
PATemperaturePIK::PATemperaturePIK(IceGrid &g, const PISMConfig &conf)
  : PAYearlyCycle(g, conf)
{
    precipitation_correction = false;
    delta_T = NULL;
}

PATemperaturePIK::~PATemperaturePIK() 
{
    delete delta_T;
}

PetscErrorCode PATemperaturePIK::init(PISMVars &vars) {
  PetscErrorCode ierr;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  ierr = verbPrintf(2, grid.com,
     "* Initializing the atmosphere model PATemperaturePIK.\n"); CHKERRQ(ierr);

  ierr = PAYearlyCycle::init(vars); CHKERRQ(ierr);

  // initialize pointers to fields the parameterization depends on:
  surfelev = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (!surfelev) SETERRQ(grid.com, 1, "ERROR: surface_altitude is not available");

  lat = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
  if (!lat) SETERRQ(grid.com, 1, "ERROR: latitude is not available");

/// Surface (annual mean and summer mean) temperature parametrization: 
  ierr = PISMOptionsIsSet("-temp_huybrechts_dewolde99", temp_huybrechts_dewolde99_set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-temp_era_interim", temp_era_interim_set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-temp_era_interim_sin", temp_era_interim_sin_set); CHKERRQ(ierr);

  if (temp_huybrechts_dewolde99_set) {
    ierr = verbPrintf(2, grid.com,
                      "    Near-surface air temperature is parameterized as in Huybrechts & De Wolde (1999).\n"); CHKERRQ(ierr);
  }else if (temp_era_interim_set) {
    ierr = verbPrintf(2, grid.com,
                      "    Near-surface air temperature is parameterized based on ERA INTERIM data.\n"); CHKERRQ(ierr);
  }else if (temp_era_interim_sin_set) {
    ierr = verbPrintf(2, grid.com,
                      "    Near-surface air temperature is parameterized based on ERA INTERIM data with a sin(lat) dependence.\n"); CHKERRQ(ierr);
  }else{
    ierr = verbPrintf(2, grid.com,
                      "    Near-surface annual mean air temperature is parameterized as in Martin et al. (2011),\n"
          "    and near-surface summer mean air temperature is computed as anomaly to the Huybrechts & De Wolde (1999) - temperature.\n"); CHKERRQ(ierr);
  }


  /// Precipitation parametrization: 
  ierr = PISMOptionsIsSet("-precip_change", precipitation_correction); CHKERRQ(ierr);

  if (precipitation_correction) {
    bool delta_T_set;
    std::string delta_T_file;

    ierr = PISMOptionsString("-precip_change",
                             "Specifies the air temperature offsets file to use with -precip_change",
                             delta_T_file, delta_T_set); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com, 
                      "    Reading delta_T data from forcing file %s for -precip_change actions ...\n",
                      delta_T_file.c_str());  CHKERRQ(ierr);

    delta_T = new Timeseries(&grid, "delta_T", config.get_string("time_dimension_name"));
    delta_T->set_units("Kelvin", "Kelvin");
    delta_T->set_dimension_units(grid.time->units_string(), "");
    delta_T->set_attr("long_name", "air temperature offsets");

    //ierr = delta_T->read(delta_T_file, grid.time->use_reference_date()); CHKERRQ(ierr);
    PIO nc(grid.com, "netcdf3", grid.get_unit_system());
    ierr = nc.open(delta_T_file, PISM_NOWRITE); CHKERRQ(ierr);
    {
      ierr = delta_T->read(nc, grid.time); CHKERRQ(ierr);
    }
    ierr = nc.close(); CHKERRQ(ierr);

  }


  ierr = PISMOptionsReal("-precip_increase_per_degree",
                           "Precipitation increase per degree of warming ",
                           precip_increase_per_degree, precip_increase_per_degree_set); CHKERRQ(ierr);

  if (precip_increase_per_degree_set){
    if (precipitation_correction == false){
      PetscPrintf(grid.com,
		  "PISM ERROR: Options -precip_increase_per_degree requires that option -precip_change <filename> is set,\n"
		  "but -precip_change is either not set or the file specification is missing.\n Aborting\n"); CHKERRQ(ierr);
      PISMEnd();
    }else{
      PetscReal precip_percentage = (precip_increase_per_degree - 1.0) * 100; // in percent per degree of warming
      ierr = verbPrintf(2, grid.com,
                    "      Precipitation is increased by %1.1f percent per degree of warming\n", precip_percentage); CHKERRQ(ierr);
    }
  }

  return 0;
}


PetscErrorCode PATemperaturePIK::precip_time_series(int i, int j, double *values) {

  for (unsigned int k = 0; k < m_ts_times.size(); k++)
    values[k] = precipitation(i,j);

  return 0;
}

// Scale present-day precipitation field
PetscErrorCode PATemperaturePIK::mean_precipitation(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = PAYearlyCycle::mean_precipitation(result); CHKERRQ(ierr);

  if ((delta_T != NULL) && precipitation_correction) {
//    // ... as in Pollard & De Conto (2012), Eqn (34b): 
//    ierr = result.scale(pow (2.0, (0.1* (*delta_T)(t + 0.5 * dt)))); CHKERRQ(ierr); // scale by 2^(0.1*DeltaT)

//    // ... using a factor based on Clausius-Clapeyron:  
    if (precip_increase_per_degree_set){
      ierr = result.scale(pow (precip_increase_per_degree, ((*delta_T)(m_t + 0.5 * m_dt)))); CHKERRQ(ierr); // scale by factor^(DeltaT), i.e., ((factor - 1.0)+ 100) PERCENT precip increase per degree
    }else{
      ierr = result.scale(pow (1.05, ((*delta_T)(m_t + 0.5 * m_dt)))); CHKERRQ(ierr); // scale by 1.05^(DeltaT), i.e., 5% precip increase per degree, DEFAULT CASE
    }
  }

  return 0;
}

//! \brief Updates mean annual and mean July near-surface air temperatures.
//! Note that the precipitation rate is time-independent and does not need
//! to be updated.
PetscErrorCode PATemperaturePIK::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;

if (lat->metadata().has_attribute("missing_at_bootstrap")) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: latitude variable was missing at bootstrap;\n"
                       "  SeaRISE-Greenland atmosphere model depends on latitude and would return nonsense!!\n");
    CHKERRQ(ierr);
    PISMEnd();
  }


  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12))
    return 0;

  m_t  = my_t;
  m_dt = my_dt;

  IceModelVec2S &h = *surfelev, &lat_degN = *lat;

  ierr = h.begin_access();   CHKERRQ(ierr);
  ierr = lat_degN.begin_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.begin_access();  CHKERRQ(ierr);
  ierr = air_temp_mean_july.begin_access();  CHKERRQ(ierr);// NOTE: Of course, this is not the mean July, but the mean SUMMER temperature! (Here only denoted as air_temp_mean_july so that the field can be interpreted correctly in the PDD-routines (which were originally meant for the Northern Hemisphere).)

  for (int i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j<grid.ys+grid.ym; ++j) {

	if (temp_huybrechts_dewolde99_set){
	  PetscReal gamma_a;
	  if (h(i,j) < 1500.0) {
	    gamma_a = -0.005102;
	  }else{
	    gamma_a = -0.014285;
	  }
	  air_temp_mean_annual(i,j) = 273.15 + 34.46 + gamma_a * h(i,j) - 0.68775 * lat_degN(i,j)*(-1.0); // = TMA, mean annual temperature in Huybrechts & DeWolde (1999)
	  air_temp_mean_july(i,j) = 273.15 + 14.81 - 0.00692 * h(i,j) - 0.27937 * lat_degN(i,j)*(-1.0); // = TMS, mean summer temperature in Huybrechts & DeWolde (1999)  

	}else if (temp_era_interim_set){  // parametrization based on multiple regression analysis of ERA INTERIM data
	  air_temp_mean_annual(i,j) = 273.15 + 29.2 - 0.0082 * h(i,j) - 0.576 * lat_degN(i,j)*(-1.0);
	  air_temp_mean_july(i,j)   = 273.15 + 16.5 - 0.0068 * h(i,j) - 0.248 * lat_degN(i,j)*(-1.0);

	}else if (temp_era_interim_sin_set){  // parametrization based on multiple regression analysis of ERA INTERIM data with sin(lat)
	  air_temp_mean_annual(i,j) = 273.15 - 2.0 -0.0082*h(i,j) + 18.4 * (sin(3.1415*lat_degN(i,j)/180)+0.8910)/(1-0.8910);
	  air_temp_mean_july(i,j)   = 273.15 + 3.2 -0.0067*h(i,j) +  8.3 * (sin(3.1415*lat_degN(i,j)/180)+0.8910)/(1-0.8910);

	}else{
	  // annual mean temperature = Martin et al. (2011) parametrization
	  // summer mean temperature = anomaly to Huybrechts & DeWolde (1999)
	  air_temp_mean_annual(i,j) = 273.15 + 30 - 0.0075 * h(i,j) - 0.68775 * lat_degN(i,j)*(-1.0);  // surface temperature parameterization as in Martin et al. 2011, Eqn. 2.0.2

	  PetscReal gamma_a;
	      if (h(i,j) < 1500.0) {
		gamma_a = -0.005102;
	      }else{
		gamma_a = -0.014285;
	      }

	  PetscReal TMA = 273.15 + 34.46 + gamma_a * h(i,j) - 0.68775 * lat_degN(i,j)*(-1.0); // = TMA, mean annual temperature in Huybrechts & DeWolde (1999)
	  PetscReal TMS = 273.15 + 14.81 - 0.00692 * h(i,j) - 0.27937 * lat_degN(i,j)*(-1.0); // = TMS, mean summer temperature in Huybrechts & DeWolde (1999)

	  air_temp_mean_july(i,j) = air_temp_mean_annual(i,j) + (TMS - TMA);   

	  //// ALTERNATIVE:
	  //// annual mean temperature = Martin et al. (2011)
	  //// summer mean temperature = Huybrechts & DeWolde (1999)
	  //      air_temp_mean_annual(i,j) = 273.15 + 30 - 0.0075 * h(i,j) - 0.68775 * lat_degN(i,j)*(-1.0);  // annual mean temperature as in Martin et al. 2011, Eqn. 2.0.2
	  //      air_temp_mean_july(i,j) = 273.15 + 14.81 - 0.00692 * h(i,j) - 0.27937 * lat_degN(i,j)*(-1.0); // = TMS, mean summer temperature in Huybrechts & DeWolde (1999) 

	} 

    }
  }

  ierr = h.end_access();   CHKERRQ(ierr);
  ierr = lat_degN.end_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.end_access();  CHKERRQ(ierr);
  ierr = air_temp_mean_july.end_access();  CHKERRQ(ierr);

  return 0;
}
