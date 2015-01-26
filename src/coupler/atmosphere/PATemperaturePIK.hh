// Copyright (C)  2008-2015 PISM Authors

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

#ifndef __PATemperaturePIK_hh
#define __PATemperaturePIK_hh

#include "PAYearlyCycle.hh"
#include "Timeseries.hh"


class PATemperaturePIK : public PAYearlyCycle {
public:
  PATemperaturePIK(IceGrid &g, const PISMConfig &conf);
   // : PAYearlyCycle(g, conf)

  virtual ~PATemperaturePIK();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(double my_t, double my_dt);
  //virtual PetscErrorCode precip_time_series(int i, int j, double *values);
  virtual PetscErrorCode mean_precipitation(IceModelVec2S &result);
protected:
  bool precipitation_correction, precip_increase_per_degree_set, temp_huybrechts_dewolde99_set, temp_era_interim_set, temp_era_interim_sin_set;
  PetscReal precip_increase_per_degree;
  Timeseries *delta_T;
  IceModelVec2S *lat, *surfelev;
};


#endif	// __PATemperaturePIK_hh
