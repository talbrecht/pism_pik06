// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#ifndef _PAYEARLYCYCLE_H_
#define _PAYEARLYCYCLE_H_

#include "PISMAtmosphere.hh"
#include "iceModelVec.hh"

//! A class containing an incomplete implementation of an atmosphere model
//! based on a temperature parameterization using mean annual and mean July
//! (mean summer) temperatures and a cosine yearly cycle. Uses a stored
//! (constant in time) precipitation field.
class PAYearlyCycle : public PISMAtmosphereModel {
public:
  PAYearlyCycle(IceGrid &g, const PISMConfig &conf);
  virtual ~PAYearlyCycle();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual void add_vars_to_output(std::string keyword, std::set<std::string> &result);
  virtual PetscErrorCode define_variables(std::set<std::string> vars, const PIO &nc, PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(std::set<std::string> vars, const PIO &nc);
  //! This method implements the parameterization.
  virtual PetscErrorCode update(double my_t, double my_dt) = 0;
  virtual PetscErrorCode mean_precipitation(IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result);
  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result);

  virtual PetscErrorCode init_timeseries(double *ts, unsigned int N);
  virtual PetscErrorCode temp_time_series(int i, int j, double *values);
  virtual PetscErrorCode precip_time_series(int i, int j, double *values);
protected:
  PISMVars *variables;
  double snow_temp_july_day;
  std::string reference, precip_filename;
  IceModelVec2S air_temp_mean_annual, air_temp_mean_july, precipitation, precip_standard;
  NCSpatialVariable air_temp_snapshot;
  std::vector<double> m_ts_times, m_cosine_cycle;
private:
  PetscErrorCode allocate_PAYearlyCycle();
};

#endif /* _PAYEARLYCYCLE_H_ */
