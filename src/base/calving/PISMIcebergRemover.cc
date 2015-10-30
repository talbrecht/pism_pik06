/* Copyright (C) 2013, 2014, 2015 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "PISMIcebergRemover.hh"
#include "connected_components.hh"
#include "Mask.hh"
#include "PISMVars.hh"

PISMIcebergRemover::PISMIcebergRemover(IceGrid &g, const PISMConfig &conf)
  : PISMComponent(g, conf) {

  PetscErrorCode ierr = allocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISM ERROR: failed to allocate PISMIcebergRemover.\n");
    PISMEnd();
  }
}

PISMIcebergRemover::~PISMIcebergRemover() {
  PetscErrorCode ierr = deallocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISM ERROR: failed to deallocate PISMIcebergRemover.\n");
    PISMEnd();
  }
}

PetscErrorCode PISMIcebergRemover::init(PISMVars &vars) {
  m_bcflag = dynamic_cast<IceModelVec2Int*>(vars.get("bcflag"));
  return 0;
}

/**
 * Use PISM's ice cover mask to update ice thickness, removing "icebergs".
 *
 * @param[in,out] pism_mask PISM's ice cover mask
 * @param[in,out] ice_thickness ice thickness
 */
PetscErrorCode PISMIcebergRemover::update(IceModelVec2Int &pism_mask,
                                          IceModelVec2S &ice_thickness) {
  PetscErrorCode ierr;
  const int
    mask_grounded_ice = 1,
    mask_floating_ice = 2;
  MaskQuery M(pism_mask);

  // prepare the mask that will be handed to the connected component
  // labeling code:
  {
    ierr = m_iceberg_mask.begin_access(); CHKERRQ(ierr);
    ierr = pism_mask.begin_access(); CHKERRQ(ierr);
    for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
      for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
        if (M.grounded_ice(i,j) == true)
          m_iceberg_mask(i,j) = mask_grounded_ice;
        else if (M.floating_ice(i,j) == true)
          m_iceberg_mask(i,j) = mask_floating_ice;
      }
    }

    // Mark icy SSA Dirichlet B.C. cells as "grounded" because we
    // don't want them removed.
    if (m_bcflag) {
      ierr = m_bcflag->begin_access(); CHKERRQ(ierr);
      for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
        for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
          if (m_bcflag->as_int(i,j) == 1 && M.icy(i,j))
            m_iceberg_mask(i,j) = mask_grounded_ice;
        }
      }
      ierr = m_bcflag->end_access(); CHKERRQ(ierr);
    }
    ierr = pism_mask.end_access(); CHKERRQ(ierr);
    ierr = m_iceberg_mask.end_access(); CHKERRQ(ierr);
  }

  // identify icebergs using serial code on processor 0:
  {
    ierr = m_iceberg_mask.put_on_proc0(m_mask_p0); CHKERRQ(ierr);

    if (grid.rank == 0) {
      double *mask;
      ierr = VecGetArray(m_mask_p0, &mask); CHKERRQ(ierr);

      cc(mask, grid.Mx, grid.My, true, mask_grounded_ice);

      ierr = VecRestoreArray(m_mask_p0, &mask); CHKERRQ(ierr);
    }

    ierr = m_iceberg_mask.get_from_proc0(m_mask_p0); CHKERRQ(ierr);
  }

  // correct ice thickness and the cell type mask using the resulting
  // "iceberg" mask:
  {
    ierr = m_iceberg_mask.begin_access(); CHKERRQ(ierr);
    ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
    ierr = pism_mask.begin_access(); CHKERRQ(ierr);

    if (m_bcflag != NULL) {
      // if SSA Dirichlet B.C. are in use, do not modify mask and ice
      // thickness at Dirichlet B.C. locations
      ierr = m_bcflag->begin_access(); CHKERRQ(ierr);
      for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
        for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
          if (m_iceberg_mask(i,j) > 0.5 && (*m_bcflag)(i,j) < 0.5) {
            ice_thickness(i,j) = 0.0;
            pism_mask(i,j)     = MASK_ICE_FREE_OCEAN;
          }
        }
      }
      ierr = m_bcflag->end_access(); CHKERRQ(ierr);
    } else {
      for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
        for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
          if (m_iceberg_mask(i,j) > 0.5) {
            ice_thickness(i,j) = 0.0;
            pism_mask(i,j)     = MASK_ICE_FREE_OCEAN;
          }
        }
      }
    }
    ierr = pism_mask.end_access(); CHKERRQ(ierr);
    ierr = ice_thickness.end_access(); CHKERRQ(ierr);
    ierr = m_iceberg_mask.end_access(); CHKERRQ(ierr);
  }

  // update ghosts of the mask and the ice thickness (then surface
  // elevation can be updated redundantly)
  ierr = pism_mask.update_ghosts(); CHKERRQ(ierr);
  ierr = ice_thickness.update_ghosts(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMIcebergRemover::allocate() {
  PetscErrorCode ierr;

  ierr = m_iceberg_mask.create(grid, "iceberg_mask", WITH_GHOSTS); CHKERRQ(ierr);
  ierr = m_iceberg_mask.allocate_proc0_copy(m_mask_p0); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMIcebergRemover::deallocate() {
  PetscErrorCode ierr;

  ierr = VecDestroy(&m_mask_p0); CHKERRQ(ierr);

  return 0;
}

void PISMIcebergRemover::add_vars_to_output(std::string, std::set<std::string> &) {
  // empty
}

PetscErrorCode PISMIcebergRemover::define_variables(std::set<std::string>, const PIO &,
                                                    PISM_IO_Type) {
  // empty
  return 0;
}

PetscErrorCode PISMIcebergRemover::write_variables(std::set<std::string>, const PIO& ) {
  // empty
  return 0;
}
