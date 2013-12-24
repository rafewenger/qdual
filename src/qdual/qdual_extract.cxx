/// \file qdual_extract.cxx
/// Subroutines for extracting dual isosurface mesh

/*
  QDUAL: Quality Dual Isosurface Generation
  Copyright (C) 2014 Arindam Bhattacharya, Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "ijkgrid_macros.h"
#include "ijkdualtable.h"
#include "ijkisopoly.txx"
#include "ijktime.txx"

#include "qdual_extract.h"

using namespace IJK;
using namespace IJKDUALTABLE;


// **************************************************
// EXTRACT ROUTINES
// **************************************************

/// Extract isosurface polytopes
/// Returns list representing isosurface polytopes
/// @param scalar_grid = scalar grid data
/// @param isovalue = isosurface scalar value
/// @param iso_poly[] = vector of isosurface polygope vertices
///   iso_simplices[numv_per_poly*ip+k] = 
///     cube containing k'th vertex of polytope ip.
void QDUAL::extract_dual_isopoly
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_poly,
 DUALISO_INFO & dualiso_info)
{
  dualiso_info.time.extract = 0;

  clock_t t0 = clock();

  // initialize output
  iso_poly.clear();

  if (scalar_grid.NumCubeVertices() < 1) { return; }

  IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, scalar_grid, VERTEX_INDEX) {

    extract_dual_isopoly_around_bipolar_edge
      (scalar_grid, isovalue, iend0, edge_dir, iso_poly);
  }

  clock_t t1 = clock();
  clock2seconds(t1-t0, dualiso_info.time.extract);
}

/// Extract isosurface polytopes
/// Returns list representing isosurface polytopes
/// @param scalar_grid = scalar grid data
/// @param isovalue = isosurface scalar value
/// @param iso_poly[] = vector of isosurface polygope vertices
///   iso_simplices[numv_per_poly*ip+k] = 
///     cube containing k'th vertex of polytope ip.
/// @param facet_vertex[i] = Edge of cube containing iso_poly[i].
void QDUAL::extract_dual_isopoly
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_poly,
 std::vector<FACET_VERTEX_INDEX> & facet_vertex,
 DUALISO_INFO & dualiso_info)
{
  dualiso_info.time.extract = 0;

  clock_t t0 = clock();

  // initialize output
  iso_poly.clear();
  facet_vertex.clear();

  if (scalar_grid.NumCubeVertices() < 1) { return; }

  IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, scalar_grid, VERTEX_INDEX) {

    extract_dual_isopoly_around_bipolar_edge
      (scalar_grid, isovalue, iend0, edge_dir, iso_poly, facet_vertex);
  }

  clock_t t1 = clock();
  clock2seconds(t1-t0, dualiso_info.time.extract);
}

