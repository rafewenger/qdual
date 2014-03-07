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
/// @param iso_poly[] = vector of isosurface polytope vertices
///   iso_poly[numv_per_poly*ip+k] = 
///      cube containing k'th vertex of polytope ip.
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


/// Extract isosurface quadrilaterals.
/// Returns list representing isosurface quadrilaterals.
/// Returns also list of directions orthogonal to quadrilaterals.
/// @param scalar_grid = scalar grid data
/// @param isovalue = isosurface scalar value
/// @param[out] iso_quad[] = vector of isosurface polygope vertices
///   iso_simplices[numv_per_poly*ip+k] = 
///     cube containing k'th vertex of polytope ip.
/// @param[out] facet_vertex[i] = Edge of cube containing iso_quad[i].
/// @param[out] orth_dir[i] = Direction orthgonal to quadrilateral i.
void QDUAL::extract_dual_isoquad
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_quad,
 std::vector<FACET_VERTEX_INDEX> & facet_vertex,
 std::vector<DIRECTION_TYPE> & orth_dir,
 DUALISO_INFO & dualiso_info)
{
  dualiso_info.time.extract = 0;

  clock_t t0 = clock();

  // initialize output
  iso_quad.clear();
  facet_vertex.clear();
  orth_dir.clear();

  if (scalar_grid.NumCubeVertices() < 1) { return; }

  IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, scalar_grid, VERTEX_INDEX) {

    std::vector<ISO_VERTEX_INDEX>::size_type prev_size = iso_quad.size();

    extract_dual_isoquad_around_bipolar_edge
      (scalar_grid, isovalue, iend0, edge_dir, iso_quad, facet_vertex);

    if (prev_size != iso_quad.size()) {
      // edge is bipolar.
      orth_dir.push_back(edge_dir);
    }

  }

  clock_t t1 = clock();
  clock2seconds(t1-t0, dualiso_info.time.extract);
}


// Specialized version of extract.
// First quad vertex is always lowest/leftmost.
void QDUAL::extract_dual_isoquad_around_bipolar_edge
(const DUALISO_SCALAR_GRID_BASE & scalar_grid, 
 const SCALAR_TYPE isovalue, 
 const VERTEX_INDEX iend0, const DIRECTION_TYPE edge_dir,
 std::vector<ISO_VERTEX_INDEX> & isoquad,
 std::vector<FACET_VERTEX_INDEX> & facet_vertex)
{
  VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);

  bool is_end0_positive = true;
  if (scalar_grid.Scalar(iend0) < isovalue) 
    { is_end0_positive = false; };

  bool is_end1_positive = true;
  if (scalar_grid.Scalar(iend1) < isovalue) 
    { is_end1_positive = false; };

  if (!is_end0_positive && is_end1_positive) {
    IJK::extract_dual_isopoly_around_edge
      (scalar_grid, iend0, iend1, edge_dir, isoquad, facet_vertex);
  }
  else if (is_end0_positive && !is_end1_positive) {
    std::vector<ISO_VERTEX_INDEX>::size_type k = isoquad.size();
    IJK::extract_dual_isopoly_around_edge
      (scalar_grid, iend0, iend1, edge_dir, isoquad, facet_vertex);
    std::swap(isoquad[k+1], isoquad[k+2]);
    std::swap(facet_vertex[k+1], facet_vertex[k+2]);
  }
}
