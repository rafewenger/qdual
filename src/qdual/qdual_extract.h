/// \file qdual_extract.h
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

#ifndef _QDUAL_EXTRACT_
#define _QDUAL_EXTRACT_

#include "ijkcube.txx"

#include "ijkdualtable.h"

#include "qdual_types.h"
#include "qdual_datastruct.h"

namespace QDUAL {

  /// Extract isosurface polytopes
  /// Returns list representing isosurface polytopes
  /// @param scalar_grid = scalar grid data
  /// @param isovalue = isosurface scalar value
  /// @param iso_poly[] = vector of isosurface polytope vertices
  ///   iso_poly[numv_per_poly*ip+k] = 
  ///      cube containing k'th vertex of polytope ip.
  void extract_dual_isopoly
    (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
     const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_poly,
     DUALISO_INFO & dualiso_info);

  /// Extract isosurface polytopes
  /// Returns list representing isosurface polytopes
  /// @param scalar_grid = scalar grid data
  /// @param isovalue = isosurface scalar value
  /// @param iso_poly[] = vector of isosurface polytope vertices
  ///   iso_simplices[numv_per_poly*ip+k] = 
  ///     cube containing k'th vertex of polytope ip.
  /// @param cube_edge[i] = Edge of cube containing iso_poly[i].
  void extract_dual_isopoly
    (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
     const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_poly,
     std::vector<FACET_VERTEX_INDEX> & facet_vertex,
     DUALISO_INFO & dualiso_info);

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
  void extract_dual_isoquad
    (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
     const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_quad,
     std::vector<FACET_VERTEX_INDEX> & facet_vertex,
     std::vector<DIRECTION_TYPE> & orth_dir,
     DUALISO_INFO & dualiso_info);

  // Specialized version of extract.
  // First quad vertex is always lowest/leftmost.
  void extract_dual_isoquad_around_bipolar_edge
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iend0, const DIRECTION_TYPE edge_dir,
   std::vector<ISO_VERTEX_INDEX> & isoquad,
   std::vector<FACET_VERTEX_INDEX> & facet_vertex);

}

#endif
