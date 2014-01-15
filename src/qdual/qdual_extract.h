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
}

#endif
