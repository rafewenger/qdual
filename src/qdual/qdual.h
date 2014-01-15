/// \file qdual.h
/// Quality dual contouring isosurface generation

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

/*!
  \mainpage QDual: Quality dual isosurface generation.

  QDUAL is a program for generating isosurfaces with quality mesh elements
  using the dual contouring algorithm.  Quadrilateral and triangle angles 
  are bounded away from 0 and 180.
*/

#ifndef _QDUAL_
#define _QDUAL_

#include <string>

#include "ijk.txx"
#include "qdual_types.h"
#include "qdual_datastruct.h"

/// ijkdual classes and routines.
namespace QDUAL {

// **************************************************
// DUAL CONTOURING 
// **************************************************

  /// Dual Contouring Algorithm.
  void dual_contouring
    (const DUALISO_DATA & dualiso_data, const SCALAR_TYPE isovalue, 
     DUAL_ISOSURFACE & dual_isosurface, DUALISO_INFO & dualiso_info);

// **************************************************
// DUAL CONTOURING
// **************************************************

  /// Dual Contouring Algorithm.
  void dual_contouring
   (const DUALISO_SCALAR_GRID_BASE & scalar_grid, 
    const SCALAR_TYPE isovalue, 
    std::vector<VERTEX_INDEX> & quad_vert, 
    std::vector<COORD_TYPE> & vertex_coord);

  /// Dual Contouring Algorithm.
  /// Represents each grid edge by a single integer.
  /// @param flag_separate_neg  If true, isosurface patches
  ///        separate negative grid vertices.
  /// @param[out] quad_vert[] Array of quadrilateral vertices.
  ///        quad_vert[iq*4+j] is index to j'th vertex of quadrilateral iq.
  /// @param[out] vertex_coord[] Array of vertex coordinates.
  ///        vertex_coord[i*dimension+j] is j'th coordinate of vertex i.
  /// @param merge_data Data structure for merging edges.  
  ///        Requires memory of size(MERGE_INDEX) for each grid edge.
  void dual_contouring
   (const DUALISO_SCALAR_GRID_BASE & scalar_grid, 
    const SCALAR_TYPE isovalue, 
    const VERTEX_POSITION_METHOD vertex_position_method,
    const bool flag_select_split,
    const bool flag_separate_neg,
    std::vector<VERTEX_INDEX> & quad_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    MERGE_DATA & merge_data, DUALISO_INFO & dualiso_info);

  // **************************************************
  // CONVERT QUADRILATERALS TO TRIANGLES
  // **************************************************

  /// Convert quadrilaterals (embedded in 3D) to triangles.
  /// Assumes input isosurface is quadrilaterals embedded in 3D.
  void convert_quad_to_tri
  (const DUALISO_DATA & dualiso_data, DUAL_ISOSURFACE & dual_isosurface);
}

#endif
