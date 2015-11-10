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
#include <memory>

#include "ijk.txx"
#include "qdual_types.h"
#include "qdual_datastruct.h"

/// qdual classes and routines.
namespace QDUAL {

// **************************************************
// DUAL CONTOURING 
// **************************************************

  /// Dual Contouring Algorithm. (main function red1)
  void quality_dual_contouring
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

  /// Extract isosurface using Dual Contouring algorithm
  /// Returns list of isosurface simplex vertices
  ///  and list of isosurface vertex coordinates
  void dual_contouring
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, 
   const DUALISO_DATA_FLAGS & dualiso_data_flags,
   std::vector<VERTEX_INDEX> & quad_vert,
   std::vector<VERTEX_INDEX> & dual_edge,
   std::vector<COORD_TYPE> & vertex_coord,
   std::vector<DUAL_ISOVERT> &iso_vlist,
   IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   MERGE_DATA & merge_data, 
   DUALISO_INFO & dualiso_info);

  // **************************************************
  // CONVERT QUADRILATERALS TO TRIANGLES
  // **************************************************

  /// Convert quadrilaterals (embedded in 3D) to triangles.
  /// Assumes input isosurface is quadrilaterals embedded in 3D.
  void convert_quad_to_tri
  (const DUALISO_DATA & dualiso_data, DUAL_ISOSURFACE & dual_isosurface);
}

#endif
