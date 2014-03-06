/// \file qdual_position.h
/// Position dual isosurface vertices

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


#ifndef _QDUAL_POSITION_
#define _QDUAL_POSITION_

#include <vector>

#include "ijkdualtable.h"

#include "qdual_types.h"
#include "qdual_datastruct.h"
#include "qdual_table.h"


namespace QDUAL {

// **************************************************
// SINGLE ISOSURFACE VERTEX IN A GRID CUBE
// **************************************************

  /// Position dual isosurface vertices in cube centers.
  void position_dual_isovertices_cube_center
    (const DUALISO_GRID & grid,
     const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord);

  /// Position dual isosurface vertices in cube centers.
  /// C++ STL vector format for array coord[].
  void position_dual_isovertices_cube_center
    (const DUALISO_GRID & grid, 
     const std::vector<ISO_VERTEX_INDEX> & vlist, 
     std::vector<COORD_TYPE> & coord);

  /// Position dual isosurface vertices in centroid 
  ///   of isosurface-edge intersections.
  void position_dual_isovertices_centroid
    (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
     const SCALAR_TYPE isovalue,
     const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord);

  /// Position dual isosurface vertices in centroid 
  ///   of isosurface-edge intersections.
  /// C++ STL vector format for array coord[].
  void position_dual_isovertices_centroid
    (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
     const SCALAR_TYPE isovalue,
     const std::vector<ISO_VERTEX_INDEX> & vlist, 
     std::vector<COORD_TYPE> & coord);


// **************************************************
// ALLOW MULTIPLE ISOSURFACE VERTICES IN A GRID CUBE
// **************************************************

  /// Position dual isosurface vertices using centroids.
  void position_dual_isovertices_centroid_multi
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   COORD_TYPE * coord);

  /// Position dual isosurface vertices using centroids.
  /// More than one vertex can be in a cube.
  /// Version with coord[] defined as an std::vector.
  void position_dual_isovertices_centroid_multi
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<DUAL_ISOVERT> & iso_vlist, 
   std::vector<COORD_TYPE> & coord);

  /// Position dual isosurface vertices near cube centers.
  /// More than one vertex can be in a cube.
  /// If cube contains multiple isosurface then vertices are positioned
  ///   near but not on cube center.
  void position_dual_isovertices_near_cube_center_multi
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   const COORD_TYPE offset,
   COORD_TYPE * coord);

  /// Position dual isosurface vertices near cube centers.
  /// More than one vertex can be in a cube.
  /// C++ STL vector format for array coord[].
  void position_dual_isovertices_near_cube_center_multi
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   const COORD_TYPE offset,
   std::vector<COORD_TYPE> & coord);

  /// Position dual isosurface vertices at random location in cube.
  /// @param seed Random generator seed.
  void position_dual_isovertices_random
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   const RANDOM_SEED_TYPE seed,
   COORD_TYPE * coord);

  /// Position dual isosurface vertices at random location in cube.
  /// @param seed Random generator seed.
  /// C++ STL vector format for array coord[].
  void position_dual_isovertices_random
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   const RANDOM_SEED_TYPE seed,
   std::vector<COORD_TYPE> & coord);

  /// Position dual isosurface vertices at random location in cube.
  /// Don't position near restricted facets.
  /// @param seed Seed for random number generator.
  /// @param flag_V1w_close If true, isosurface vertices with degree dimension
  ///                       are placed 1/3 from the corner grid vertex.
  void position_dual_isovertices_random
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const QDUAL_TABLE & qdual_table,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   const RANDOM_SEED_TYPE seed,
   const bool flag_V1w_close,
   COORD_TYPE * coord);

  /// Position dual isosurface vertices at random location in cube.
  /// Don't position near restricted facets.
  /// C++ STL vector format for array coord[].
  void position_dual_isovertices_random
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const QDUAL_TABLE & qdual_table,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   const RANDOM_SEED_TYPE seed,
   const bool flag_V1w_close,
   std::vector<COORD_TYPE> & coord);

};

#endif
