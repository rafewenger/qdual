/// \file ijkdual_position.cxx
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

#include "ijkcoord.txx"
#include "ijkinterpolate.txx"
#include "ijkisopoly.txx"

#include "ijkdualtable.h"

#include "qdual_position.h"

using namespace IJK;
using namespace IJKDUALTABLE;


// **************************************************
// SINGLE ISOSURFACE VERTEX IN A GRID CUBE
// **************************************************

/// Position dual isosurface vertices in cube centers
void QDUAL::position_dual_isovertices_cube_center
(const DUALISO_GRID & grid,
 const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
{
  const int dimension = grid.Dimension();

  for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
    VERTEX_INDEX iv = vlist[i];

    grid.ComputeCubeCenterCoord(iv, coord+i*dimension);
  }
}

/// Position dual isosurface vertices in cube centers.
/// C++ STL vector format for array coord[].
void QDUAL::position_dual_isovertices_cube_center
(const DUALISO_GRID & grid,
 const std::vector<ISO_VERTEX_INDEX> & vlist, 
 std::vector<COORD_TYPE> & coord)
{
  const int dimension = grid.Dimension();

  coord.resize(vlist.size()*dimension);
  position_dual_isovertices_cube_center(grid, vlist, &(coord.front()));
}

///   of isosurface-edge intersections
void QDUAL::position_dual_isovertices_centroid
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
{
  const int dimension = scalar_grid.Dimension();
  GRID_COORD_TYPE grid_coord[dimension];
  COORD_TYPE vcoord[dimension];
  COORD_TYPE coord0[dimension];
  COORD_TYPE coord1[dimension];
  COORD_TYPE coord2[dimension];

  for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
    VERTEX_INDEX iv = vlist[i];

    int num_intersected_edges = 0;
    set_coord(dimension, 0.0, vcoord);

    for (int edge_dir = 0; edge_dir < dimension; edge_dir++)
      for (int k = 0; k < scalar_grid.NumFacetVertices(); k++) {
        VERTEX_INDEX iend0 = scalar_grid.FacetVertex(iv, edge_dir, k);
        VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);

        SCALAR_TYPE s0 = scalar_grid.Scalar(iend0);
        bool is_end0_positive = true;
        if (s0 < isovalue)
          { is_end0_positive = false; };

        SCALAR_TYPE s1 = scalar_grid.Scalar(iend1);
        bool is_end1_positive = true;
        if (s1 < isovalue)
          { is_end1_positive = false; };

        if (is_end0_positive != is_end1_positive) {

          scalar_grid.ComputeCoord(iend0, coord0);
          scalar_grid.ComputeCoord(iend1, coord1);

          linear_interpolate_coord
            (dimension, s0, coord0, s1, coord1, isovalue, coord2);

          add_coord(dimension, vcoord, coord2, vcoord);

          num_intersected_edges++;
        }

      }

    if (num_intersected_edges > 0) {
      multiply_coord
        (dimension, 1.0/num_intersected_edges, vcoord, coord+i*dimension);
    }
    else {
      scalar_grid.ComputeCubeCenterCoord(iv, coord+i*dimension);
    }

  }
}


/// Position dual isosurface vertices in centroid 
///   of isosurface-edge intersections
void QDUAL::position_dual_isovertices_centroid
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist, 
 std::vector<COORD_TYPE> & coord)

{
  const int dimension = scalar_grid.Dimension();

  coord.resize(vlist.size()*dimension);
  position_dual_isovertices_centroid
    (scalar_grid, isovalue, vlist, &(coord.front()));
}


// **************************************************
// ALLOW MULTIPLE ISOSURFACE VERTICES IN A GRID CUBE
// **************************************************

/// Position dual isosurface vertices using centroids.
void QDUAL::position_dual_isovertices_centroid_multi
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 COORD_TYPE * coord)
{
  const int dimension = scalar_grid.Dimension();
  const int num_cube_vertices = scalar_grid.NumCubeVertices();
  COORD_TYPE vcoord[dimension];
  COORD_TYPE coord0[dimension];
  COORD_TYPE coord1[dimension];
  COORD_TYPE coord2[dimension];
  IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);

  for (ISO_VERTEX_INDEX i = 0; i < iso_vlist.size(); i++) {
    VERTEX_INDEX icube = iso_vlist[i].cube_index;
    FACET_VERTEX_INDEX ipatch = iso_vlist[i].patch_index;
    TABLE_INDEX it = iso_vlist[i].table_index;

    int num_intersected_edges = 0;
    IJK::set_coord(dimension, 0.0, vcoord);

    for (int ie = 0; ie < cube.NumEdges(); ie++) {
      if (isodual_table.IsBipolar(it, ie)) {
        if (isodual_table.IncidentIsoVertex(it, ie) == ipatch) {
          int k0 = cube.EdgeEndpoint(ie, 0);
          int k1 = cube.EdgeEndpoint(ie, 1);
          VERTEX_INDEX iend0 = scalar_grid.CubeVertex(icube, k0);
          VERTEX_INDEX iend1 = scalar_grid.CubeVertex(icube, k1);
          SCALAR_TYPE s0 = scalar_grid.Scalar(iend0);
          SCALAR_TYPE s1 = scalar_grid.Scalar(iend1);

          scalar_grid.ComputeCoord(iend0, coord0);
          scalar_grid.ComputeCoord(iend1, coord1);

          linear_interpolate_coord
            (dimension, s0, coord0, s1, coord1, isovalue, coord2);

          add_coord(dimension, vcoord, coord2, vcoord);

          num_intersected_edges++;
        }
      }
    }

    if (num_intersected_edges > 0) {
      multiply_coord
        (dimension, 1.0/num_intersected_edges, vcoord, coord+i*dimension);
    }
    else {
      scalar_grid.ComputeCubeCenterCoord(icube, coord+i*dimension);
    }
  }
}

/// Position dual isosurface vertices using centroids.
void QDUAL::position_dual_isovertices_centroid_multi
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const std::vector<DUAL_ISOVERT> & iso_vlist, 
 std::vector<COORD_TYPE> & coord)
{
  const int dimension = scalar_grid.Dimension();

  coord.resize(iso_vlist.size()*dimension);
  position_dual_isovertices_centroid_multi
    (scalar_grid, isodual_table,isovalue, iso_vlist,
     &(coord.front()));
}

/// Position dual isosurface vertices near cube centers.
/// More than one vertex can be in a cube.
/// If cube contains multiple isosurface then vertices are positioned
///   near but not on cube center.
void QDUAL::position_dual_isovertices_near_cube_center_multi
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 const COORD_TYPE offset,
 COORD_TYPE * coord)
{
  const int dimension = scalar_grid.Dimension();
  COORD_TYPE vcoord[dimension];
  IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);
  IJK::UNIT_CUBE<int,int,int> unit_cube(dimension);
  bool intersects_facet[2*dimension];

  for (VERTEX_INDEX i = 0; i < iso_vlist.size(); i++) {
    VERTEX_INDEX icube = iso_vlist[i].cube_index;
    FACET_VERTEX_INDEX ipatch = iso_vlist[i].patch_index;
    TABLE_INDEX it = iso_vlist[i].table_index;

    if (isodual_table.NumIsoVertices(it) == 1) {
      scalar_grid.ComputeCubeCenterCoord(icube, coord+i*dimension);
    }
    else {

      scalar_grid.ComputeCubeCenterCoord(icube, vcoord);

      // Set intersects facet to false.
      for (int d = 0; d < dimension; d++) {
        intersects_facet[d] = false;
        intersects_facet[2*d] = false;
      }

      for (int ie = 0; ie < cube.NumEdges(); ie++) {
        if (isodual_table.IsBipolar(it, ie)) {
          if (isodual_table.IncidentIsoVertex(it, ie) == ipatch) {
            int k0 = cube.EdgeEndpoint(ie, 0);
            int k1 = cube.EdgeEndpoint(ie, 1);

            for (int d = 0; d < dimension; d++) {
              int c = unit_cube.VertexCoord(k0, d);
              if (c == unit_cube.VertexCoord(k1, d)) {
                if (c > 0) 
                  { intersects_facet[2*d] = true; }
                else
                  { intersects_facet[d] = true; }
              }
            }
          }
        }
      }

      for (int d = 0; d < dimension; d++) {

        if (intersects_facet[d] != intersects_facet[2*d]) {
          if (intersects_facet[d]) 
            { vcoord[d] -= offset; }
          else
            { vcoord[d] += offset; }
        }
      }
      IJK::copy_coord(dimension, vcoord, coord+i*dimension);
    }
  }
}

/// Position dual isosurface vertices near cube centers.
/// More than one vertex can be in a cube.
/// C++ STL vector format for array coord[].
void QDUAL::position_dual_isovertices_near_cube_center_multi
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 const COORD_TYPE offset,
 std::vector<COORD_TYPE> & coord)
{
  const int dimension = scalar_grid.Dimension();

  coord.resize(iso_vlist.size()*dimension);
  position_dual_isovertices_near_cube_center_multi
    (scalar_grid, isodual_table, isovalue, iso_vlist, offset, &(coord.front()));
}
