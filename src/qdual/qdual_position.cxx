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
  ARRAY<GRID_COORD_TYPE> grid_coord(dimension);
  ARRAY<COORD_TYPE> vcoord(dimension);
  ARRAY<COORD_TYPE> coord0(dimension);
  ARRAY<COORD_TYPE> coord1(dimension);
  ARRAY<COORD_TYPE> coord2(dimension);

  for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
    VERTEX_INDEX iv = vlist[i];

    int num_intersected_edges = 0;
    set_coord(dimension, 0.0, vcoord.Ptr());

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

          scalar_grid.ComputeCoord(iend0, coord0.Ptr());
          scalar_grid.ComputeCoord(iend1, coord1.Ptr());

          linear_interpolate_coord
            (dimension, s0, coord0.Ptr(), s1, coord1.Ptr(), isovalue, coord2.Ptr());

          add_coord(dimension, vcoord.Ptr(), coord2.Ptr(), vcoord.Ptr());

          num_intersected_edges++;
        }

      }

    if (num_intersected_edges > 0) {
      multiply_coord
        (dimension, 1.0/num_intersected_edges, vcoord.Ptr(), coord+i*dimension);
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
  ARRAY<COORD_TYPE> vcoord(dimension);
  ARRAY<COORD_TYPE> coord0(dimension);
  ARRAY<COORD_TYPE> coord1(dimension);
  ARRAY<COORD_TYPE> coord2(dimension);

  IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);

  for (ISO_VERTEX_INDEX i = 0; i < iso_vlist.size(); i++) {
    VERTEX_INDEX icube = iso_vlist[i].cube_index;
    FACET_VERTEX_INDEX ipatch = iso_vlist[i].patch_index;
    TABLE_INDEX it = iso_vlist[i].table_index;

    int num_intersected_edges = 0;
    IJK::set_coord(dimension, 0.0, vcoord.Ptr());

    for (int ie = 0; ie < cube.NumEdges(); ie++) {
      if (isodual_table.IsBipolar(it, ie)) {
        if (isodual_table.IncidentIsoVertex(it, ie) == ipatch) {
          int k0 = cube.EdgeEndpoint(ie, 0);
          int k1 = cube.EdgeEndpoint(ie, 1);
          VERTEX_INDEX iend0 = scalar_grid.CubeVertex(icube, k0);
          VERTEX_INDEX iend1 = scalar_grid.CubeVertex(icube, k1);
          SCALAR_TYPE s0 = scalar_grid.Scalar(iend0);
          SCALAR_TYPE s1 = scalar_grid.Scalar(iend1);

          scalar_grid.ComputeCoord(iend0, coord0.Ptr());
          scalar_grid.ComputeCoord(iend1, coord1.Ptr());

          linear_interpolate_coord
            (dimension, s0, coord0.Ptr(), s1, coord1.Ptr(), isovalue, coord2.Ptr());

          add_coord(dimension, vcoord.Ptr(), coord2.Ptr(), vcoord.Ptr());

          num_intersected_edges++;
        }
      }
    }

    if (num_intersected_edges > 0) {
      multiply_coord
        (dimension, 1.0/num_intersected_edges, vcoord.Ptr(), coord+i*dimension);
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
  ARRAY<COORD_TYPE> vcoord(dimension);
  IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);
  IJK::UNIT_CUBE<int,int,int> unit_cube(dimension);
  ARRAY<bool> intersects_facet(2*dimension);
 

  for (VERTEX_INDEX i = 0; i < iso_vlist.size(); i++) {
    VERTEX_INDEX icube = iso_vlist[i].cube_index;
    FACET_VERTEX_INDEX ipatch = iso_vlist[i].patch_index;
    TABLE_INDEX it = iso_vlist[i].table_index;

    if (isodual_table.NumIsoVertices(it) == 1) {
      scalar_grid.ComputeCubeCenterCoord(icube, coord+i*dimension);
    }
    else {

      scalar_grid.ComputeCubeCenterCoord(icube, vcoord.Ptr());

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
      IJK::copy_coord(dimension, vcoord.Ptr(), coord+i*dimension);
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


// **************************************************
// GENERATE RANDOM POSITION (FOR TESTING)
// **************************************************

/// Position dual isosurface vertices at random location in cube.
void QDUAL::position_dual_isovertices_random
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 const RANDOM_SEED_TYPE seed,
 COORD_TYPE * coord)
{
  const int dimension = scalar_grid.Dimension();
  ARRAY<COORD_TYPE> cube_coord(dimension);
  const int NUM_INTERVALS(1000);

  srand(seed);

  for (ISO_VERTEX_INDEX i = 0; i < iso_vlist.size(); i++) {
    VERTEX_INDEX icube = iso_vlist[i].cube_index;

    scalar_grid.ComputeCoord(icube, cube_coord.Ptr());

    for (int d = 0; d < dimension; d++) {
      COORD_TYPE x = rand()%(NUM_INTERVALS+1);
      x = x/NUM_INTERVALS;
      coord[i*dimension+d] = cube_coord[d] + x;
    }
  }
}

/// Position dual isosurface vertices at random location in cube.
/// C++ STL vector format for array coord[].
void QDUAL::position_dual_isovertices_random
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 const RANDOM_SEED_TYPE seed,
 std::vector<COORD_TYPE> & coord)
{
  const int dimension = scalar_grid.Dimension();

  coord.resize(iso_vlist.size()*dimension);
  position_dual_isovertices_random
    (scalar_grid, iso_vlist, seed, &(coord.front()));
}

/// Position dual isosurface vertices at random location in cube.
/// Don't position near restricted facets.
void QDUAL::position_dual_isovertices_random
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const QDUAL_TABLE & qdual_table,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 const RANDOM_SEED_TYPE seed,
 const bool flag_V1w_close,
 const int num_intervals,
 COORD_TYPE * coord)
{
  const int dimension = scalar_grid.Dimension();
  ARRAY<COORD_TYPE> cube_coord(dimension);
  COORD_TYPE xmin, xmax;
  QDUAL_TABLE::DIR_BITS mask;

  srand(seed);

  for (ISO_VERTEX_INDEX i = 0; i < iso_vlist.size(); i++) {
    VERTEX_INDEX icube = iso_vlist[i].cube_index;
    FACET_VERTEX_INDEX ipatch = iso_vlist[i].patch_index;
    TABLE_INDEX it = iso_vlist[i].table_index;
    QDUAL_TABLE::DIR_BITS connect_dir = qdual_table.connectDir(it, ipatch);

    scalar_grid.ComputeCoord(icube, cube_coord.Ptr());

    int num_connect = 0;
    for (int k = 0; k < 2*dimension; k++) {
      mask = (1 << k);
      if ((mask & connect_dir) != 0) { num_connect++; }
    }

    for (int d = 0; d < dimension; d++) {
      COORD_TYPE x = rand()%(num_intervals+1);

      // Vertices are not on the boundary of the cubes.
      xmin = 0.0001;
      xmax = 0.9999;

      if (num_connect == dimension) {
        mask = (1 << d);
        if ((mask & connect_dir) != 0) { xmax = 1.0/3.0; }
        else {
          mask = (1 << d+dimension);
          if ((mask & connect_dir) != 0) { xmin = 2.0/3.0; }
        }
      }
      else {
        mask = (1 << d);
        if ((mask & connect_dir) != 0) { xmax = 2.0/3.0; }
        mask = (1 << (d+dimension));
        if ((mask & connect_dir) != 0) { xmin = 1.0/3.0; }
      }

      x = (x/num_intervals)*(xmax-xmin) + xmin;
      coord[i*dimension+d] = cube_coord[d] + x;
    }
  }
}

/// Position dual isosurface vertices at random location in cube.
/// Don't position near restricted facets.
/// C++ STL vector format for array coord[].
void QDUAL::position_dual_isovertices_random
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const QDUAL_TABLE & qdual_table,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 const RANDOM_SEED_TYPE seed,
 const bool flag_V1w_close,
 const int num_intervals,
 std::vector<COORD_TYPE> & coord)
{
  const int dimension = scalar_grid.Dimension();

  coord.resize(iso_vlist.size()*dimension);
  position_dual_isovertices_random
    (scalar_grid, qdual_table, iso_vlist, seed, flag_V1w_close, num_intervals,
     &(coord.front()));
}
