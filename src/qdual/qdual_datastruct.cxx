/// \file qdual_datastruct.cxx
/// qdual data structures

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

#include <assert.h>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "qdual_types.h"
#include "qdual_datastruct.h"

using namespace IJK;
using namespace QDUAL;


// **************************************************
// DUAL CONTOURING ISOSURFACE CLASS
// **************************************************

void DUAL_ISOSURFACE::Clear()
{
  isopoly_vert.clear();
  vertex_coord.clear();
}

// **************************************************
// CLASS DUALISO_DATA_FLAGS
// **************************************************

// Initialize DUALISO_DATA_FLAGS
void DUALISO_DATA_FLAGS::Init()
{
  interpolation_type = LINEAR_INTERPOLATION;
  vertex_position_method = CENTROID_EDGE_ISO;
  flag_separate_neg = true;
  flag_split_non_manifold = true;
  flag_select_split = false;
  allow_multiple_iso_vertices = true;
  max_small_magnitude = 0.0001;
  use_triangle_mesh = false;
  use_quad_tri_mesh = false;
  quad_tri_method = UNDEFINED_TRI;
  flag_NO_collapse = false;
  flag_no_restriction_AB=false;
  flag_no_restriciton_B=false;
  flag_no_restriciton_C=false;
  flag_collapse_info = false;
  flag_collapse_debug = false;
  qdual_epsilon = 0.33;
  flag_move_vertices = false;
  flag_cap_col = true; 
  flag_V1w_close = true;
  flag_delete_isolate = false;

  random_seed = 10;
  random_num_intervals = 1000;
}


/// Set
void DUALISO_DATA_FLAGS::Set(const DUALISO_DATA_FLAGS & data_flags)
{
  *this = data_flags;
}

// **************************************************
// CLASS DUALISO_DATA
// **************************************************

// Initialize DUALISO_DATA
void DUALISO_DATA::Init()
{
  is_scalar_grid_set = false;
}

void DUALISO_DATA::FreeAll()
{
  is_scalar_grid_set = false;
}


// Copy scalar grid
void DUALISO_DATA::CopyScalarGrid(const DUALISO_SCALAR_GRID_BASE & scalar_grid2)
{
  scalar_grid.Copy(scalar_grid2);
  is_scalar_grid_set = true;
}

// Subsample scalar grid
void DUALISO_DATA::SubsampleScalarGrid
(const DUALISO_SCALAR_GRID_BASE & scalar_grid2, const int subsample_resolution)
{
  scalar_grid.Subsample(scalar_grid2, subsample_resolution);
  is_scalar_grid_set = true;
}

// Supersample scalar grid
void DUALISO_DATA::SupersampleScalarGrid
(const DUALISO_SCALAR_GRID_BASE & scalar_grid2, const int supersample_resolution)
{
  scalar_grid.Supersample(scalar_grid2, supersample_resolution);
  is_scalar_grid_set = true;
}

// Copy, subsample or supersample scalar grid.
void DUALISO_DATA::SetScalarGrid
(const DUALISO_SCALAR_GRID_BASE & scalar_grid2, 
 const bool flag_subsample, const int subsample_resolution,
 const bool flag_supersample, const int supersample_resolution)
{
  PROCEDURE_ERROR error("DUALISO_DATA::SetScalarGrid");

  if (flag_subsample && flag_supersample) {
    error.AddMessage("Scalar grid cannot both be subsampled and supersampled.");
    throw error;
  }
  
  if (flag_subsample) {
    // subsample grid
    SubsampleScalarGrid(scalar_grid2, subsample_resolution);
  }
  else if (flag_supersample) {
    // supersample grid
    SupersampleScalarGrid(scalar_grid2, supersample_resolution);
  }
  else {
    CopyScalarGrid(scalar_grid2);
  };
}

/// Set type of interpolation
void DUALISO_DATA::SetInterpolationType
(const INTERPOLATION_TYPE interpolation_type)
{
  this->interpolation_type = interpolation_type;
}

/// Set isosurface vertex position method.
void DUALISO_DATA::SetVertexPositionMethod    
(const VERTEX_POSITION_METHOD vertex_position_method)
{
  this->vertex_position_method = vertex_position_method;
}

/// Set flag for allowing multiple isosurface vertices in a single cube.
void DUALISO_DATA::SetAllowMultipleIsoVertices(const bool flag)
{
  this->allow_multiple_iso_vertices = flag;
}

void DUALISO_DATA::SetConvertQuadToTri
(const QUAD_TRI_METHOD method)
{
  this->use_triangle_mesh = true;
  this->quad_tri_method = method;
}

void DUALISO_DATA::UnsetConvertQuadToTri()
{
  this->use_triangle_mesh = false;
}


/// Check data structure
bool DUALISO_DATA::Check(IJK::ERROR & error) const
{
  IJK::ERROR error2;

  if (!IsScalarGridSet()) {
    error.AddMessage("Scalar grid is not set.");
    return(false);
  }

  return(true);
}


// **************************************************
// DUALISO TIME
// **************************************************

QDUAL::DUALISO_TIME::DUALISO_TIME()
// constructor
{
  Clear();
}

void QDUAL::DUALISO_TIME::Clear()
{
  preprocessing = 0.0;
  extract = 0.0;
  merge = 0.0;
  position = 0.0;
  total = 0.0;
}

void QDUAL::DUALISO_TIME::Add(const DUALISO_TIME & dualiso_time)
{
  preprocessing += dualiso_time.preprocessing;
  extract += dualiso_time.extract;
  merge += dualiso_time.merge;
  position += dualiso_time.position;
  total += dualiso_time.total;
}

// **************************************************
// INFO CLASSES
// **************************************************

QDUAL::GRID_INFO::GRID_INFO()
{
  Clear();
}

void QDUAL::GRID_INFO::Clear()
{
  num_cubes = 0;
}

void QDUAL::SCALAR_INFO::Init(const int dimension)
{
  this->dimension = 0;
  SetDimension(dimension);

  Clear();
}

void QDUAL::SCALAR_INFO::FreeAll()
{
  dimension = 0;
}

void QDUAL::SCALAR_INFO::Clear()
{
  num_non_empty_cubes = 0;
  num_bipolar_edges = 0;
}

void QDUAL::SCALAR_INFO::SetDimension(const int dimension)
{
  FreeAll();

  this->dimension = dimension;

  Clear();
}

void QDUAL::SCALAR_INFO::Copy(const SCALAR_INFO & info)
{
  Init(info.Dimension());
  num_non_empty_cubes = info.num_non_empty_cubes;
  num_bipolar_edges = info.num_bipolar_edges;
}

/// Copy assignment.
const SCALAR_INFO &  QDUAL::SCALAR_INFO::operator =
(const SCALAR_INFO & right)
{
  if (&right != this) {
    FreeAll();
    Copy(right);
  }

  return *this;
}

QDUAL::SCALAR_INFO::~SCALAR_INFO()
{
  dimension = 0;
  Clear();
}

QDUAL::MULTI_ISOV_INFO::MULTI_ISOV_INFO()
{
  Clear();
}

void QDUAL::MULTI_ISOV_INFO::Clear()
{
  num_cubes_single_isov = 0;
  num_cubes_multi_isov = 0;
  num_non_manifold_split = 0;
  num_1_2_change = 0;
}

QDUAL::DUALISO_INFO::DUALISO_INFO()
{
  Clear();
}

QDUAL::DUALISO_INFO::DUALISO_INFO(const int dimension):scalar(dimension)
{
  Clear();
}

void QDUAL::DUALISO_INFO::Clear()
{
  grid.Clear();
  scalar.Clear();
  time.Clear();
}

// **************************************************
// MERGE DATA
// **************************************************

void QDUAL::MERGE_DATA::Init
(const int dimension, const AXIS_SIZE_TYPE * axis_size,
 const MERGE_INDEX num_obj_per_vertex, const MERGE_INDEX num_obj_per_edge)
{
  this->num_obj_per_vertex = num_obj_per_vertex;
  this->num_obj_per_edge = num_obj_per_edge;
  this->num_obj_per_grid_vertex = 
    dimension*num_obj_per_edge + num_obj_per_vertex;
  compute_num_grid_vertices(dimension, axis_size, num_vertices);
  num_edges = dimension*num_vertices;
  vertex_id0 = num_obj_per_edge*num_edges;
  MERGE_INDEX num_obj = 
    num_obj_per_vertex*num_vertices + num_obj_per_edge*num_edges;
  INTEGER_LIST<MERGE_INDEX,MERGE_INDEX>::Init(num_obj);
}

bool QDUAL::MERGE_DATA::Check(ERROR & error) const
{
  if (MaxNumInt() < 
      NumObjPerVertex()*NumVertices() + NumObjPerEdge()*NumEdges()) {
    error.AddMessage("Not enough allocated memory.");
    return(false);
  };

  return(true);
}

