/// \file qdual_types.h
/// Type definitions for qdual.

/*
  QDUAL: Quality Dual Isosurface Generation
  Copyright (C) 2014-2015 Arindam Bhattacharya, Rephael Wenger

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

#ifndef _QDUAL_TYPES_
#define _QDUAL_TYPES_

#include <string>

#include "ijk.txx"
#include "ijkscalar_grid.txx"
#include "ijkmerge.txx"


namespace NamedConstants{
	static const int DIM3 = 3;
	static const int NUM_CUBE_FACETS = 6;
	static const int NUM_CUBE_FACET_VERT = 4;
	static const int NUM_CUBE_VERT = 8;
	static const int NUM_CUBE_EDGES = 12;
	static const int MAX_NUM_ISOVERT_PER_CUBE = 4;
	static const int NO_ISOVERTEX_IN_CUBE = -1;
	static const int VERT_PER_QUAD = 4;
	static const int VERT_PER_TRI = 3;
}


namespace QDUAL {


// **************************************************
// SCALAR TYPES
// **************************************************

  typedef float SCALAR_TYPE;     ///< Scalar value type.
  typedef float COORD_TYPE;      ///< Isosurface vertex coordinate type.
  typedef int VERTEX_INDEX;      ///< Grid vertex index type.
  typedef float GRADIENT_TYPE;   ///< Gradient coordinate type.
  typedef int GRID_COORD_TYPE;   ///< Grid vertex coordinate type.
  typedef int AXIS_SIZE_TYPE;    ///< Axis size type.
  typedef int ISO_VERTEX_INDEX;  ///< Isosurface vertex index type.
  typedef int MERGE_INDEX;       ///< Merge index type.
  typedef unsigned short DEGREE_TYPE;   ///< DEGREE_PER_VERTEX
  typedef unsigned char DIRECTION_TYPE; ///< Direction type.
  typedef unsigned int RANDOM_SEED_TYPE; ///< Random seed type.
  typedef int QUAD_INDEX; /// iso_quad index type


  /// Edge index type.
  /// Vertex and edge indices must have the same type.
  typedef VERTEX_INDEX EDGE_INDEX;

  /// Facet vertex index type.
  typedef unsigned char FACET_VERTEX_INDEX;

// **************************************************
// ARRAY TYPES
// **************************************************

  typedef std::vector<COORD_TYPE> COORD_ARRAY;   ///< Grid coordinate array.
  typedef std::vector<VERTEX_INDEX>              /// Vertex index array.
    VERTEX_INDEX_ARRAY; 
  typedef std::vector<SCALAR_TYPE> SCALAR_ARRAY; ///< Scalar array.

// **************************************************
// ENUMERATED TYPES
// **************************************************

  /// Interpolation type.
  /// - LINEAR_INTERPOLATION: Determine the location of an isosurface vertex
  ///   using linear interpolation on the endpoints of the cube, pyramid
  ///   or simplex containing the isosurface vertex.
  /// - MULTILINEAR_INTERPOLATION: Determine the location of an 
  ///   isosurface vertex using multilinear interpolation on the cube vertices.
  typedef enum { LINEAR_INTERPOLATION, MULTILINEAR_INTERPOLATION }
  INTERPOLATION_TYPE;

  /// Isosurface vertex position method.
  /// CUBE_CENTER: Position isosurface vertices at cube centers.
  /// CENTROID_EDGE_ISO: Position isosurface vertices at the centroid
  ///                    of the edge isosurface intersections.
  /// RANDOM_POS: Generate random position in cube. (For testing.)
  typedef enum { CUBE_CENTER, CENTROID_EDGE_ISO, RANDOM_POS } 
  VERTEX_POSITION_METHOD;

  /// Quadrilateral triangulation method.
  typedef enum { UNDEFINED_TRI, UNIFORM_TRI, SPLIT_MAX_ANGLE }
  QUAD_TRI_METHOD;

// **************************************************
// CLASSES
// **************************************************

  typedef IJK::BOX<VERTEX_INDEX> GRID_BOX;  ///< Grid box type.

  /// Representation of grid edge by lower/leftmost endpoint and direction.
  struct GRID_EDGE {
    VERTEX_INDEX endpoint0;
    int direction;
  };
}

#endif
