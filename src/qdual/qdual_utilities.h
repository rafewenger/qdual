/// \file qdual_utilities.h
/// QDUAL utilities.
/// Various functions for processing data.

/*
  QDUAL: Quality Dual Isosurface Generation
  Copyright (C) 2015 Arindam Bhattacharya, Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _QDUAL_
#define _QDUAL_

#include <string>
#include <memory>

#include "ijk.txx"
#include "qdual_types.h"
#include "qdual_datastruct.h"

namespace QDUAL {

	/// Compute degree of each vertex
	/// @param num_vert_per_poly Number of vertices per polygon.
	/// @param poly_vert List of non degenerate polygons.
	/// @pre Polygons are not degenerate.
	void compute_degree_per_vertex(
		const int vert_per_poly,
		std::vector<VERTEX_INDEX> & poly_vert, // only non degenerate quads
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist		
		);

  /// Reset degree per vertex to 0.
	void reset_degree_per_vertex
    (std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist);
}

#endif
