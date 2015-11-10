/// \file qdual_utilities.cxx
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

#include "qdual_utilities.h"

// Compute degree of each vertex
// @param num_vert_per_poly Number of vertices per polygon.
// @param poly_vert List of non degenerate polygons.
// @pre Polygons are not degenerate.
void QDUAL::compute_degree_per_vertex(
	const int vert_per_poly,
	std::vector<VERTEX_INDEX> & poly_vert, // only non degenerate quads
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist		
	)
{
	const int num_poly = poly_vert.size()/vert_per_poly;
	for (int q=0; q<num_poly; q++)
	{
		for (int d=0;d<vert_per_poly;d++){
			VERTEX_INDEX v = poly_vert[vert_per_poly*q+d];
			iso_vlist[v].ver_degree++;
		}
	}
}

void QDUAL::reset_degree_per_vertex
(std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist)
{
  for (int i = 0; i < iso_vlist.size(); i++) 
    { iso_vlist[i].ver_degree = 0; }
}

