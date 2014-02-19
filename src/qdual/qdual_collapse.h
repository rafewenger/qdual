#ifndef _QDUAL_COLLAPSE_
#define _QDUAL_COLLAPSE_

#include <string>
#include <vector>
#include <utility>
#include <cmath> 
#include <algorithm>
#include "ijk.txx"
#include "ijkmesh.txx"
#include "qdual_types.h"
#include "qdual_datastruct.h"

#include "ijkdualtable.h"
#include "ijkdualtable_ambig.h"

using namespace QDUAL;

namespace QCOLLAPSE{

	enum COLLAPSE_TYPE {FACET, EDGE, VERTEX};
	//store information regarding collapse
	class COLLAPSE_INFO{
		int num_facet_collapse;
		int num_edge_collapse;
		int num_vertex_collapse;
		std::vector<std::pair <int, int> > edge_pairs;
	};
	// Dual Collapse
	void dual_collapse(
		const DUALISO_DATA & dualiso_data,
		const DUALISO_SCALAR_GRID_BASE & scalar_grid,
		std::vector<VERTEX_INDEX> & quad_vert,
		const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
		const std::vector<COORD_TYPE> & vertex_coord,
		const float epsilon
		);


	// Triangulate quads based on their angles
	void triangulate_quad_angle_based(
		std::vector<VERTEX_INDEX> & non_degen_quad_vert, // only non degenerate quads
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
		const std::vector<COORD_TYPE> & vertex_coord
		);

}
#endif // !_QDUAL_COLLAPSE_
