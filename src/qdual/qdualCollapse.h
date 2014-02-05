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
	//store information regarding collapse
	class COLLAPSE_INFO{
		int num_facet_collapse;
		int num_edge_collapse;
		int num_vertex_collapse;
		std::vector<std::pair <int, int> > edge_pairs;
	};
	// Dual Collapse
	void dual_collapse(
		const DUALISO_SCALAR_GRID_BASE & scalar_grid,
		const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
		std::vector<VERTEX_INDEX> & quad_vert,
		const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
		const std::vector<COORD_TYPE> & vertex_coord,
		const float epsilon
		);
}
#endif // !_QDUAL_COLLAPSE_
