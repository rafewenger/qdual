#ifndef _QDUAL_COLLAPSE_
#define _QDUAL_COLLAPSE_

#include <string>

#include "ijk.txx"
#include "qdual_types.h"
#include "qdual_datastruct.h"

#include "ijkdualtable.h"
#include "ijkdualtable_ambig.h"

using namespace QDUAL;

namespace QCOLLAPSE{
	// FACET associated with a cube index and an axis
	// the facet is orthogonal to the axis
	typedef  std::pair <int,int> FACET;

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
