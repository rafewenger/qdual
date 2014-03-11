#ifndef _QDUAL_COLLAPSE_
#define _QDUAL_COLLAPSE_

#include <string>
#include <vector>
#include <utility>
#include <cmath> 
#include <algorithm>

#include "ijk.txx"
#include "ijkmesh.txx"
#include "ijktime.txx"
#include "qdual_types.h"
#include "qdual_datastruct.h"
#include "qdual_restrictions.h"

#include "ijkdualtable.h"
#include "ijkdualtable_ambig.h"


namespace QCOLLAPSE{

	enum COLLAPSE_TYPE {FACET, EDGE, VERTEX};
	//store information regarding collapse
	class COLLAPSE_DEBUG{
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
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
		std::vector<COORD_TYPE> & vertex_coord,
        const std::vector<DIRECTION_TYPE> & orth_dir,
		const float epsilon,
		IJK::ARRAY<VERTEX_INDEX> &collapse_map,
		DUALISO_INFO & dualiso_info
		);

	//Setup the COLLAPSE MAP
	//param 1 : collapse map
	//param 2 : size of the map
	void setup_collapse_map(
		IJK::ARRAY<VERTEX_INDEX> &collapse_map,
		const int num_vertex);

	//Update vertex from the collapse map
	int find_vertex(
		IJK::ARRAY<VERTEX_INDEX> &collapse_map,
		int endpt);
    
	///Delete isolated vertices
	void delIsolated(
		std::vector<VERTEX_INDEX> & quad_vert,
		vector<VERTEX_INDEX> isolatedList,
		const DUALISO_SCALAR_GRID_BASE & scalar_grid,
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
		DUALISO_INDEX_GRID & first_isov,
		IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
        bool printInfo
		);
}
#endif // !_QDUAL_COLLAPSE_
