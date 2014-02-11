#ifndef _QDDUAL_REMOVE_DEGENERATE_
#define _QDDUAL_REMOVE_DEGENERATE_

#include "qdual_types.h"
#include "qdual_datastruct.h"
#include "ijk.txx"
#include "ijkmesh.txx"
using namespace QDUAL;

namespace QTRIANGULATE{
	// Triangulate the non degenerate quads into tris. 
	// Returns an updated set of quads and tris. 
	void triangulate_non_degen_quads (
		std::vector<VERTEX_INDEX> & quad_vert,
		std::vector<VERTEX_INDEX> & tri_vert,
		const std::vector<COORD_TYPE> & vertex_coord);

	// Triangulate all quads
	void triangulate_quads (
		std::vector<VERTEX_INDEX> & quad_vert,
		std::vector<VERTEX_INDEX> & tri_vert,
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
		const std::vector<COORD_TYPE> & vertex_coord);

	// Remove degfenerate quads
	void remove_degenerate_quads(
		std::vector<VERTEX_INDEX> & quad_vert,
		std::vector<VERTEX_INDEX> & non_degen_quad_vert, // only non degenerate quads
		std::vector<VERTEX_INDEX> & tri_vert,
		const std::vector<COORD_TYPE> & vertex_coord
		);


	// Triangulate quads based on their angles
	void triangulate_quad_angle_based(
		std::vector<VERTEX_INDEX> & non_degen_quad_vert, // only non degenerate quads
		std::vector<VERTEX_INDEX> & tri_vert,
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
		const std::vector<COORD_TYPE> & vertex_coord
		);
}
#endif // !_QDDUAL_REMOVE_DEGENERATE_
