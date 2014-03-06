#ifndef _QDDUAL_REMOVE_DEGENERATE_
#define _QDDUAL_REMOVE_DEGENERATE_
#include <unordered_map>
#include "qdual_types.h"
#include "qdual_datastruct.h"
#include "ijk.txx"
#include "ijkmesh.txx"
#include "qdual_table.h"
using namespace QDUAL;

namespace QTRIANGULATE{
	// Triangulate the non degenerate quads into tris. 
	// Returns an updated set of quads and tris. 
	// Returns a bool which "has_degen_quads"
	bool triangulate_non_degen_quads (
		std::vector<VERTEX_INDEX> & quad_vert,
		std::vector<VERTEX_INDEX> & tri_vert,
		const std::vector<COORD_TYPE> & vertex_coord);

	// Triangulate all quads
	void triangulate_quads (
		const DUALISO_SCALAR_GRID_BASE & scalar_grid,
		std::vector<VERTEX_INDEX> & quad_vert,
		std::vector<VERTEX_INDEX> & tri_vert,
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
		const std::vector<COORD_TYPE> & vertex_coord,
		QDUAL_TABLE & qdual_table,
		IJK::BOOL_GRID<DUALISO_GRID> &boundary_grid,
		const std::vector<DIRECTION_TYPE> &orth_dir);

	// Remove degfenerate quads
	void remove_degenerate_quads(
		std::vector<VERTEX_INDEX> & quad_vert,
		std::vector<VERTEX_INDEX> & non_degen_quad_vert, // only non degenerate quads
		std::vector<VERTEX_INDEX> & tri_vert,
		const std::vector<COORD_TYPE> & vertex_coord
		);


	// Triangulate quads based on their angles
	void triangulate_quad_angle_based(
		const DUALISO_SCALAR_GRID_BASE & scalar_grid,
		std::vector<VERTEX_INDEX> & non_degen_quad_vert, // only non degenerate quads
		std::vector<VERTEX_INDEX> & tri_vert,
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
		const std::vector<COORD_TYPE> & vertex_coord,
		IJK::BOOL_GRID<DUALISO_GRID> &boundary_grid,
		QDUAL_TABLE & qdual_table,
		std::unordered_map<VERTEX_INDEX,VERTEX_INDEX> & diagonalMap,
		const std::vector<DIRECTION_TYPE> &orth_dir
		);

	//Check if the vertex is a boundary 
	//param 1 index into isovlist
	//returns flag_boundary true if in the boundary
	void isBoundaryIsoVertex(
		const int vertex, //index into isovlist
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
		IJK::BOOL_GRID<DUALISO_GRID> &boundary_grid,
		QDUAL_TABLE & qdual_table,
		bool & flag_boundary
		);

	// Compute degree of each vertex
	// Only for non degenerate poly
	// param 1 : num vertex in the poly.
	// param 2 : non degenerate polys
	void compute_degree_per_vertex(
		const int vert_per_poly,
		std::vector<VERTEX_INDEX> & poly_vert, // only non degenerate quads
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist		
		);
	void reset_degree_per_vertex(
		const int vert_per_poly,
		std::vector<VERTEX_INDEX> & poly_vert, // only non degenerate quads
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist		
		);
}
#endif // !_QDDUAL_REMOVE_DEGENERATE_
