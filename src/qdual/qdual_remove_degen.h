#ifndef _QDDUAL_REMOVE_DEGENERATE_
#define _QDDUAL_REMOVE_DEGENERATE_
#include <unordered_map>
#include <vector>
#include "qdual_types.h"
#include "qdual_datastruct.h"
#include "ijk.txx"
#include "ijkmesh.txx"
#include "qdual_table.h"
using namespace QDUAL;

namespace QTRIANGULATE{

	//Returns true if degen quads exists
	//Returns the quadvert without the degen QUADS.
	//Returns the track_quad_indices. 
	bool triangulate_degenerate_quads (
		std::vector<VERTEX_INDEX> & quad_vert,
		std::vector<VERTEX_INDEX> & dual_edge,
		std::vector<VERTEX_INDEX> & tri_vert,
		const std::vector<COORD_TYPE> & vertex_coord,
		std::unordered_map<QUAD_INDEX, QUAD_INDEX> &track_quad_indices);

	// Remove degenerate quads
	// Also sets the track_quad_indices.
	void remove_degenerate_quads(
		std::vector<VERTEX_INDEX> & quad_vert,
		std::vector<VERTEX_INDEX> & dual_edge,
		std::vector<VERTEX_INDEX> & non_degen_quad_vert, // only non degenerate quads
		std::vector<VERTEX_INDEX> & non_degen_quad_dual_edge,
		std::vector<VERTEX_INDEX> & tri_vert,
		const std::vector<COORD_TYPE> & vertex_coord,
		//returns
		std::unordered_map<QUAD_INDEX, QUAD_INDEX > &track_quad_indices
		);


	// Triangulate quads based on their angles
	void triangulate_quad_angle_based(
		const DUALISO_SCALAR_GRID_BASE & scalar_grid,
		std::vector<VERTEX_INDEX> & non_degen_quad_vert, // only non degenerate quads
		const std::vector<VERTEX_INDEX> & orgQuadvert,
		std::vector<VERTEX_INDEX> & tri_vert,
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
		const std::vector<COORD_TYPE> & vertex_coord,
		IJK::BOOL_GRID<DUALISO_GRID> &boundary_grid,
		QDUAL_TABLE & qdual_table,
		std::unordered_map<VERTEX_INDEX,VERTEX_INDEX> & diagonalMap,
		const std::vector<DIRECTION_TYPE> &orth_dir,
		std::unordered_map<QUAD_INDEX, QUAD_INDEX > & track_quad_indices,
		IJK::ARRAY<VERTEX_INDEX> &collapse_map,
		bool printInfo
		);

	void triangulate_quad_angle_based(
		const DUALISO_SCALAR_GRID_BASE & scalar_grid,
		std::vector<VERTEX_INDEX> & quadVert, // only non degenerate quads
		std::vector<VERTEX_INDEX> & dual_edge,
		const std::vector<VERTEX_INDEX> & origQuadVert,
		std::vector<VERTEX_INDEX> & tri_vert,
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
		const std::vector<COORD_TYPE> & vertex_coord,
		IJK::BOOL_GRID<DUALISO_GRID> &boundary_grid,
		QDUAL_TABLE & qdual_table,
		//std::unordered_map<VERTEX_INDEX,VERTEX_INDEX> & diagonalMap,
		const std::vector<DIRECTION_TYPE> &orth_dir,
		std::unordered_map< QUAD_INDEX, QUAD_INDEX > & track_quad_indices,
		IJK::ARRAY<VERTEX_INDEX> & origCollapse_map,
		bool printInfo
		);

	//Check if the vertex is a boundary 
	//Param "vertex" index into isovlist
	//returns flag_boundary true if in the boundary
	void isBoundaryIsoVertex(
		const int vertex, //index into isovlist
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
		IJK::BOOL_GRID<DUALISO_GRID> &boundary_grid,
		QDUAL_TABLE & qdual_table,
		bool & flag_boundary
		);


	//Setup the hashmap to find quads which share diagonals
	void hashQuadsDual2GridEdge(
		std::unordered_map<VERTEX_INDEX,VERTEX_INDEX> & diagonalMap,
		std::vector<VERTEX_INDEX> & quad_vert,
		const std::vector<DIRECTION_TYPE> &orth_dir,
		std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
		const std::vector<COORD_TYPE> & vertex_coord);


	/// version B of 
	//Setup the hashmap to find quads which share diagonals
	void hashQuadsDual2GridEdgeB(
		std::unordered_map<VERTEX_INDEX,VERTEX_INDEX> & diagonalMap,
		const std::vector<VERTEX_INDEX> &dual_edge);
}
#endif // !_QDDUAL_REMOVE_DEGENERATE_
