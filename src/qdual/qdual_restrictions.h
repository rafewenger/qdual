#ifndef _QDUAL_RESTRICTIONS_
#define _QDUAL_RESTRICTIONS_
#include <memory>
#include <vector>

#include "qdual_types.h"
#include "qdual_datastruct.h"
#include "ijk.txx"
#include "ijkmesh.txx"
#include "ijkgrid_macros.h"
#include "qdual_table.h"


using namespace QDUAL;
using namespace std;

class RESTRICTED_EDGE
{
public:
	VERTEX_INDEX end0;
	int edge_dir;
};

void compute_restrictions_BList(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,	
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table
	);
void compute_restrictions_CList(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const float e,
	const bool moveVertex,
	const SCALAR_TYPE isovalue,
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	const std::vector<COORD_TYPE> & vertex_coord);

// Compute the sep vert
// note this is a pre-req for compute_restriction_CList
void compute_sep_vert(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	QDUAL_TABLE & qdual_table);

// Compute the sep edges
void compute_sep_edges(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	QDUAL_TABLE & qdual_table);


// Main functions to set restrictions

void set_restrictions(
	const DUALISO_DATA & dualiso_data,
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	std::vector<COORD_TYPE> & vertex_coord,
	DUALISO_INFO & dualiso_info,
	IJK::BOOL_GRID<DUALISO_GRID> &boundary_grid,
	vector<VERTEX_INDEX> &isolatedList);

#endif // !_QDUAL_RESTRICTIONS_
