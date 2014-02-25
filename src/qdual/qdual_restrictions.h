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

class RestricitonInfo
{
public:
	int restriction_BList_size;
	int restriction_CList_size;

	vector<VERTEX_INDEX> restricted_vertex_info;
	vector< pair<VERTEX_INDEX, int> > restricted_edges_info;

RestricitonInfo()
{
	restriction_CList_size = 0;
	restriction_BList_size = 0;
}
~RestricitonInfo(){};

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



// Main functions to set restrictions

void set_restrictions(
	const DUALISO_DATA & dualiso_data,
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	const std::vector<COORD_TYPE> & vertex_coord,
	std::unique_ptr<RestricitonInfo> &rs_info);

#endif // !_QDUAL_RESTRICTIONS_
