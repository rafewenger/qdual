#include "qdual_restrictions.h"
#include "ijkscalar_grid.txx"


using namespace std;
using namespace IJK;
using namespace NamedConstants;

bool moveAwayFromFacetF(
	const int vertex,
	const int f, 
	const float e, //epsilon
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	std::vector<COORD_TYPE> & vertex_coord
	)
{
	bool isMoved = false;
	int orthoDir = (f + DIM3)%DIM3;
	float vCoordOrthoDir = vertex_coord[DIM3*vertex+ orthoDir];
	float facetCoordOrthoDir = iso_vlist[vertex].cube_coord[orthoDir];

	if (f > DIM3-1)
	{
		facetCoordOrthoDir++;
		if ((facetCoordOrthoDir-vCoordOrthoDir)<e)
		{
			vertex_coord[DIM3*vertex+orthoDir] = facetCoordOrthoDir-e;
			isMoved = true;
		}
	}
	else
	{
		if ((vCoordOrthoDir-facetCoordOrthoDir)<e)
		{
			vertex_coord[DIM3*vertex+orthoDir] = facetCoordOrthoDir+e;
			isMoved = true;
		}
	}
	return isMoved;
}
// Move vertices away from vertex
// vertex: Index into vertex coord
// e: distance to move vertex by
void moveAwayFromVertex(
	const int vertex,
	const float e, //epsilon
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	std::vector<COORD_TYPE> & vertex_coord, 
	DUALISO_INFO & dualiso_info
	)
{
	for (int orthoDir = 0; orthoDir < NUM_CUBE_FACETS; orthoDir++)
	{
		bool  isMoved = moveAwayFromFacetF(vertex,orthoDir,e,iso_vlist, vertex_coord);
		if (isMoved)
			dualiso_info.mv_info.moveFromVertices++;
	}
}

//Move away from the edge given by the facets f1,f2
//of the cube containing the "vertex"
void moveAwayFromEdge(
	const int vertex,
	const int f1,
	const int f2,
	const float e,//epsilon
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	std::vector<COORD_TYPE> & vertex_coord,
	DUALISO_INFO & dualiso_info
	)
{
	bool isMoved = moveAwayFromFacetF(vertex, f1,
		e, iso_vlist, vertex_coord);
	if (isMoved)
		dualiso_info.mv_info.moveFromEdges++;

	isMoved =  moveAwayFromFacetF(vertex, f2,
		e, iso_vlist, vertex_coord);
	if (isMoved)
		dualiso_info.mv_info.moveFromEdges++;
}

//Subfunction to compute_restrictions_BList (the faster version)
bool doesLoopSurroundEdge (	
	const VERTEX_INDEX end0,
	const int edge_dir,
	const SCALAR_TYPE isovalue,
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const DUALISO_INDEX_GRID & first_isov,
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	const IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	const QDUAL_TABLE & qdual_table,
	const int c0, //index of the cube  left and below the edge
	const int d1, 
	const int d2
	)
{

	static int d1coef[NUM_CUBE_FACET_VERT] = {0,1,0,1};
	static int d2coef[NUM_CUBE_FACET_VERT] = {0,0,1,1};

	VERTEX_INDEX end1 = scalar_grid.NextVertex(end0,edge_dir);
	SCALAR_TYPE s_end0, s_end1;
	s_end0 = scalar_grid.Scalar(end0);
	s_end1 = scalar_grid.Scalar(end1);

	if (((s_end0 > isovalue) && ( s_end1 < isovalue)) ||
		((s_end0 < isovalue) && ( s_end1 > isovalue)))
	{	return false; }
	int num_connect=0;
	for (int i = 0; i<NUM_CUBE_FACET_VERT; i++)
	{
		VERTEX_INDEX c = c0 + d1coef[i] * scalar_grid.AxisIncrement(d1) +
			d2coef[i]*scalar_grid.AxisIncrement(d2);

		VERTEX_INDEX indx_iso_vlist = first_isov.Scalar(c);

		if (indx_iso_vlist == NO_ISOVERTEX_IN_CUBE)
		{
			return false;
		}
		bool flag_connect = false;

		VERTEX_INDEX f1 = (d1+ (1-d1coef[i])* DIM3) % NUM_CUBE_FACETS;
		VERTEX_INDEX f2 = (d2+ (1-d2coef[i])* DIM3) % NUM_CUBE_FACETS;

		int table_ind = iso_vlist[indx_iso_vlist].table_index;
		int num_iso_verts = isodual_table.NumIsoVertices(table_ind);
		//for each isosurface vertex in cube c
		for (int j = 0; j < num_iso_verts; j++)
		{
			IJKDUALTABLE::TABLE_INDEX it = iso_vlist[indx_iso_vlist+j].table_index;
			FACET_VERTEX_INDEX ip = iso_vlist[indx_iso_vlist+j].patch_index;
			QDUAL_TABLE::DIR_BITS edge_flag = qdual_table.connectDir(it, ip);
			bool cond1 = (edge_flag & (1<<f1));
			bool cond2 = (edge_flag & (1<<f2));

			if (cond1 && cond2)
			{
				flag_connect = true;
			}
		}
		if (flag_connect) 
		{
			num_connect++;
		}

	}
	if(num_connect >=3)
		return true;
	else 
		return false;
}

// With added Check 
bool check_edge_has_square_isosurface_path(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const DUALISO_INDEX_GRID & first_isov,
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	const IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	const QDUAL_TABLE & qdual_table,
	const int c0,
	const int d1, 
	const int d2
	)
{
	static int d1coef[NUM_CUBE_FACET_VERT] = {0,1,0,1};
	static int d2coef[NUM_CUBE_FACET_VERT] = {0,0,1,1};
	int num_connect=0;
	/// Check if the edge is a probable candidate
	for (int i = 0; i<NUM_CUBE_FACET_VERT; i++)
	{
		VERTEX_INDEX c = c0 + d1coef[i] * scalar_grid.AxisIncrement(d1) +
			d2coef[i]*scalar_grid.AxisIncrement(d2);
		VERTEX_INDEX indx_iso_vlist = first_isov.Scalar(c);
		if (indx_iso_vlist == NO_ISOVERTEX_IN_CUBE)
		{
			return false;
		}
	}

	for (int i = 0; i<NUM_CUBE_FACET_VERT; i++)
	{
		VERTEX_INDEX c = c0 + d1coef[i] * scalar_grid.AxisIncrement(d1) +
			d2coef[i]*scalar_grid.AxisIncrement(d2);
		VERTEX_INDEX indx_iso_vlist = first_isov.Scalar(c);
		bool flag_connect = false;

		VERTEX_INDEX f1 = (d1+ (1-d1coef[i])* DIM3) % NUM_CUBE_FACETS;
		VERTEX_INDEX f2 = (d2+ (1-d2coef[i])* DIM3) % NUM_CUBE_FACETS;

		int table_ind = iso_vlist[indx_iso_vlist].table_index;
		int num_iso_verts = isodual_table.NumIsoVertices(table_ind);
		//for each isosurface vertex in cube c
		for (int j = 0; j < num_iso_verts; j++)
		{
			IJKDUALTABLE::TABLE_INDEX it = iso_vlist[indx_iso_vlist+j].table_index;
			FACET_VERTEX_INDEX ip = iso_vlist[indx_iso_vlist+j].patch_index;
			QDUAL_TABLE::DIR_BITS edge_flag = qdual_table.connectDir(it, ip);
			bool cond1 = (edge_flag & (1<<f1));
			bool cond2 = (edge_flag & (1<<f2));

			if (cond1 && cond2)
			{
				flag_connect = true;
			}
		}
		if (flag_connect) 
		{
			num_connect++;
		}

	}
	if(num_connect >=3)
		return true;
	else 
		return false;
}

//// Square isosurface paths // **** WORKING SLOW CODE **** 
void compute_restrictions_BList_slow(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	vector<RESTRICTED_EDGE> &restricted_edges, 
	const bool print_info)
{
	// Restriction B list is the list of positive and negative interior grid edges
	// surrounded by a square isourface path
	VERTEX_INDEX end0=0;

	for (int edge_dir = 0; edge_dir < DIM3; edge_dir++) {

		int d1 = (edge_dir+1)% DIM3;
		int d2 = (edge_dir+2)%DIM3;
		VERTEX_INDEX c0_increment = 
			scalar_grid.AxisIncrement(d1) + scalar_grid.AxisIncrement(d2);

		IJK_FOR_EACH_INTERIOR_GRID_EDGE_IN_DIRECTION
			(end0, edge_dir, scalar_grid, VERTEX_INDEX) {

				// c0 is the index of the cube left and below to the edge
				VERTEX_INDEX c0 = end0 - c0_increment;
				VERTEX_INDEX index_iso_vlist_c0 = first_isov.Scalar(c0);
				if (index_iso_vlist_c0 == NO_ISOVERTEX_IN_CUBE ) 
				{
					continue;
				}

				VERTEX_INDEX end1 = scalar_grid.NextVertex(end0,edge_dir);
				SCALAR_TYPE s_end0, s_end1;
				s_end0 = scalar_grid.Scalar(end0);
				s_end1 = scalar_grid.Scalar(end1);

				if (((s_end0 > isovalue) && ( s_end1 < isovalue)) ||
					((s_end0 < isovalue) && ( s_end1 > isovalue)))
				{	continue; }

				bool flag_edge = 
					check_edge_has_square_isosurface_path
					(scalar_grid, first_isov, iso_vlist, isodual_table, 
					qdual_table, c0, d1, d2);

				if (flag_edge) {
					RESTRICTED_EDGE re;
					re.end0 = end0;
					re.edge_dir = edge_dir;
					restricted_edges.push_back(re);
				}
		}
	}
}

// Faster implementation of computing restrictions Blist 
void compute_restrictions_BList(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	vector<RESTRICTED_EDGE> &restricted_edges, 
	const bool print_info)
{
	restricted_edges.clear();
	VERTEX_INDEX iv;

	IJK_FOR_EACH_INTERIOR_GRID_VERTEX
		(iv, scalar_grid, VERTEX_INDEX)
	{
		VERTEX_INDEX indx_iso_vlist_iv = first_isov.Scalar(iv);
		if (indx_iso_vlist_iv != -1)
		{
			for (int edge_dir = 0; edge_dir < DIM3; edge_dir++)
			{
				int d1 = (edge_dir+1)% DIM3;
				int d2 = (edge_dir+2)%DIM3;
				VERTEX_INDEX c0_increment = 
					scalar_grid.AxisIncrement(d1) + scalar_grid.AxisIncrement(d2);

				VERTEX_INDEX c0 = iv - c0_increment;				
				VERTEX_INDEX index_iso_vlist_c0 = first_isov.Scalar(c0);

				if (index_iso_vlist_c0 == NO_ISOVERTEX_IN_CUBE ) 
				{
					continue;
				}
				bool flag_edge = 
					doesLoopSurroundEdge
					(iv, edge_dir, isovalue, scalar_grid, first_isov, iso_vlist, isodual_table, 
					qdual_table, c0, d1, d2);
				if (flag_edge) {
					RESTRICTED_EDGE re;
					re.end0 = iv;
					re.edge_dir = edge_dir;
					restricted_edges.push_back(re);
				}
			}
		}
	}

	for (int edgeDir = 0; edgeDir < DIM3; edgeDir++)
	{
		IJK_FOR_EACH_VERTEX_IN_GRID_FACET_INTERIOR
			(iv, edgeDir, scalar_grid, VERTEX_INDEX)
		{
			int d1 = (edgeDir+1)% DIM3;
			int d2 = (edgeDir+2)%DIM3;
			VERTEX_INDEX c0_increment = 
				scalar_grid.AxisIncrement(d1) + scalar_grid.AxisIncrement(d2);
			VERTEX_INDEX c0 = iv - c0_increment;
			bool flag_edge = 
				doesLoopSurroundEdge
				(iv, edgeDir, isovalue, scalar_grid, first_isov, iso_vlist, isodual_table, 
				qdual_table, c0, d1, d2);
			if (flag_edge) {
				RESTRICTED_EDGE re;
				re.end0 = iv;
				re.edge_dir = edgeDir;
				restricted_edges.push_back(re);
			}
		}
	}
}




// Return true if the vertex is surrounded by a box
bool  check_vertex_has_box_around (
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const int iv,
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	int  & count
	)
{
	QDUAL_TABLE::DIR_BITS edgeFlag=0;
	VERTEX_INDEX iv0 = iv - scalar_grid.CubeVertexIncrement(NUM_CUBE_VERT-1);


	for (int j = 0; j < NUM_CUBE_VERT; j++)
	{
		VERTEX_INDEX icube = scalar_grid.CubeVertex(iv0,j);
		VERTEX_INDEX indx_iso_vlist = first_isov.Scalar(icube);
		if ((indx_iso_vlist == NO_ISOVERTEX_IN_CUBE))
			return false;
		else
		{
			int table_ind = iso_vlist[indx_iso_vlist].table_index;
			int num_iso_verts = isodual_table.NumIsoVertices(table_ind);

			//for each isosurface vertex in cube c
			for (int k = 0; k < num_iso_verts; k++)
			{
				if ( iso_vlist[indx_iso_vlist+k].sep_vert == iv)
				{
					count++;
					IJKDUALTABLE::TABLE_INDEX it = iso_vlist[indx_iso_vlist+k].table_index;
					unsigned char ip = iso_vlist[indx_iso_vlist+k].patch_index;
					QDUAL_TABLE::DIR_BITS w_edge_flag = qdual_table.connectDir(it,ip);
					edgeFlag = edgeFlag | w_edge_flag;
				}
			}
		}
	}
	if (count >=3 && edgeFlag == 63)
	{
		return true;
	}
	return false;
}


// Return true if the vertex is surrounded by a box
bool  check_vertex_has_box_around (
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const int iv,
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	vector<VERTEX_INDEX> &isolatedList
	)
{
	QDUAL_TABLE::DIR_BITS edgeFlag=0;

	int count = 0;
	VERTEX_INDEX iv0 = iv - scalar_grid.CubeVertexIncrement(NUM_CUBE_VERT-1);


	for (int j = 0; j < NUM_CUBE_VERT; j++)
	{
		VERTEX_INDEX icube = scalar_grid.CubeVertex(iv0,j);
		VERTEX_INDEX indx_iso_vlist = first_isov.Scalar(icube);

		if (indx_iso_vlist != NO_ISOVERTEX_IN_CUBE)
		{
			int table_ind = iso_vlist[indx_iso_vlist].table_index;
			int num_iso_verts = isodual_table.NumIsoVertices(table_ind);

			//for each isosurface vertex in cube c
			for (int k = 0; k < num_iso_verts; k++)
			{
				if ( iso_vlist[indx_iso_vlist+k].sep_vert == iv)
				{
					count++;
					IJKDUALTABLE::TABLE_INDEX it = iso_vlist[indx_iso_vlist+k].table_index;
					unsigned char ip = iso_vlist[indx_iso_vlist+k].patch_index;
					QDUAL_TABLE::DIR_BITS w_edge_flag = qdual_table.connectDir(it,ip);
					edgeFlag = edgeFlag | w_edge_flag;
				}
			}
		}
	}
	if (count == 8 && edgeFlag == 63)
	{
		isolatedList.push_back(iv);
	}
	if (count >=3 && edgeFlag == 63)
	{
		return true;
	}
	return false;
}

// compute the vertex list which have isosurface boxes around it 
void compute_restrictions_CList(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	const std::vector<COORD_TYPE> & vertex_coord,
	vector<VERTEX_INDEX> &restriction_Clist,
	const IJK::BOOL_GRID<DUALISO_GRID> &boundary_grid,
	vector<VERTEX_INDEX> &isolatedList)
{
	int numv =  scalar_grid.NumVertices();

	for (VERTEX_INDEX iv = 0; iv < numv; iv++)
	{
		if (!boundary_grid.Scalar(iv)) 
		{
			VERTEX_INDEX indx_iso_vlist_iv = first_isov.Scalar(iv);
			if (indx_iso_vlist_iv != -1)
			{
				int count=0;	
				bool flag_box = check_vertex_has_box_around 
					(scalar_grid, iv, iso_vlist, isodual_table,
					first_isov, qdual_table, count);

				if (flag_box){
					// restricted vertex
					restriction_Clist.push_back(iv);
					if (count==8)
					{
						isolatedList.push_back(iv);
					}
				}
			}
		}
	}

}

void compute_sep_vert(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	QDUAL_TABLE & qdual_table)
{
	GRID_COORD_TYPE * c = new GRID_COORD_TYPE[3];

	IJKDUALTABLE::TABLE_INDEX it;
	unsigned char ip;
	int i=0;
	for (i = 0; i < iso_vlist.size(); i++)
	{
		it = iso_vlist[i].table_index;
		ip = iso_vlist[i].patch_index;
		VERTEX_INDEX icube = iso_vlist[i].cube_index;
		VERTEX_INDEX icorner = qdual_table.V1(it,ip);

		if (icorner !=255 )
		{
			iso_vlist[i].sep_vert = scalar_grid.CubeVertex(icube, icorner);
		}
		else
		{
			iso_vlist[i].sep_vert = scalar_grid.NumVertices();
		}
	}
}

void compute_sep_edges(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	QDUAL_TABLE & qdual_table)
{
	GRID_COORD_TYPE * c = new GRID_COORD_TYPE[3];

	IJKDUALTABLE::TABLE_INDEX it;
	unsigned char ip;
	int i=0;
	for (i = 0; i < iso_vlist.size(); i++)
	{
		it = iso_vlist[i].table_index;
		ip = iso_vlist[i].patch_index;
		VERTEX_INDEX icube = iso_vlist[i].cube_index;
    GRID_EDGE cube_edge = qdual_table.V2(it,ip);

		if (cube_edge.endpoint0 !=255 )
		{
      VERTEX_INDEX endpoint0 =
        scalar_grid.CubeVertex(icube, cube_edge.endpoint0);
			iso_vlist[i].sep_edge = endpoint0*DIM3 + cube_edge.direction;
		}
		else
		{
			iso_vlist[i].sep_edge = scalar_grid.NumVertices()*DIM3;
		}
	}
}


//
void set_restrictionsA(
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	QDUAL_TABLE & qdual_table
	)
{

	for (int i = 0; i < iso_vlist.size(); i++)
	{
		for (int j=0; j < NUM_CUBE_FACETS-1; j++)
		{
			IJKDUALTABLE::TABLE_INDEX it = iso_vlist[i].table_index;
			unsigned char ip = iso_vlist[i].patch_index;
			QDUAL_TABLE::DIR_BITS edge_flag = qdual_table.connectDir(it,ip);
			if ((edge_flag & (1<<j)) != 0 )
			{
				int f2 = (j+DIM3) % NUM_CUBE_FACETS;
				iso_vlist[i].restricted_facets = iso_vlist[i].restricted_facets | ((1<<f2));
			}
		}
	}
}

// Set Restrictions B
void set_restrictionsB(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	std::vector<COORD_TYPE> & vertex_coord,
	const float e, // how far to move vertices away from edges and vertices.
	const SCALAR_TYPE isovalue,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	bool print_info,
	bool move_vertex,
	bool move_vertex2, // move vertices by epsilon by 2
	DUALISO_INFO & dualiso_info)
{
	vector<RESTRICTED_EDGE> restricted_edges;
	compute_restrictions_BList( scalar_grid, isovalue, iso_vlist, 
		isodual_table, first_isov, qdual_table, restricted_edges, print_info);
	//store information 
	dualiso_info.rs_info.restriction_BList_size = restricted_edges.size(); 
	//dualiso_info.rs_info.restricted_edges_info = restricted_edges;

	int d1coef[NUM_CUBE_FACET_VERT]={0,1,0,1};
	int d2coef[NUM_CUBE_FACET_VERT]= {0,0,1,1};

	for (int i=0;i < restricted_edges.size(); i++)
	{
		VERTEX_INDEX end0 = restricted_edges[i].end0;
		int edge_dir = restricted_edges[i].edge_dir;

		int d1 = (edge_dir + 1) % DIM3;
		int d2 = (edge_dir + 2) % DIM3;
		VERTEX_INDEX c0 = end0 - scalar_grid.AxisIncrement(d1) - scalar_grid.AxisIncrement(d2);


		VERTEX_INDEX  end1 = scalar_grid.NextVertex(end0, edge_dir);
		if(print_info)
		{
			GRID_COORD_TYPE * coord = new GRID_COORD_TYPE[3];
			scalar_grid.ComputeCoord(end0,coord);
			cout <<"\nRestricted edge "<<end0<<"("<< coord[0]<<" "<<coord[1]<<" "<<coord[2]<<") ";
			scalar_grid.ComputeCoord(end1,coord);
			cout <<end1<<" ("<< coord[0]<<" "<<coord[1]<<" "<<coord[2]<<")"<<endl;
		}
		for (int j=0; j < NUM_CUBE_FACET_VERT; j++)
		{
			VERTEX_INDEX c = c0 + d1coef[j]*scalar_grid.AxisIncrement(d1) 
				+ d2coef[j]*scalar_grid.AxisIncrement(d2);

			VERTEX_INDEX f1 = (d1 + (1-d1coef[j])*DIM3)% NUM_CUBE_FACETS;
			VERTEX_INDEX f2 = (d2 + (1-d2coef[j])*DIM3)% NUM_CUBE_FACETS;

			bool flag_connect = false;

			VERTEX_INDEX indx_iso_vlist = first_isov.Scalar(c);
			if (indx_iso_vlist != -1){
				int table_ind = iso_vlist[indx_iso_vlist].table_index;
				int num_iso_verts = isodual_table.NumIsoVertices(table_ind);
				//for each isosurface vertex in cube c
				for (int k = 0; k < num_iso_verts; k++)
				{
					IJKDUALTABLE::TABLE_INDEX it = iso_vlist[indx_iso_vlist+k].table_index;
					FACET_VERTEX_INDEX ip = iso_vlist[indx_iso_vlist+k].patch_index;
					QDUAL_TABLE::DIR_BITS edge_flag = qdual_table.connectDir(it, ip);
					bool cond1 = (edge_flag & (1<<f1));
					bool cond2 = (edge_flag & (1<<f2));
					if ( cond1 && cond2 )
					{
						if (print_info){
							cout <<"	org vertex: "<<vertex_coord[3*(indx_iso_vlist+k)]<<" "
								<<vertex_coord[3*(indx_iso_vlist+k)+1]<<" "
								<<vertex_coord[3*(indx_iso_vlist+k)+2]<<" ";
						}

						if (move_vertex)
						{
							moveAwayFromEdge(indx_iso_vlist+k, f1, f2,
								e, iso_vlist, vertex_coord, dualiso_info);
						}
						else if (move_vertex2)//move by epsilon/2.0
						{
							moveAwayFromEdge(indx_iso_vlist+k, f1, f2,
								e/2.0, iso_vlist, vertex_coord, dualiso_info);
						}

						if (print_info){
							cout <<"	new vertex: "<<vertex_coord[3*(indx_iso_vlist+k)]<<" "
								<<vertex_coord[3*(indx_iso_vlist+k)+1]<<" "
								<<vertex_coord[3*(indx_iso_vlist+k)+2]<<" \n";
						}
						iso_vlist[indx_iso_vlist+k].restricted_facets = 
							(iso_vlist[indx_iso_vlist+k].restricted_facets | ((1<<f1) | (1<<f2)));
					}
				}

			}
		}
	}
}


// Set restrictions C 
void set_restrictionsC(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const float e, //epsilon to move vertices.
	const bool moveVertex,
	const bool moveVertex2, // move vertices by epsilon/2
	const SCALAR_TYPE isovalue,
	const std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	std::vector<COORD_TYPE> & vertex_coord,
	DUALISO_INFO & dualiso_info,
	const IJK::BOOL_GRID<DUALISO_GRID> &boundary_grid,
	vector<VERTEX_INDEX> &isolatedList)
{
	vector<VERTEX_INDEX> restriction_Clist;

	compute_restrictions_CList( scalar_grid, isovalue, iso_vlist, isodual_table, first_isov,
		qdual_table, vertex_coord, restriction_Clist, boundary_grid, isolatedList);

	dualiso_info.rs_info.restriction_CList_size = restriction_Clist.size();
	dualiso_info.rs_info.restricted_vertex_info = restriction_Clist;

	for (int v = 0; v < restriction_Clist.size(); v++)
	{
		VERTEX_INDEX iv0 = restriction_Clist[v] - scalar_grid.CubeVertexIncrement(NUM_CUBE_VERT-1);

		for (int j = 0; j < NUM_CUBE_VERT; j++)
		{
			VERTEX_INDEX c = scalar_grid.CubeVertex(iv0,j);

			VERTEX_INDEX indx_iso_vlist = first_isov.Scalar(c);
			if (indx_iso_vlist !=-1)
			{
				int table_ind = iso_vlist[indx_iso_vlist].table_index;
				int num_iso_verts = isodual_table.NumIsoVertices(table_ind);
				//for each isosurface vertex in cube c
				for (int k = 0; k < num_iso_verts; k++)
				{
					if ( iso_vlist[indx_iso_vlist+k].sep_vert == restriction_Clist[v])
					{
						iso_vlist[indx_iso_vlist+k].flag_restrictionC = true;
						if (moveVertex)
						{
							moveAwayFromVertex(indx_iso_vlist+k,
								e, iso_vlist, vertex_coord, dualiso_info);
						}
						else if (moveVertex2)
						{
							moveAwayFromVertex(indx_iso_vlist+k,
								e/2.0, iso_vlist, vertex_coord, dualiso_info);
						}
					}
				}
			}
		}
	}

	const int num_quads = quad_vert.size()/NUM_CUBE_FACET_VERT;
	for (int i = 0; i < num_quads; i++)
	{
		for (int j = 0; j < VERT_PER_QUAD; j++)
		{
			VERTEX_INDEX v = quad_vert[NUM_CUBE_FACET_VERT*i +j];
			if (iso_vlist[v].flag_restrictionC)
			{
				iso_vlist[v].restricted_facets = 255;
			}
		}
	}
}


// set the restrictions

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
	vector<VERTEX_INDEX> &isolatedList)
{
	if (!dualiso_data.flag_no_restriction_AB)
	{
		set_restrictionsA(iso_vlist, qdual_table);
		if (!dualiso_data.flag_no_restriction_B)
		{
			set_restrictionsB(scalar_grid, vertex_coord, dualiso_data.qdual_epsilon, isovalue, iso_vlist, isodual_table,
				first_isov, qdual_table, dualiso_data.flag_collapse_debug,
				dualiso_data.flag_move_vertices, dualiso_data.flag_move_vertices2, dualiso_info);
		}
	}
	if(!dualiso_data.flag_no_restriction_C)
	{
		set_restrictionsC(scalar_grid, dualiso_data.qdual_epsilon, 
			dualiso_data.flag_move_vertices, dualiso_data.flag_move_vertices2, isovalue, quad_vert,
			iso_vlist, isodual_table, first_isov,qdual_table, 
			vertex_coord, dualiso_info, boundary_grid, isolatedList);
	}
}



