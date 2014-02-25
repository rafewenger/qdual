#include "qdual_restrictions.h"
#include "ijkscalar_grid.txx"


using namespace std;
using namespace IJK;
using namespace NamedConstants;




void check_edge_has_square_isosurface_path(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const DUALISO_INDEX_GRID & first_isov,
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	const IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	const QDUAL_TABLE & qdual_table,
	const int c0,
	const int d1, 
	const int d2, 
	int & num_connect
	)
{
	int d1coef[NUM_CUBE_FACET_VERT]={0,1,0,1};
	int d2coef[NUM_CUBE_FACET_VERT]= {0,0,1,1};

	for (int i = 0; i<NUM_CUBE_FACET_VERT; i++)
	{
		VERTEX_INDEX c = c0 + d1coef[i] * scalar_grid.AxisIncrement(d1) +
			d2coef[i]*scalar_grid.AxisIncrement(d2);
		VERTEX_INDEX indx_iso_vlist = first_isov.Scalar(c);
		bool flag_connect = false;

		if (indx_iso_vlist != NO_ISOVERTEX_IN_CUBE)
		{
			GRID_COORD_TYPE * c_coord = new GRID_COORD_TYPE[3];
			scalar_grid.ComputeCoord(c,c_coord);

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
				//if (((edge_flag & (1<<f1)) !=0 ) && (edge_flag & (1<<f2) != 0))
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
	}
}

// Square isosurface paths
void compute_restrictions_BList(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	vector< pair<VERTEX_INDEX, int> > &restricted_edges, 
	const bool print_info)
{


	// Restriction B list is the list of positive and negative interior grid edges
	// surrounded by a square isourface path
	VERTEX_INDEX end0=0;
	int edge_dir;
	IJK_FOR_EACH_INTERIOR_GRID_EDGE(end0, edge_dir, scalar_grid, VERTEX_INDEX)
	{
		VERTEX_INDEX end1 = scalar_grid.NextVertex(end0,edge_dir);
	

		SCALAR_TYPE s_end0, s_end1;
		s_end0 = scalar_grid.Scalar(end0);
		s_end1 = scalar_grid.Scalar(end1);

		if (((s_end0 > isovalue) && ( s_end1 < isovalue)) ||
			((s_end0 < isovalue) && ( s_end1 > isovalue)))
			continue;

		int num_connect = 0;
		int d1 = (edge_dir+1)% DIM3;
		int d2 = (edge_dir+2)%DIM3;
		// c0 is the index of the cube left and below to the edge
		VERTEX_INDEX c0 = end0 - scalar_grid.AxisIncrement(d1) - scalar_grid.AxisIncrement(d2);
		VERTEX_INDEX index_iso_vlist_c0 = first_isov.Scalar(c0);
		if (index_iso_vlist_c0 == NO_ISOVERTEX_IN_CUBE )
			continue;
		else
		{
			check_edge_has_square_isosurface_path(scalar_grid, first_isov,
				iso_vlist, isodual_table, qdual_table, c0, d1, d2, num_connect);
			if (num_connect >=3)
			{
				restricted_edges.push_back(make_pair(end0, edge_dir));
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
	QDUAL_TABLE & qdual_table
	)
{
	QDUAL_TABLE::DIR_BITS edgeFlag=0;
	int count =0;
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
				//cout <<" isovlist index "<< indx_iso_vlist+k << " cube index= "
				//	<< iso_vlist[indx_iso_vlist+k].cube_index << endl;
				if ( iso_vlist[indx_iso_vlist+k].sep_vert == iv)
				{
					count++;
					IJKDUALTABLE::TABLE_INDEX it = iso_vlist[indx_iso_vlist+k].table_index;
					unsigned char ip = iso_vlist[indx_iso_vlist+k].patch_index;
					QDUAL_TABLE::DIR_BITS w_edge_flag = qdual_table.connectDir(it,ip);
					edgeFlag = edgeFlag | w_edge_flag;

					COORD_TYPE *a = new COORD_TYPE[3];
					scalar_grid.ComputeCoord(iso_vlist[indx_iso_vlist+k].cube_index ,a);
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

// compute the vertex list which have isosurface boxes around it 
void compute_restrictions_CList(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	const std::vector<COORD_TYPE> & vertex_coord,
	vector<VERTEX_INDEX> &restriction_Clist)
{
	IJK::BOOL_GRID<DUALISO_GRID> boundary_grid;
	boundary_grid.SetSize(scalar_grid);
	boundary_grid.SetAll(false);
	compute_boundary_grid(boundary_grid);
	int numv =  scalar_grid.NumVertices();

	for (VERTEX_INDEX iv = 0; iv < numv; iv++)
	{
		if (!boundary_grid.Scalar(iv)) 
		{
			VERTEX_INDEX indx_iso_vlist_iv = first_isov.Scalar(iv);
			if (indx_iso_vlist_iv ==-1)
				continue;
			else
			{
				bool flag_box = check_vertex_has_box_around 
					(scalar_grid, iv, iso_vlist, isodual_table,
					first_isov, qdual_table);

				if (flag_box)
					restriction_Clist.push_back(iv);
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
	for (int i = 0; i < iso_vlist.size(); i++)
	{
		it = iso_vlist[i].table_index;
		ip = iso_vlist[i].patch_index;
		VERTEX_INDEX icube = iso_vlist[i].cube_index;
		VERTEX_INDEX icorner = qdual_table.V1(it,ip);
		if (icorner !=255 )
		{
			iso_vlist[i].sep_vert = scalar_grid.CubeVertex(icube, icorner);
			scalar_grid.ComputeCoord(i, c);

			//cout <<"debug vertex "<< i 
			//	<< " ["<< c[0]<<","<<c[1]<<","<<c[2]<<"]";

			scalar_grid.ComputeCoord(iso_vlist[i].sep_vert , c);

			//cout <<" separates " << iso_vlist[i].sep_vert <<" [ " 
			//	<< " ["<< c[0]<<","<<c[1]<<","<<c[2]<<"]"<<endl;
		}
		else
		{
			/*iso_vlist[i].sep_vert = scalar_grid.NumVertices();
			cout <<"debug vertex "<< i << " seperates NOTHING " <<  iso_vlist[i].sep_vert << endl;  */
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

void set_restrictionsB(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	bool print_info,
	std::unique_ptr<RestrictionInfo> &rs_info)
{
	
	vector< pair<VERTEX_INDEX, int> > restricted_edges;

	compute_restrictions_BList( scalar_grid, isovalue, iso_vlist, 
		isodual_table, first_isov, qdual_table, restricted_edges, print_info);
	//store information 
	rs_info->restriction_BList_size = restricted_edges.size();
	rs_info->restricted_edges_info = restricted_edges;
	
	int d1coef[NUM_CUBE_FACET_VERT]={0,1,0,1};
	int d2coef[NUM_CUBE_FACET_VERT]= {0,0,1,1};

	for (int i=0;i < restricted_edges.size(); i++)
	{
		VERTEX_INDEX end0 = restricted_edges[i].first;
		GRID_COORD_TYPE * coord = new GRID_COORD_TYPE[3];
		scalar_grid.ComputeCoord(end0,coord);

		int edge_dir = restricted_edges[i].second;
		int d1 = (edge_dir + 1) % DIM3;
		int d2 = (edge_dir + 2) % DIM3;
		VERTEX_INDEX c0 = end0 - scalar_grid.AxisIncrement(d1) - scalar_grid.AxisIncrement(d2);


		for (int j=0; j < NUM_CUBE_FACET_VERT; j++)
		{
			VERTEX_INDEX c = c0 + d1coef[j]*scalar_grid.AxisIncrement(d1) 
				+ d2coef[j]*scalar_grid.AxisIncrement(d2);
			scalar_grid.ComputeCoord(c,coord);
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
						unsigned char test_bef = iso_vlist[indx_iso_vlist+k].restricted_facets;
						iso_vlist[indx_iso_vlist+k].restricted_facets = 
							(iso_vlist[indx_iso_vlist+k].restricted_facets | ((1<<f1) | (1<<f2)));
						unsigned char test_after = iso_vlist[indx_iso_vlist+k].restricted_facets;
					}

				}

			}
		}
	}
}



void set_restrictionsC(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	const std::vector<COORD_TYPE> & vertex_coord,
	std::unique_ptr<RestrictionInfo> &rs_info)
{
	vector<VERTEX_INDEX> restriction_Clist;
	compute_restrictions_CList( scalar_grid, isovalue, iso_vlist, isodual_table, first_isov,
		qdual_table, vertex_coord, restriction_Clist);
	//store info
	
	rs_info->restriction_CList_size = restriction_Clist.size();
	rs_info->restricted_vertex_info = restriction_Clist;

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
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	const std::vector<COORD_TYPE> & vertex_coord,
	std::unique_ptr<RestrictionInfo> &rs_info)
{
	if (!dualiso_data.flag_no_restriction_AB)
	{
		set_restrictionsA(iso_vlist, qdual_table);
		if (!dualiso_data.flag_no_restriciton_B)
		{
			set_restrictionsB(scalar_grid, isovalue, iso_vlist, isodual_table,
				first_isov, qdual_table, dualiso_data.flag_collapse_debug, rs_info);
		}
	}
	if(!dualiso_data.flag_no_restriciton_C)
	{
		set_restrictionsC(scalar_grid, isovalue, quad_vert,
			iso_vlist, isodual_table, first_isov,qdual_table, vertex_coord, rs_info);
	}
}