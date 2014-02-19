#include "qdual_restrictions.h"
#include "ijkscalar_grid.txx"

using namespace std;
using namespace IJK;
using namespace NamedConstants;

// Square isosurface paths
void compute_restrictions_BList(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	DUALISO_INDEX_GRID & first_isov,
	QDUAL_TABLE & qdual_table,
	vector< pair<VERTEX_INDEX, int> > &restricted_edges)
{
	int d1coef[NUM_CUBE_FACET_VERT]={0,1,0,1};
	int d2coef[NUM_CUBE_FACET_VERT]= {0,0,1,1};

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

		int num_connect =0;
		int d1 = (edge_dir+1)% DIM3;
		int d2 = (edge_dir+2)%DIM3;
		// c0 is th eindec of the cube left and below
		VERTEX_INDEX c0 = end0 - scalar_grid.AxisIncrement(d1) - scalar_grid.AxisIncrement(d2);

		for (int i = 0; i < NUM_CUBE_FACET_VERT; i++)
		{
			VERTEX_INDEX c = c0 + d1coef[i] * scalar_grid.AxisIncrement(d1) +
				d2coef[i]*scalar_grid.AxisIncrement(d2);

			VERTEX_INDEX f1 = ((d1+ (1-d1coef[i]))* DIM3) % NUM_CUBE_FACETS;
			VERTEX_INDEX f2 = ((d2+ (1-d2coef[i]))* DIM3) % NUM_CUBE_FACETS;
			bool flag_connect = false;
			VERTEX_INDEX indx_iso_vlist = first_isov.Scalar(c);

			if (indx_iso_vlist != -1){

				int table_ind = iso_vlist[indx_iso_vlist].table_index;
				int num_iso_verts = isodual_table.NumIsoVertices(table_ind);
				//for each isosurface vertex in cube c
				for (int j = 0; j < num_iso_verts; j++)
				{
					IJKDUALTABLE::TABLE_INDEX it = iso_vlist[indx_iso_vlist+j].table_index;
					FACET_VERTEX_INDEX ip = iso_vlist[indx_iso_vlist+j].patch_index;
					QDUAL_TABLE::DIR_BITS edge_flag = qdual_table.connectDir(it, ip);
					if (((edge_flag & (1<<f1)) !=0 ) && (edge_flag & (1<<f2) != 0))
					{
						flag_connect = true;
					}
				}
			}
			if (flag_connect) 
			{
				num_connect++;
			}
		}
		if (num_connect >=3)
		{
			//cout << end0 <<" - " << scalar_grid.NextVertex(end0,edge_dir) <<" is restricted "<<endl;
			//COORD_TYPE *a = new COORD_TYPE[3];
			//scalar_grid.ComputeCoord(end0,a);
			//cout <<" -- " << a[0]<<","<<a[1]<<","<<a[2];
			//scalar_grid.ComputeCoord(scalar_grid.NextVertex(end0,edge_dir),a);
			//cout <<" -- " << a[0]<<","<<a[1]<<","<<a[2] << endl;
			restricted_edges.push_back(make_pair(end0, edge_dir));
		}
	}
}

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
	for (VERTEX_INDEX iv = 0; iv < numv; iv++) {
		if (!boundary_grid.Scalar(iv)) {
			QDUAL_TABLE::DIR_BITS edgeFlag=0;
			int count =0;
			VERTEX_INDEX iv0 = iv - scalar_grid.CubeVertexIncrement(NUM_CUBE_VERT-1);
			for (int j = 0; j < NUM_CUBE_VERT; j++)
			{
				VERTEX_INDEX icube = scalar_grid.CubeVertex(iv0,j);
				VERTEX_INDEX indx_iso_vlist = first_isov.Scalar(icube);
				if (indx_iso_vlist != -1)
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
		}
		else
		{
			iso_vlist[i].sep_vert = scalar_grid.NumVertices();
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
	QDUAL_TABLE & qdual_table)
{
	
	vector< pair<VERTEX_INDEX, int> > restricted_edges;
	compute_restrictions_BList( scalar_grid, isovalue, iso_vlist, 
		isodual_table, first_isov, qdual_table, restricted_edges);

	int d1coef[NUM_CUBE_FACET_VERT]={0,1,0,1};
	int d2coef[NUM_CUBE_FACET_VERT]= {0,0,1,1};

	for (int i=0;i < restricted_edges.size(); i++)
	{
		VERTEX_INDEX end0 = restricted_edges[i].first;
		int edge_dir = restricted_edges[i].second;
		int d1 = (edge_dir + 1) % DIM3;
		int d2 = (edge_dir + 2) % DIM3;
		VERTEX_INDEX c0 = end0 - scalar_grid.AxisIncrement(d1) - scalar_grid.AxisIncrement(d2);


		for (int j=0; j < NUM_CUBE_FACET_VERT-1; j++)
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
					if (     (edge_flag & (1<<f1) != 0) && ((edge_flag & (1<f2))!=0) )
					{
						iso_vlist[indx_iso_vlist+k].restricted_facets = 
							(iso_vlist[indx_iso_vlist+k].restricted_facets | (1<<f1) | (1<<f2));
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
	const std::vector<COORD_TYPE> & vertex_coord)
{
	vector<VERTEX_INDEX> restriction_Clist;
	compute_restrictions_CList( scalar_grid, isovalue, iso_vlist, isodual_table, first_isov,
		qdual_table, vertex_coord, restriction_Clist);

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
					//cout <<" isovlist index "<< indx_iso_vlist+k << " cube index= "
					//	<< iso_vlist[indx_iso_vlist+k].cube_index << endl;
					if ( iso_vlist[indx_iso_vlist+k].sep_vert == v)
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
		for (int j = 0; j < NUM_CUBE_FACET_VERT; j++)
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
	const std::vector<COORD_TYPE> & vertex_coord)
{
	if (!dualiso_data.flag_no_restriction_AB)
	{
		set_restrictionsA(iso_vlist, qdual_table);
		if (!dualiso_data.flag_no_restriciton_B)
		{
			set_restrictionsB(scalar_grid, isovalue, iso_vlist, isodual_table,
				first_isov, qdual_table);
		}
	}
	if(!dualiso_data.flag_no_restriciton_C)
	{
		set_restrictionsC(scalar_grid, isovalue, quad_vert,
			iso_vlist, isodual_table, first_isov,qdual_table, vertex_coord);
	}
}