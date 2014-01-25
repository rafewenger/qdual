#include <cmath> 
#include <algorithm>
#include "qdualCollapse.h"

using namespace QCOLLAPSE;
using namespace std;


//Check if patch is simple 
//Return True if the vert is simple
inline bool vert_simple(const unsigned char patch_ind){
	if ((int)patch_ind == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//For all directions
void epsilon_close_facet(
	const std::vector<COORD_TYPE> & vertex_coord,
	const VERTEX_INDEX &endPt1,
	const VERTEX_INDEX &v, // vertex index defining the facets 
	const float epsilon,
	IJK::ARRAY<COORD_TYPE> &vcoord, 
	vector<FACET>  &vec_f)
{	
	for (int d=0; d<3; d++)
	{
		if ( std::abs(vertex_coord[3*endPt1+d] - vcoord[d] ) < epsilon )
		{
			FACET f;
			f = make_pair (v,d);
			vec_f.push_back(f);
		}
	}
}

//For particular direction.
void epsilon_close_facet(
	const std::vector<COORD_TYPE> & vertex_coord,
	const VERTEX_INDEX &endPt,
	const VERTEX_INDEX &v, // vertex index defining the facets 
	const float epsilon,
	IJK::ARRAY<COORD_TYPE> &vcoord,
	const int d, // the facet being checked will be orthogonal to this direction
	vector<FACET>  &vec_f)
{	
	if ( std::abs(vertex_coord[3*endPt+d] - vcoord[d] ) < epsilon )
	{
		FACET f;
		f = make_pair (v,d);
		vec_f.push_back(f);
	}
}

//Function returns TRUE if endPt  is close to the FACET
//the facet is defined in vec_f.
bool is_epsilon_close_facet(
	const std::vector<COORD_TYPE> & vertex_coord,
	const VERTEX_INDEX &endPt, 
	const float epsilon,
	const int d,
	IJK::ARRAY<COORD_TYPE> &vcoord)
{
	cout <<" base coord in vec_f (" << vcoord[0] <<" "<<vcoord[1]<<" "<< vcoord[2] <<")"<< endl; 
	cout <<" direction "<< d <<endl;
	if ( std::abs(vertex_coord[3*endPt+d] - vcoord[d] ) < epsilon )
		return true;
	else
		return false;
}


//Is the collapse permitted ?
//Returns TRUE if it is.
//***Currently always returns true.
bool is_permitted_collapse()
{
	return true;
}

//Collapse across facets
void collapse_across_facets(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const std::vector<DUAL_ISOVERT> & iso_vlist,
	std::vector<VERTEX_INDEX> & quad_vert,
	const std::vector<COORD_TYPE> & vertex_coord,
	const float epsilon)
{
	int num_quads = quad_vert.size()/4;
	int dimension = scalar_grid.Dimension();
	// DEBUG
	cout <<"dimesion " << dimension << endl;
	std::cout <<"number of quad verts "<<quad_vert.size() << " num quads " << num_quads <<endl;
	//for (int q=0; q < num_quads; q++)
	//{
	//	cout << q << " ["
	//		<< quad_vert[q*4] << " "
	//		<< quad_vert[q*4+1] << " "
	//		<< quad_vert[q*4+2] << " "
	//		<< quad_vert[q*4+3] << "]"<< endl;
	//}

	//for (int q=0; q < num_quads; q++)
	//{
	//	cout << q << " vertex index ("
	//		<< quad_vert[q*4] << ") info (cube index " << iso_vlist[quad_vert[q*4]].cube_index
	//		<< ") (Patch index " << (int)iso_vlist[quad_vert[q*4]].patch_index
	//		<< ") (Table index " << iso_vlist[quad_vert[q*4]].table_index <<")"<<endl;
	//}
	// DEBUG end

	int count = 4;
	int v1,v2;
	for (int q=0; q<10; q++)
	{

		cout <<"\nvertices are ["
			<< quad_vert[q*4] << " "
			<< quad_vert[q*4+1] << " "
			<< quad_vert[q*4+2] << " "
			<< quad_vert[q*4+3] << "]"<< endl;

		//for each edge
		for ( v1 = count-1, v2 = 0; v2 < count; v1 = v2++ ) {
			// DEBUG
			/*	cout <<"edges are "
			<< quad_vert[q*4+v1] << " "
			<< quad_vert[q*4+v2] << endl;*/

			int endPt1 = quad_vert[q*4+v1];
			int endPt2 = quad_vert[q*4+v2];
			// DEBUG 
			cout <<"\nendpoint1 = "<< endPt1 <<" endpt2 = "<< endPt2 <<endl;

			if ( vert_simple(iso_vlist[endPt2].patch_index ) 
				&& !vert_simple( iso_vlist[quad_vert[endPt1]].patch_index ))
			{
				swap(endPt1, endPt2);
				//DEBUG
				cout <<"post swap endpoint1 = "<< endPt1 <<" endpt2 = "<< endPt2 <<endl;
			}
			else if (iso_vlist[quad_vert[endPt2]].cube_index < iso_vlist[quad_vert[endPt1]].cube_index)
			{
				swap(endPt1, endPt2);
			}
			// DEBUG
			cout << vertex_coord[3*endPt1+0] <<" " << vertex_coord[3*endPt1+1] << " " << vertex_coord[3*endPt1+2] << endl;
			cout << vertex_coord[3*endPt2+0] <<" " << vertex_coord[3*endPt2+1] << " " << vertex_coord[3*endPt2+2] << endl;

			VERTEX_INDEX v = iso_vlist[endPt1].cube_index;
			IJK::ARRAY<COORD_TYPE> vcoord(dimension);

			scalar_grid.ComputeCoord( iso_vlist[endPt1].cube_index, vcoord.Ptr());
			// DEBUG
			cout <<" vcoord (" << vcoord[0] <<" "<<vcoord[1]<<" "<< vcoord[2] <<")"<< endl; 

			vector<FACET>  vec_f; // set of facets that endpt1 is close to.
			//using the all direction version 
			//epsilon_close_facet(vertex_coord, endPt1, v, epsilon, vcoord, vec_f);
			
			//the third and the sixth parameters together define a facet.
			epsilon_close_facet(vertex_coord, endPt1, v, epsilon, vcoord, 0, vec_f);
			epsilon_close_facet(vertex_coord, endPt1, v, epsilon, vcoord, 1, vec_f);
			epsilon_close_facet(vertex_coord, endPt1, v, epsilon, vcoord, 2, vec_f);

			VERTEX_INDEX v_axis; // the other 3 vertices which set up all the facets of the cube, cube_index
			for (int d=0; d<dimension; d++)
			{
				VERTEX_INDEX v_axis =  scalar_grid.NextVertex(v,d);
				scalar_grid.ComputeCoord( v_axis, vcoord.Ptr());
				epsilon_close_facet(vertex_coord, endPt1, v_axis, epsilon, vcoord, d, vec_f);
			}

						// DEBUG 
			if (vec_f.size() > 0){
				for (int i=0;i < vec_f.size(); i++)
				{
					scalar_grid.ComputeCoord( vec_f[i].first, vcoord.Ptr());
					cout << vec_f[i].first 
						<<" vcoord [" << vcoord[0] <<" "<<vcoord[1]<<" "<< vcoord[2]<<"] "
						<<" - "<< vec_f[i].second <<",";
				}
				cout <<endl;
			}
			int closeto_nm_facets=vec_f.size();
			if ( closeto_nm_facets == 1) //close to a facet but not to an edge
			{
				scalar_grid.ComputeCoord( vec_f[0].first, vcoord.Ptr());
				bool is_close = is_epsilon_close_facet(vertex_coord, endPt2, epsilon, vec_f[0].second, vcoord);
				if (is_close)
				{
					bool permitted_collapse = is_permitted_collapse();
					if (permitted_collapse)
					{
						cout <<"**** Collapse "<<endl;
					}
				}
			}

		}
	}
}

// dual collapse main function
void QCOLLAPSE::dual_collapse(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
	std::vector<VERTEX_INDEX> & quad_vert,
	const std::vector<DUAL_ISOVERT> & iso_vlist, 
	const std::vector<COORD_TYPE> & vertex_coord,
	const float epsilon
	)
{
	collapse_across_facets(scalar_grid, iso_vlist, quad_vert, vertex_coord, epsilon);
}