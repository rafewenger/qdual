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




//Function returns TRUE if endPt  is close to the FACET
//the facet is defined in vec_f.
bool is_epsilon_close_facet(
	const std::vector<COORD_TYPE> & vertex_coord,
	const VERTEX_INDEX &endPt, 
	const float epsilon,
	const int d,
	IJK::ARRAY<COORD_TYPE> &vcoord)
{
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

//Find the closest facet in the 'd' direction
//param 'facet_base_coord' returns the facet base coord in direction d
void find_closest_facet 
	(const int DIM3,
	//IJK::ARRAY<COORD_TYPE> &endPt1_coord,
	const COORD_TYPE * endPt1_coord,
	const int d,
	IJK::ARRAY<GRID_COORD_TYPE> &facet_base_coord)
{

	int d1 = (d+1)%3;
	int d2 = (d+2)%3;
	facet_base_coord[d] =  floor(endPt1_coord[d] + 0.5); // temporary round function.
	facet_base_coord[d1] = floorf(endPt1_coord[d1]);
	facet_base_coord[d2] = floorf(endPt1_coord[d2]);
	// DEBUG
	//cout <<"d "<< d <<" facet base "<<facet_base_coord[0]<<","<<facet_base_coord[1]<<","<<facet_base_coord[2]<<" ";
}




// Compute distance to a particular facet
void distance_to_facet (
	const int DIM3,
	//IJK::ARRAY<COORD_TYPE> &endPt_coord,
	const COORD_TYPE * endPt_coord,
	const int d,
	float &dist
	)
{
	IJK::ARRAY<GRID_COORD_TYPE> facet_base_coord(DIM3,0); //base coord of the facet in dir d
	find_closest_facet(DIM3, endPt_coord, d, facet_base_coord);
	dist = std::abs(endPt_coord[d]-facet_base_coord[d]);
}

//Find the closest edge in the "d" direction
//And the distance to the closest edge
void find_closest_edge_and_distance 
	(const int DIM3,
	const COORD_TYPE * endPt1_coord,
	const int d,
	IJK::ARRAY<GRID_COORD_TYPE> &edge_base_coord,
	float & close_dist)
{
	int d1 = (d+1)%3;
	int d2 = (d+2)%3;
	edge_base_coord[d] =  floor(endPt1_coord[d]); 
	edge_base_coord[d1] = floor(endPt1_coord[d1] + 0.5);
	edge_base_coord[d2] = floor(endPt1_coord[d2] + 0.5);

	float dist1=0, dist2=0;
	//distance_to_facet ( DIM3, endPt1_coord, d1, dist1);
	//distance_to_facet ( DIM3, endPt1_coord, d2, dist2);

	dist1 = abs(edge_base_coord[d1]-endPt1_coord[d1]);
	dist2 = abs(edge_base_coord[d2]-endPt1_coord[d2]);
	//return the larger distance
	if (dist1 < dist2)
		close_dist = dist2;
	else
		close_dist = dist1;

	// DEBUG
	//cout <<"d "<< d <<" edge base "<< edge_base_coord[0] << "," << edge_base_coord[1] 
	//<<","<<edge_base_coord[2]<<" ";
	//cout <<" endpt [" << endPt1_coord[0] <<","<< endPt1_coord[1]<<","<<endPt1_coord[2]<<"]"<<endl;
	//cout <<"dir "<<d1<<" dist "<< dist1 << " dir "<<d2<<" dist "<< dist2 <<endl;
}

//Check if the point is 'epsilon' close to facet 'f' in direction 'd'
//The facet 'f' is defined by the params 'facet_base_coord' and 'd/closest_facet'
bool is_epsilon_close(
	const COORD_TYPE * endPt1_coord,
	const IJK::ARRAY<GRID_COORD_TYPE> &facet_base_coord,
	const int d,
	const float epsilon)
{
	//cout << " dist ["<<std::abs(endPt1_coord[d] - facet_base_coord[d])<<"]" ;
	if (std::abs(endPt1_coord[d] - facet_base_coord[d]) < epsilon )
	{
		//cout <<"d "<<d<<" epsilon close\n";
		//cout <<" epsilon close endPt1[d] - facet_base[d] d= "<<d<<"  val = "<< std::abs(endPt1_coord[d] - facet_base_coord[d]) <<endl;
		return true;
	}
	else
	{
		//cout <<"\n";
		return false;
	}
}


//Setup the COLLAPSE_MAP
void setup_collapse_edges(
	const std::vector<COORD_TYPE> & vertex_coord,
	IJK::ARRAY<VERTEX_INDEX> &collapse_map)
{
	int num_vertex = vertex_coord.size();
	for (int v=0; v<num_vertex; v++)
	{
		collapse_map[v]=v;
	}
}

//Collapse edge, endPt2 is mapped to endPt1.
void update_collapse_edges(
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	const VERTEX_INDEX endPt1,
	const VERTEX_INDEX endPt2)
{
	collapse_map[endPt2]=endPt1;
}

//implements find part of the union find operation
int find_vertex_recursive(IJK::ARRAY<VERTEX_INDEX> &collapse_map, int endpt)
{
	if (collapse_map[endpt]!= endpt)
	{
		// DEBUG
		//cout <<"ce["<<endpt<<"]="<<collapse_map[endpt]<<"!="<<endpt<< endl;
		collapse_map[endpt]=find_vertex_recursive(collapse_map, collapse_map[endpt]);
	}
	return collapse_map[endpt];
}


//Implements find vertex iteratively
int find_vertex(IJK::ARRAY<VERTEX_INDEX> &collapse_map, int endpt)
{
	VERTEX_INDEX temp = endpt;
	while (collapse_map[temp]!=temp)
	{
		int temp2 = collapse_map[temp];
		collapse_map[temp]=collapse_map[temp2];
		temp = temp2;
	}
	return temp;
}



//Collapse across facets
void collapse_across_facets(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const std::vector<DUAL_ISOVERT> & iso_vlist,
	std::vector<VERTEX_INDEX> & quad_vert,
	const std::vector<COORD_TYPE> & vertex_coord,
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	const float epsilon)
{
	const int num_quads = quad_vert.size()/4;
	const int dimension = scalar_grid.Dimension();
	const int DIM3 = 3;

	int count = 4;
	int v1,v2;
	for (int q=0; q<num_quads; q++)
	{
		//for each edge
		for ( v1 = count-1, v2 = 0; v2 < count; v1 = v2++ ) {


			int endPt1 = quad_vert[q*4+v1];
			int endPt2 = quad_vert[q*4+v2];

			endPt1 = find_vertex(collapse_map, endPt1);
			endPt2 = find_vertex(collapse_map, endPt2);
			if (endPt1!=endPt2)
			{
				if ( vert_simple(iso_vlist[endPt2].patch_index ) 
					&& !vert_simple( iso_vlist[quad_vert[endPt1]].patch_index ))
				{
					swap(endPt1, endPt2);
				}
				else if (iso_vlist[quad_vert[endPt2]].cube_index < iso_vlist[quad_vert[endPt1]].cube_index)
				{
					swap(endPt1, endPt2);
				}


				const COORD_TYPE * endPt1_coord = & (vertex_coord[DIM3*endPt1]);
				const COORD_TYPE * endPt2_coord = & (vertex_coord[DIM3*endPt2]);

				int num_close=0; //number of close facets.
				int closest_facet=0;// direction of the closest facet
				IJK::ARRAY<GRID_COORD_TYPE> facet_base_coord(dimension,0); //base coord of the facet in dir d
				IJK::ARRAY<GRID_COORD_TYPE> closest_facet_base_coord(dimension,0);

				//check which facets is endPt1 closest to.
				for (int d=0; d<dimension;d++)
				{	
					find_closest_facet(DIM3, endPt1_coord, d, facet_base_coord);
					if ( is_epsilon_close(endPt1_coord, facet_base_coord, d, epsilon))
					{
						num_close++;
						closest_facet = d;
						closest_facet_base_coord[0] = facet_base_coord[0];
						closest_facet_base_coord[1] = facet_base_coord[1];
						closest_facet_base_coord[2] = facet_base_coord[2];
					}
				}

				//cout <<" num close " << num_close <<endl;
				if(num_close == 1)
				{
					//find the closest facet to endpt2 in direction closest_facet.
					find_closest_facet(DIM3, endPt2_coord, closest_facet, facet_base_coord);
					if ( is_epsilon_close(endPt2_coord, facet_base_coord, closest_facet, epsilon))
					{
						if (scalar_grid.ComputeVertexIndex(&facet_base_coord[0])
							== scalar_grid.ComputeVertexIndex(&closest_facet_base_coord[0]))
						{
							if (is_permitted_collapse)
							{ 
								cout <<"**** collapse **** "<< endPt2 <<" to " << endPt1 <<endl;
								update_collapse_edges(collapse_map,endPt1, endPt2);
							}
						}
					}
				}
				else 
				{
					
					//check which facets is endPt2 closest to.
					num_close=0;
					for (int d=0; d<dimension;d++)
					{	
						find_closest_facet(DIM3, endPt2_coord, d, facet_base_coord);
						if ( is_epsilon_close(endPt2_coord, facet_base_coord, d, epsilon))
						{
							num_close++;
							closest_facet = d;
							closest_facet_base_coord[0] = facet_base_coord[0];
							closest_facet_base_coord[1] = facet_base_coord[1];
							closest_facet_base_coord[2] = facet_base_coord[2];
						}
					}
					//cout <<"num close " << num_close <<endl;
					if (num_close==1)//endPt2 is close to a facet but not to an edge.
					{
						find_closest_facet(DIM3, endPt1_coord, closest_facet, facet_base_coord);	
						if ( is_epsilon_close(endPt1_coord, facet_base_coord, closest_facet, epsilon))
						{
							if (scalar_grid.ComputeVertexIndex(&facet_base_coord[0])
								== scalar_grid.ComputeVertexIndex(&closest_facet_base_coord[0]))
							{
								if (is_permitted_collapse)
								{ 
									cout <<"veretx " <<scalar_grid.ComputeVertexIndex(&facet_base_coord[0])
										<<" = "<<scalar_grid.ComputeVertexIndex(&closest_facet_base_coord[0])<<endl;
									cout <<"**** collapse **** "<< endPt1 <<" to " << endPt2 <<endl;
									update_collapse_edges(collapse_map,endPt2, endPt1);

								}
							}
						}
					}
				}
			}//if, only if the endPts are not already contracted		
		} // for each edge end 
	}
}



//Collapse across edges
void collapse_across_edges(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const std::vector<DUAL_ISOVERT> & iso_vlist,
	std::vector<VERTEX_INDEX> & quad_vert,
	const std::vector<COORD_TYPE> & vertex_coord,
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	const float epsilon)
{
	const int num_quads = quad_vert.size()/4;
	const int dimension = scalar_grid.Dimension();
	const int DIM3 = 3;

	int count = 4;
	int v1,v2;
	for (int q=0; q<num_quads; q++)
	{
		//for each edge
		for ( v1 = count-1, v2 = 0; v2 < count; v1 = v2++ ) {


			int endPt1 = quad_vert[q*4+v1];
			int endPt2 = quad_vert[q*4+v2];

			endPt1 = find_vertex(collapse_map, endPt1);
			endPt2 = find_vertex(collapse_map, endPt2);

			if (endPt1!=endPt2)
			{
				if ( vert_simple(iso_vlist[endPt2].patch_index ) 
					&& !vert_simple( iso_vlist[quad_vert[endPt1]].patch_index ))
				{
					swap(endPt1, endPt2);
				}
				else if (iso_vlist[quad_vert[endPt2]].cube_index < iso_vlist[quad_vert[endPt1]].cube_index)
				{
					swap(endPt1, endPt2);
				}
				
				const COORD_TYPE * endPt1_coord = & (vertex_coord[DIM3*endPt1]);
				const COORD_TYPE * endPt2_coord = & (vertex_coord[DIM3*endPt2]);

				int num_close=0; //number of close edge.
				int closest_edge_dir=0;// direction of the closest edge
				IJK::ARRAY<GRID_COORD_TYPE> edge_base_coord(dimension,0); //base coord of the edge in dir d
				IJK::ARRAY<GRID_COORD_TYPE> closest_edge_base_coord(dimension,0);

				float closest_distance=0;
				//check which edge is endPt1 closest to.
				for (int d=0; d<dimension;d++)
				{	
					find_closest_edge_and_distance(DIM3, endPt1_coord, d, edge_base_coord, closest_distance);
					if (closest_distance < epsilon)
					{
						num_close++;
						closest_edge_dir=d;
						closest_edge_base_coord[0]=edge_base_coord[0];
						closest_edge_base_coord[1]=edge_base_coord[1];
						closest_edge_base_coord[2]=edge_base_coord[2];
					}
				}

				if(num_close == 1)
				{
					//find the closest edge and distance for the endpt2
					find_closest_edge_and_distance(DIM3, endPt2_coord, closest_edge_dir,
						edge_base_coord, closest_distance);
					if (closest_distance < epsilon)
					{
						if (scalar_grid.ComputeVertexIndex(&edge_base_coord[0])
							== scalar_grid.ComputeVertexIndex(&closest_edge_base_coord[0]))
							if (is_permitted_collapse)
							{ 
								cout <<"**** collapse **** "<< endPt2 <<" to " << endPt1 <<endl;
								update_collapse_edges(collapse_map,endPt1, endPt2);
							}
					}
				}
				else
				{
					
					//edge is close to a vertex 
					num_close=0;
					for (int d=0; d<dimension;d++)
					{	
						find_closest_edge_and_distance(DIM3, endPt2_coord, d, edge_base_coord, closest_distance);
						if (closest_distance < epsilon)
						{		
							num_close++;
							closest_edge_dir=d;
							closest_edge_base_coord[0]=edge_base_coord[0];
							closest_edge_base_coord[1]=edge_base_coord[1];
							closest_edge_base_coord[2]=edge_base_coord[2];
						}
					}
					if(num_close == 1)
					{
						//find the closest edge and distance for the endpt2
						find_closest_edge_and_distance(DIM3, endPt1_coord, closest_edge_dir,
							edge_base_coord, closest_distance);
						if (closest_distance < epsilon)
						{
							if (scalar_grid.ComputeVertexIndex(&edge_base_coord[0])
								== scalar_grid.ComputeVertexIndex(&closest_edge_base_coord[0]))
								if (is_permitted_collapse)
								{ 
									cout <<"**** collapse **** "<< endPt1 << " to " << endPt2 <<endl;
									update_collapse_edges(collapse_map,endPt2, endPt1);
								}
						}
					}

				}
			} //if, only if the endPts are not already contracted		
		} // for each edge end 
	}
}



//Find the closest grid vertex to endpt_coord
void find_closest_grid_vertex(const int DIM3,
							  const COORD_TYPE * endPt_coord,
							  IJK::ARRAY<GRID_COORD_TYPE> &vertex_base_coord,
							  float & close_dist)
{
	std::vector<float> dists;
	float dist=0;
	for (int d=0;d<DIM3;d++)
	{
		vertex_base_coord[d]= floorf(endPt_coord[d]+0.5);
		dist = abs(vertex_base_coord[d]-endPt_coord[d]);
		dists.push_back(dist);
	}
	close_dist = *std::max_element(dists.begin(), dists.end());
}


//collapse across vertices
void collapse_across_vertices(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const std::vector<DUAL_ISOVERT> & iso_vlist,
	std::vector<VERTEX_INDEX> & quad_vert,
	const std::vector<COORD_TYPE> & vertex_coord,
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	const float epsilon)
{
	const int num_quads = quad_vert.size()/4;
	const int dimension = scalar_grid.Dimension();
	const int DIM3 = 3;
	int count = 4;
	int v1, v2;
	for (int q=0; q<num_quads; q++)
		//for each edge
			for ( v1 = count-1, v2 = 0; v2 < count; v1 = v2++ ) {

				int endPt1 = quad_vert[q*4+v1];
				int endPt2 = quad_vert[q*4+v2];
				endPt1 = find_vertex(collapse_map, endPt1);
				endPt2 = find_vertex(collapse_map, endPt2);
				if (endPt1!=endPt2)
				{
					if ( vert_simple(iso_vlist[endPt2].patch_index ) 
						&& !vert_simple( iso_vlist[quad_vert[endPt1]].patch_index ))
					{
						swap(endPt1, endPt2);
					}
					else if (iso_vlist[quad_vert[endPt2]].cube_index < iso_vlist[quad_vert[endPt1]].cube_index)
					{
						swap(endPt1, endPt2);
					}
					const COORD_TYPE * endPt1_coord =  & (vertex_coord[DIM3*endPt1]);
					const COORD_TYPE * endPt2_coord = & (vertex_coord[DIM3*endPt2]);

					float closest_distance=0;
					IJK::ARRAY<GRID_COORD_TYPE> vertex_base_coord(dimension,0);
					IJK::ARRAY<GRID_COORD_TYPE> closest_vertex_base_coord(dimension,0);
					find_closest_grid_vertex(DIM3, endPt1_coord, vertex_base_coord, closest_distance);
					if (closest_distance < epsilon)
					{
						for (int d=0;d<DIM3;d++)
							closest_vertex_base_coord[d]=vertex_base_coord[d];

						find_closest_grid_vertex(DIM3, endPt2_coord, vertex_base_coord, closest_distance);
						if (closest_distance < epsilon)
						{
							if (scalar_grid.ComputeVertexIndex(&vertex_base_coord[0])
								== scalar_grid.ComputeVertexIndex(&closest_vertex_base_coord[0]))
								if (is_permitted_collapse)
								{ 
									//cout <<"**** collapse **** "<<endl;
									update_collapse_edges(collapse_map,endPt1, endPt2);
								}
						}
					}
				}
			}//for each edge
}

//update the collapses in "collapse_map" in "quad_vert"
void update_quads(
	const IJK::ARRAY<VERTEX_INDEX> & collapse_map,
	std::vector<VERTEX_INDEX> & quad_vert
	)
{
	for (int v=0;v<quad_vert.size();v++)
	{
		int k = quad_vert[v];
		quad_vert[v]=collapse_map[k];
	}
}

// dual collapse main function
// NOTE:quads are reordered at the start and again at the end.
void QCOLLAPSE::dual_collapse(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	std::vector<VERTEX_INDEX> & quad_vert,
	const std::vector<DUAL_ISOVERT> & iso_vlist, 
	const std::vector<COORD_TYPE> & vertex_coord,
	const float epsilon
	)
{

	//set up the vertex collapse map
	int num_vertex = vertex_coord.size();
	IJK::ARRAY<VERTEX_INDEX> collapse_map(num_vertex,0);
	//COLLAPSE_MAP collapse_map;
	//setup collapse_map.
	setup_collapse_edges(vertex_coord, collapse_map);
	//Reordering QuadVert 
	IJK::reorder_quad_vertices(quad_vert);
	collapse_across_facets(scalar_grid, iso_vlist, quad_vert, vertex_coord, collapse_map, epsilon);
	collapse_across_edges(scalar_grid, iso_vlist, quad_vert, vertex_coord, collapse_map, epsilon);
	collapse_across_vertices(scalar_grid, iso_vlist, quad_vert, vertex_coord, collapse_map, epsilon);
	update_quads(collapse_map, quad_vert);
	//reorder the quads back to the original.
	IJK::reorder_quad_vertices(quad_vert);
}



