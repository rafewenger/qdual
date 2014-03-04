#include "qdual_collapse.h"
#include "qdual_remove_degen.h" // for the count degree function 

using namespace QCOLLAPSE;
using namespace std;
using namespace NamedConstants;

//Double comparison
template <typename TY>
bool AreSame(TY a, TY b)
{
	return abs(a - b) < 0.00001;
}

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


inline  void print_collapse_info (
	string str, 
	const int endPt1,
	const int endPt2,
	GRID_COORD_TYPE * base_coord)
{
	cout <<str<<" | "<< endPt2 <<" mapped to " << endPt1<<"; ";
	cout <<"facet base coord "<< base_coord[0]<<","
		<<base_coord[1]<<","<<base_coord[2]<<endl;
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


//Is the endpt permitted to collapse across facet, edge, vertex ?
//Returns TRUE if it is.

bool is_permitted_collapse(
	const std::vector<DUAL_ISOVERT> & iso_vlist,
	const VERTEX_INDEX v, // is the index into iso_vlist
	const GRID_COORD_TYPE * base_coord,
	const int d, // direction of the closest facet, for VERTEX this does not matter.
	COLLAPSE_TYPE c)
{
	int	closest_cube_facet_1 = 0;
	int	closest_cube_facet_2 = 0;
	switch (c){
	case FACET:
		{
			if ( AreSame<GRID_COORD_TYPE>(iso_vlist[v].cube_coord[d], base_coord[d]))
			{
				closest_cube_facet_1 = d;

			}
			else
			{
				closest_cube_facet_1 = d+DIM3;
			}

			bool cond1 = (iso_vlist[v].restricted_facets & (1<<closest_cube_facet_1));
			if (cond1)
			{
				return false;
			}
			else
			{
				return true;
			}
			break;
		}
	case EDGE:
		{
			int d1 = (d+1)%DIM3;
			int d2 = (d+2)%DIM3;
			unsigned char edge_mask = 0;
			if (AreSame<GRID_COORD_TYPE> (base_coord[d1], iso_vlist[v].cube_coord[d1]))
				edge_mask = (edge_mask | 1<< d1);
			else
				edge_mask = (edge_mask | 1<< (d1+DIM3));
			if (AreSame<GRID_COORD_TYPE> (base_coord[d2], iso_vlist[v].cube_coord[d2]))
				edge_mask = (edge_mask | 1<< d2);
			else
				edge_mask = (edge_mask | 1<< (d2+DIM3));

			if (iso_vlist[v].restricted_facets & edge_mask)
				return false;

			break;
		}
	case VERTEX:
		{
			unsigned char vertex_mask = 0;
			for (int f = 0; f < DIM3; f++)
			{
				if (AreSame<GRID_COORD_TYPE> (base_coord[f], iso_vlist[v].cube_coord[f]))
					vertex_mask = vertex_mask | (1<<f);
				else
					vertex_mask = vertex_mask | (1<<(f+DIM3));
			}
			if (iso_vlist[v].restricted_facets & vertex_mask){
				return false;
			}
			break;
		}
	}
	return true;
}

//Find the closest facet in the 'd' direction
//param 'facet_base_coord' returns the facet base coord in direction d
void find_closest_facet 
	(const int DIM3,
	const COORD_TYPE * endPt1_coord,
	const int d,
	IJK::ARRAY<GRID_COORD_TYPE> &facet_base_coord)
{
	int d1 = (d+1)%3;
	int d2 = (d+2)%3;
	facet_base_coord[d] =  (int)floor(endPt1_coord[d] + 0.5); // temporary round function.
	facet_base_coord[d1] = (int)floorf(endPt1_coord[d1]);
	facet_base_coord[d2] = (int)floorf(endPt1_coord[d2]);
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
	edge_base_coord[d] =  (int)floor(endPt1_coord[d]); 
	edge_base_coord[d1] = (int)floor(endPt1_coord[d1] + 0.5);
	edge_base_coord[d2] = (int)floor(endPt1_coord[d2] + 0.5);

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

}

//Check if the point is 'epsilon' close to facet 'f' in direction 'd'
//The facet 'f' is defined by the params 'facet_base_coord' and 'd/closest_facet'
bool is_epsilon_close(
	const COORD_TYPE * endPt1_coord,
	const IJK::ARRAY<GRID_COORD_TYPE> &facet_base_coord,
	const int d,
	const float epsilon)
{
	if (std::abs(endPt1_coord[d] - facet_base_coord[d]) < epsilon )
	{
		return true;
	}
	else
	{
		return false;
	}
}


//Setup the COLLAPSE MAP
//param 1 : collapse map
//param 2 : size of the map
void QCOLLAPSE::setup_collapse_map(
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	const int num_vertex)
{
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
		collapse_map[endpt]=find_vertex_recursive(collapse_map, collapse_map[endpt]);
	}
	return collapse_map[endpt];
}


//Implements find vertex iteratively
int QCOLLAPSE::find_vertex(IJK::ARRAY<VERTEX_INDEX> &collapse_map, int endpt)
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
	const float epsilon,
	const bool print_info,
	DUALISO_INFO & dualiso_info)
{
	const int num_quads = quad_vert.size()/4;
	const int dimension = scalar_grid.Dimension();

	COORD_TYPE * a = new COORD_TYPE[3];
	int v1,v2;
	for (int q=0; q<num_quads; q++)
	{
		//for each edge
		for ( v1 = VERT_PER_QUAD-1, v2 = 0; v2 < VERT_PER_QUAD; v1 = v2++ ) {
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


				if(num_close == 1)
				{
					//find the closest facet to endpt2 in direction closest_facet.
					find_closest_facet(DIM3, endPt2_coord, closest_facet, facet_base_coord);

					if ( is_epsilon_close(endPt2_coord, facet_base_coord, closest_facet, epsilon))
					{
						if (scalar_grid.ComputeVertexIndex(&facet_base_coord[0])
							== scalar_grid.ComputeVertexIndex(&closest_facet_base_coord[0]))
						{
							if (is_permitted_collapse (iso_vlist, endPt2,  &(facet_base_coord[0]),
								closest_facet,  FACET))
							{ 
								if (print_info)
									print_collapse_info("collapse across facets [permitted]", 
									endPt1, endPt2, &(facet_base_coord[0]));

								dualiso_info.col_info.permitted_facet_restriction++;
								update_collapse_edges(collapse_map,endPt1, endPt2);
							}
							else
							{
								dualiso_info.col_info.not_permitted_facet_restriction++;
								if (print_info)
									print_collapse_info("collapse across facets [NOT permitted]", 
									endPt1, endPt2, &(facet_base_coord[0]));
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

					if (num_close==1)//endPt2 is close to a facet but not to an edge.
					{
						find_closest_facet(DIM3, endPt1_coord, closest_facet, facet_base_coord);	
						if ( is_epsilon_close(endPt1_coord, facet_base_coord, closest_facet, epsilon))
						{
							if (scalar_grid.ComputeVertexIndex(&facet_base_coord[0])
								== scalar_grid.ComputeVertexIndex(&closest_facet_base_coord[0]))
							{
								if (is_permitted_collapse (iso_vlist, endPt1,  &(facet_base_coord[0]),
									closest_facet, FACET))
								{ 

									if (print_info)
										print_collapse_info("collapse across facets [permitted]", 
										endPt2, endPt1, &(facet_base_coord[0]));

									dualiso_info.col_info.permitted_facet_restriction++;
									update_collapse_edges(collapse_map,endPt2, endPt1);

								}
								else{
									dualiso_info.col_info.not_permitted_facet_restriction++;
									if (print_info)
										print_collapse_info("collapse across facets [NOT permitted]", 
										endPt2, endPt1, &(facet_base_coord[0]));
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
	const float epsilon,
	const bool print_info,
	DUALISO_INFO & dualiso_info)
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

					if (closest_distance < epsilon)// end pt2 is also epislon close to the edge
					{
						if (scalar_grid.ComputeVertexIndex(&edge_base_coord[0])
							== scalar_grid.ComputeVertexIndex(&closest_edge_base_coord[0]))
						{
							if (is_permitted_collapse(iso_vlist, endPt2,  &(edge_base_coord[0]),
								closest_edge_dir, EDGE))
							{ 
								if (print_info){
									print_collapse_info("collapse across edges [permitted]", 
										endPt1, endPt2, &(edge_base_coord[0]));
									cout <<"endpt2 "<< endPt2<<"("<< endPt2_coord[0]<<","<< endPt2_coord[1]<<","<< endPt2_coord[2]<<"\n";
									cout <<"endpt1 "<< endPt1<<"("<< endPt1_coord[0]<<","<< endPt1_coord[1]<<","<< endPt1_coord[2]<<"\n";
									cout <<"edge base "<< edge_base_coord[0]<<","
										<< edge_base_coord[1]<<"," << edge_base_coord[2] << endl;
									cout <<"edge dir "<< closest_edge_dir <<endl;
								}

								dualiso_info.col_info.permitted_edge_restriction++;
								update_collapse_edges(collapse_map,endPt1, endPt2);
							}
							else
							{
								dualiso_info.col_info.not_permitted_edge_restriction++;
								if (print_info){
									print_collapse_info("collapse across edges [NOT permitted]", 
										endPt1, endPt2, &(edge_base_coord[0]));
									int v=0;
									v=scalar_grid.ComputeVertexIndex(&edge_base_coord[0]);
									cout <<"edgebase vertex index "<<v<<endl;

									cout <<"endpt2 "<< endPt2<<"("<< endPt2_coord[0]<<","<< endPt2_coord[1]<<","<< endPt2_coord[2]<<"\n";
									cout <<"endpt1 "<< endPt1<<"("<< endPt1_coord[0]<<","<< endPt1_coord[1]<<","<< endPt1_coord[2]<<"\n";
									cout <<"edge base "<< edge_base_coord[0]<<","
										<< edge_base_coord[1]<<"," << edge_base_coord[2] << endl;
									cout <<"edge dir "<< closest_edge_dir <<endl;
								}

							}
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
								if (is_permitted_collapse(iso_vlist, endPt1,  &(edge_base_coord[0]),
									closest_edge_dir, EDGE))
								{ 
									if (print_info){
										print_collapse_info("collapse across edges [permitted]", 
											endPt2, endPt1, &(edge_base_coord[0]));
										int v=0;
										v=scalar_grid.ComputeVertexIndex(&edge_base_coord[0]);
										cout <<"edgebase vertex index "<<v<<endl;

										cout <<"endpt1 "<< endPt1<<"("<< endPt1_coord[0]<<","<< endPt1_coord[1]<<","<< endPt1_coord[2]<<"\n";

										cout <<"edge base "<< edge_base_coord[0]<<","
											<< edge_base_coord[1]<<"," << edge_base_coord[2] << endl;
										cout <<"edge dir "<< closest_edge_dir <<endl;
									}

									dualiso_info.col_info.permitted_edge_restriction++;
									update_collapse_edges(collapse_map,endPt2, endPt1);
								}
								else
								{
									dualiso_info.col_info.not_permitted_edge_restriction++;
									if (print_info){
										print_collapse_info("collapse across edges [ NOT permitted]", 
											endPt2, endPt1, &(edge_base_coord[0]));
										cout <<"endpt1 "<< endPt1<<"("<< endPt1_coord[0]<<","<< endPt1_coord[1]<<","<< endPt1_coord[2]<<"\n";
										cout <<"endpt2 "<< endPt2<<"("<< endPt2_coord[0]<<","<< endPt2_coord[1]<<","<< endPt2_coord[2]<<"\n";
										cout <<"edge base "<< edge_base_coord[0]<<","
											<< edge_base_coord[1]<<"," << edge_base_coord[2] << endl;
										cout <<"edge dir "<< closest_edge_dir <<endl;
									}
								}
						}
					}

				}
			} //if, only if the endPts are not already contracted		
		} // for each edge end 
	}
}



//Find the closest grid vertex to endpt_coord
void find_closest_grid_vertex(
	const int DIM3,
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
	const float epsilon,
	const bool print_info,
	DUALISO_INFO & dualiso_info)
{
	const int num_quads = quad_vert.size()/VERT_PER_QUAD;
	const int dimension = scalar_grid.Dimension();


	int v1, v2;
	for (int q=0; q<num_quads; q++)
		//for each edge
			for ( v1 = NUM_CUBE_FACET_VERT-1, v2 = 0; v2 < NUM_CUBE_FACET_VERT; v1 = v2++ ) {

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
							{	
								bool condEndPt1 = is_permitted_collapse( iso_vlist, endPt1, &(vertex_base_coord[0]),
									1, VERTEX);
								bool condEndPt2 = is_permitted_collapse( iso_vlist, endPt2, &(vertex_base_coord[0]),
									1, VERTEX);
								if ( condEndPt1 && condEndPt2 )
								{ 
									if (print_info)
										print_collapse_info("collapse across vertex [permitted]", 
										endPt1, endPt2, &(vertex_base_coord[0]));

									dualiso_info.col_info.permitted_vertex_restriction++;
									update_collapse_edges(collapse_map,endPt1, endPt2);
								}
								else
								{
									dualiso_info.col_info.not_permitted_vertex_restriction++;
									if (print_info){
										print_collapse_info("collapse across vertex [NOT permitted]", 
											endPt1, endPt2, &(vertex_base_coord[0]));

										cout <<"endPt1 condition  "<< condEndPt1 <<"  endPt2 condition "<< condEndPt2 <<endl;
									}
								}
							}
						}
					}
				}
			}//for each edge
}

void intersectQuadAndEdgeAndMoveVert
	(
	const int v,
	const int q, //qth quad
	std::vector<VERTEX_INDEX> & quad_vert,
	IJK::ARRAY<GRID_COORD_TYPE> &edgeBase,
	int & edgeDir,
	std::vector<COORD_TYPE> & vertex_coord,
	DUALISO_INFO & dualiso_info
	)
{
    const int d1 = (edgeDir+1+DIM3)%DIM3;
	const int d2 = (edgeDir+2+DIM3)%DIM3;
    vertex_coord[DIM3*v+d1]=edgeBase[d1];
	vertex_coord[DIM3*v+d2]=edgeBase[d2];
    COORD_TYPE c = 0;
	for (int i = 0; i < VERT_PER_QUAD; i++)
	{
        int qi = quad_vert[VERT_PER_QUAD*q+i];
        c=c+ vertex_coord[DIM3*qi+edgeDir];
	}
    vertex_coord[DIM3*v+edgeDir]=c/4.0;
}


// Collapse vertices of cap quad q to a vertex on edge e
void collapseToEdge(
	const int v,
	const int q, //qth quad
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	const std::vector<DUAL_ISOVERT> & iso_vlist,
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	IJK::ARRAY<GRID_COORD_TYPE> &edgeBase,
	int & edgeDir,
	DUALISO_INFO & dualiso_info
	)
{
	for (int i = 0; i < VERT_PER_QUAD; i++)
	{
		int v2 = quad_vert[VERT_PER_QUAD*q+i];
		if (iso_vlist[v2].ver_degree == 3
			&& v2 != v)
		{
			dualiso_info.col_info.permitted_edge_restriction++;
			update_collapse_edges(collapse_map, v, v2);
		}
	}
	intersectQuadAndEdgeAndMoveVert(v, q, quad_vert, edgeBase,
		edgeDir, vertex_coord, dualiso_info );
}

void capQuadComputation(
	const VERTEX_INDEX d1, 
	const VERTEX_INDEX d2, // is the diagonal vertex which determines the cube.
	IJK::ARRAY<COORD_TYPE> & d1Coord,
	IJK::ARRAY<COORD_TYPE> & d2Coord,
	const int q, //qth quad
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const float epsilon,
	const int notDegreeThree,
	const std::vector<DUAL_ISOVERT> & iso_vlist,
	IJK::ARRAY<GRID_COORD_TYPE> &edgeBase,
	int & edgeDir,
	bool &capQuad,
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	DUALISO_INFO & dualiso_info,
    bool printInfo
	)
{

	if (iso_vlist[d1].sep_vert == iso_vlist[d2].sep_vert
		&& iso_vlist[d1].ver_degree == 3 && iso_vlist[d2].ver_degree == 3)
	{ 
		//CAP QUAD 
		float closestDistance=0;
		IJK::ARRAY<GRID_COORD_TYPE> d2EdgeBase(DIM3,0);
		IJK::ARRAY<GRID_COORD_TYPE> d1EdgeBase(DIM3,0);

		find_closest_edge_and_distance(DIM3, &(d2Coord[0]), edgeDir, d2EdgeBase, closestDistance);
		if(closestDistance < epsilon)
		{
			if (scalar_grid.ComputeVertexIndex(&d2EdgeBase[0])
				== scalar_grid.ComputeVertexIndex(&edgeBase[0]))
			{
				if (is_permitted_collapse(iso_vlist, d2,  &(d2EdgeBase[0]),
					edgeDir, EDGE))
				{
					collapseToEdge(d2, q, quad_vert, vertex_coord, iso_vlist
						,collapse_map, edgeBase, edgeDir, dualiso_info);
				
					dualiso_info.cp_info.moved2Edge++;
					if ( printInfo)
					{
						cout <<"capquad "
							<<quad_vert[VERT_PER_QUAD*q] <<" "
							<<quad_vert[VERT_PER_QUAD*q+1] <<" "
                            <<quad_vert[VERT_PER_QUAD*q+2] <<" "
                            <<quad_vert[VERT_PER_QUAD*q+3] << endl;
                        cout <<"edgebase "
                            <<edgeBase[0]<<" "<<edgeBase[1]<<" "<<edgeBase[2]<<endl;
                        cout <<"dir "<< edgeDir <<endl;
					}
				}
			}
		}
		else
		{
			find_closest_edge_and_distance(DIM3, &(d1Coord[0]), edgeDir, d1EdgeBase, closestDistance);
			if (closestDistance <   epsilon)
			{
				if (scalar_grid.ComputeVertexIndex(&d1EdgeBase[0])
					== scalar_grid.ComputeVertexIndex(&edgeBase[0]))
				{
					if (is_permitted_collapse(iso_vlist, d1,  &(d1EdgeBase[0]),
						edgeDir, EDGE))
					{
						collapseToEdge(d1, q, quad_vert, vertex_coord, iso_vlist
							,collapse_map, edgeBase, edgeDir, dualiso_info);
					
						dualiso_info.cp_info.moved2Edge++;
						if ( printInfo)
						{
							cout <<"capquad "
								<<quad_vert[VERT_PER_QUAD*q] <<" "
								<<quad_vert[VERT_PER_QUAD*q+1] <<" "
								<<quad_vert[VERT_PER_QUAD*q+2] <<" "
								<<quad_vert[VERT_PER_QUAD*q+3] << endl;
							cout <<"edgebase "
								<<edgeBase[0]<<" "<<edgeBase[1]<<" "<<edgeBase[2]<<endl;
							cout <<"dir "<< edgeDir <<endl;
						}
					}
				}
			}
		}
		dualiso_info.cp_info.numCapQuad++;
		capQuad = true;

	}
	else
		capQuad = false;
}

/// Compute the edge dual to the quad q
void computeEdgeDual2q(
	const VERTEX_INDEX q, // q th quad
	const VERTEX_INDEX C,
	const std::vector<DUAL_ISOVERT> & iso_vlist,
	const std::vector<DIRECTION_TYPE> & orth_dir,
	IJK::ARRAY<GRID_COORD_TYPE> &edgeBase,
	int & d)
{
	edgeBase[0] = iso_vlist[C].cube_coord[0];
	edgeBase[1] = iso_vlist[C].cube_coord[1];
	edgeBase[2] = iso_vlist[C].cube_coord[2];
	d = (int)orth_dir[q];
}

/// Collapse caps
void collapse_caps (
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	const std::vector<DIRECTION_TYPE> & orth_dir,
	const float epsilon,
	std::vector<DUAL_ISOVERT> & iso_vlist,
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	DUALISO_INFO & dualiso_info, 
    const bool printInfo
	)
{
	const int numQuads = quad_vert.size()/4;
	const int notDegreeThree = scalar_grid.NumVertices(); // Notdegree3

	QTRIANGULATE::compute_degree_per_vertex(VERT_PER_QUAD, quad_vert, iso_vlist);

	for (int q = 0; q < numQuads; q++)
	{
		// the diagonals are given by AB and CD
		VERTEX_INDEX A = quad_vert[VERT_PER_QUAD*q];
		VERTEX_INDEX B = quad_vert[VERT_PER_QUAD*q+2];
		VERTEX_INDEX C = quad_vert[VERT_PER_QUAD*q+1];
		VERTEX_INDEX D = quad_vert[VERT_PER_QUAD*q+3];
		//update vertices
		A = find_vertex(collapse_map,A);
		B = find_vertex(collapse_map,B);
		C = find_vertex(collapse_map,C);
		D = find_vertex(collapse_map,D);
		bool flagCapQuad = false;


		IJK::ARRAY<GRID_COORD_TYPE> edgeBase(DIM3,0);
		IJK::ARRAY<COORD_TYPE> BCoord(DIM3,0);
		IJK::ARRAY<COORD_TYPE> ACoord(DIM3,0);
		
		for (int d=0;d<DIM3;d++)
		{
            ACoord[d]=vertex_coord[DIM3*A+d];
            BCoord[d]=vertex_coord[DIM3*B+d];
		}
		int edgeDir = -1;

		computeEdgeDual2q( q, B, iso_vlist, orth_dir, edgeBase, edgeDir);

		capQuadComputation(A,B, ACoord, BCoord, q, quad_vert, vertex_coord,
			scalar_grid, epsilon, notDegreeThree, iso_vlist, 
			edgeBase, edgeDir, flagCapQuad, collapse_map, dualiso_info, printInfo);


		capQuadComputation(C,D, ACoord,BCoord, q, quad_vert, vertex_coord,
			scalar_grid, epsilon, notDegreeThree, iso_vlist, 
			edgeBase, edgeDir, flagCapQuad, collapse_map, dualiso_info, printInfo);
	}


}

//Update the collapses in "collapse_map" in "quad_vert"
void update_quads(
	const IJK::ARRAY<VERTEX_INDEX> & collapse_map,
	std::vector<VERTEX_INDEX> & quad_vert
	)
{
	const int quad_vert_size = quad_vert.size();
	for (int v=0; v<quad_vert_size; v++)
	{
		int k = quad_vert[v];
		quad_vert[v]=collapse_map[k];
	}
}

// dual collapse main function
// NOTE:quads are reordered at the start and again at the end.
void QCOLLAPSE::dual_collapse(
	const DUALISO_DATA & dualiso_data,
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<DUAL_ISOVERT> & iso_vlist, 
	std::vector<COORD_TYPE> & vertex_coord,
	const std::vector<DIRECTION_TYPE> &orth_dir,
	const float epsilon,
	DUALISO_INFO & dualiso_info
	)
{

	//set up the vertex collapse map
	int num_vertex = vertex_coord.size();
	IJK::ARRAY<VERTEX_INDEX> collapse_map(num_vertex,0);
	//setup collapse_map.
	setup_collapse_map(collapse_map, num_vertex);
	//Reordering QuadVert 
	IJK::reorder_quad_vertices(quad_vert);

    if (dualiso_data.flag_cap_col)
	collapse_caps(scalar_grid, quad_vert, vertex_coord, orth_dir, 
		epsilon, iso_vlist, collapse_map, dualiso_info, dualiso_data.flag_collapse_debug);
	
	collapse_across_facets(scalar_grid, iso_vlist, quad_vert,
	vertex_coord, collapse_map, epsilon, dualiso_data.flag_collapse_debug, dualiso_info);

	collapse_across_edges(scalar_grid, iso_vlist, quad_vert,
	vertex_coord, collapse_map, epsilon, dualiso_data.flag_collapse_debug, dualiso_info);

	collapse_across_vertices(scalar_grid, iso_vlist, quad_vert, 
	vertex_coord, collapse_map, epsilon, dualiso_data.flag_collapse_debug, dualiso_info);
	
	update_quads(collapse_map, quad_vert);
	//reorder the quads back to the original.
	IJK::reorder_quad_vertices(quad_vert);
}

void QCOLLAPSE::delIsolated(
    std::vector<VERTEX_INDEX> & quad_vert,
	vector<VERTEX_INDEX> isolatedList,
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	DUALISO_INDEX_GRID & first_isov,
    IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
    bool printInfo
	)
{
    const int sizeDelIsolated = isolatedList.size();
	for (int v = 0; v < sizeDelIsolated; v++)
	{
        VERTEX_INDEX iv0 = isolatedList[v] - scalar_grid.CubeVertexIncrement(NUM_CUBE_VERT-1);
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
					if ( iso_vlist[indx_iso_vlist+k].sep_vert == isolatedList[v])
					{
                        iso_vlist[indx_iso_vlist+k].flag_isolated = true;
					}
				}
			}
		}
	}

    vector<VERTEX_INDEX>  local_quad_vert;
	const int num_quads = quad_vert.size()/NUM_CUBE_FACET_VERT;
	for (int i = 0; i < num_quads; i++)
	{
		bool flag_temp_isolated = false;
		for (int j = 0; j < VERT_PER_QUAD; j++)
		{
			VERTEX_INDEX v = quad_vert[NUM_CUBE_FACET_VERT*i +j];
			if (iso_vlist[v].flag_isolated)
			{
				flag_temp_isolated = true;
			}
		}

		if (!flag_temp_isolated)
		{
			for (int j = 0; j < VERT_PER_QUAD; j++)
			{
				VERTEX_INDEX v = quad_vert[NUM_CUBE_FACET_VERT*i +j];
				local_quad_vert.push_back(v);
			}    
		}
	}
   if (printInfo)
   {
    cout <<"Isolated list size: "<< sizeDelIsolated <<endl;
    cout <<"Number of quads decimated: "<< (quad_vert.size() - local_quad_vert.size()) / 4 <<endl;
   }
	quad_vert.clear();
    quad_vert = local_quad_vert;
}


