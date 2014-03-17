#include "qdual_collapse.h"
#include "qdual_remove_degen.h" // for the count degree function 

using namespace QCOLLAPSE;
using namespace std;
using namespace NamedConstants;

//Comparison
template <typename TY>
bool AreSame(TY a, TY b)
{
	return abs(a - b) < 0.00001;
}

// Compare two coords
template <typename TY>
bool compareCoords(IJK::ARRAY<TY> &coord1,
				   IJK::ARRAY<TY> &coord2)
{
	for (int i = 0; i < DIM3; i++)
	{
		if (!AreSame<TY>(coord1[i], coord2[i]))
		{
			return false;
		}
	}
	return true;
}


inline bool vert_simple (
	const IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	const int table_index
	)
{
	if (isodual_table.NumIsoVertices(table_index) == 1)
	{
		return true;
	}
	else
		return false;
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

//SubFunction to compute sigma called by
//CollapseToFacet
inline void computeSigma( const COORD_TYPE w1, const COORD_TYPE w2, 
						 const COORD_TYPE alpha, float &sigma, bool printInfo)
{
	COORD_TYPE wDiff=w2-w1;
	COORD_TYPE zero=0.0;
	if (AreSame<COORD_TYPE>(wDiff, zero)){
		sigma=0.5;
		if(printInfo)
			cout <<"sigma "<<sigma<< " wdiff "<<wDiff<<endl;
		return ;
	}
	else{
		sigma=(alpha-w1)/wDiff;
		if(printInfo)
		{
			cout <<"w1 "<< w1 <<" w2 "<< w2 <<" alpha "<<alpha <<endl;
			cout <<"sigma "<<sigma<< " wdiff "<<wDiff<<endl;
			if(sigma<0.0 || sigma >1.0)
			{
				cout <<"sigma is outside range"<< sigma<<endl;
			}
		}
	}

}
// Collapse to Facet
void collapseToFacet(
	const VERTEX_INDEX endPt1,
	const COORD_TYPE* endPt1Coord,
	const COORD_TYPE* endPt2Coord,
	IJK::ARRAY<GRID_COORD_TYPE> &facetBaseCoord,
	const int faceDir,
	std::vector<COORD_TYPE> & vertex_coord,
	std::vector<DUAL_ISOVERT> & iso_vlist,
	bool debug
	)
{
	float sigma=0.0;
	computeSigma(endPt1Coord[faceDir], endPt2Coord[faceDir],
		facetBaseCoord[faceDir], sigma, debug);

	if(debug)
	{
		cout <<"\nCollapsing to facet"<<endl;
		cout <<"endPT1 "<< endPt1Coord[0]<<" "<<endPt1Coord[1]<<" "<<endPt1Coord[2]<<endl;
		cout <<"endPT2 "<< endPt2Coord[0]<<" "<<endPt2Coord[1]<<" "<<endPt2Coord[2]<<endl;
		cout <<"facetBaseCoord "<<facetBaseCoord[0]<<" "<<facetBaseCoord[1]<<" "<<facetBaseCoord[2]<<endl;
		cout <<"faceDir "<< faceDir<<endl;
	}
	for (int d = 0; d < DIM3; d++)
	{
		vertex_coord[DIM3*endPt1+d]=endPt1Coord[d]*(1-sigma) + endPt2Coord[d]*sigma;
	}
	if(debug)
	{
		cout <<"new version of endPt1 "<<endl;
		cout <<vertex_coord[DIM3*endPt1+0]<<" "<<vertex_coord[DIM3*endPt1+1]<<" "<<vertex_coord[DIM3*endPt1+2]<<endl;
	}
	iso_vlist[endPt1].flag_fixed = true;
}
//Collapse across facets version B
void collapse_across_facetsB
	(	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	std::vector<DUAL_ISOVERT> & iso_vlist,
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	const std::vector<DIRECTION_TYPE> & orth_dir,
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	const float epsilon,
	const bool printInfo,
	DUALISO_INFO & dualiso_info)

{
	const int numQuads=quad_vert.size()/VERT_PER_QUAD;
	const int dimension=scalar_grid.Dimension();
	int v1,v2;
	IJK::ARRAY<GRID_COORD_TYPE> facetBaseCoord(dimension,0);
	IJK::ARRAY<GRID_COORD_TYPE> closestFacetBaseCoord(dimension,0);
	//foreach quad
	for (QUAD_INDEX q = 0; q < numQuads; q++){
		//foreach edge
		for (v1 = VERT_PER_QUAD-1, v2 = 0; v2 < VERT_PER_QUAD; v1 = v2++)
		{
			QUAD_INDEX endPt1=quad_vert[q*VERT_PER_QUAD+v1];
			QUAD_INDEX endPt2=quad_vert[q*VERT_PER_QUAD+v2];
			endPt1=find_vertex(collapse_map, endPt1);
			endPt2=find_vertex(collapse_map, endPt2);
			if(endPt1!=endPt2){
				const COORD_TYPE* endPt1Coord=&(vertex_coord[DIM3*endPt1]);
				const COORD_TYPE* endPt2Coord=&(vertex_coord[DIM3*endPt2]);
				//DEBUG  GREEN
				if(printInfo){
					cout <<"\n collapse across facets quad "<< q<<"  endpts "<< endPt1 <<" "<<endPt2 <<endl;
					cout <<"endpt coords ("<<endPt1<<") "<< endPt1Coord[0] <<" "<< endPt1Coord[1]<<" "<<endPt1Coord[2]<<endl;
					cout <<" ("<<endPt2<<") "<< endPt2Coord[0]<<" "<<endPt2Coord[1]<<" "<<endPt2Coord[2]<<endl;
					cout <<"is Endpt "<< endPt2 <<"  fixed ? "
						<< iso_vlist[endPt2].flag_fixed <<endl;
					cout <<"is Endpt "<< endPt1 <<"  fixed ? "
						<< iso_vlist[endPt1].flag_fixed <<endl;
				}
				int numClose=0;
				int closestFacetDir;
				bool endPt1CloseToAFacet = true;
				for (int d = 0; d < dimension; d++)
				{
					find_closest_facet(DIM3, endPt1Coord, d, facetBaseCoord);
					if(is_epsilon_close(endPt1Coord, facetBaseCoord, d, epsilon))
					{
						//DEBUG
						if(printInfo){
							cout <<"tempclosestFacetBaseCoord to endPt1 ("<< closestFacetBaseCoord[0]<<" "
								<<closestFacetBaseCoord[1]<<" "
								<<closestFacetBaseCoord[2]<<") ";
							cout <<"numclose "<< numClose <<endl;
							cout <<"closestDir "<< closestFacetDir<<endl;
						}
						numClose++;
						if (numClose > 1)
						{
							endPt1CloseToAFacet = false;
							break;
						}
						closestFacetDir=d;
						closestFacetBaseCoord[0]=facetBaseCoord[0];
						closestFacetBaseCoord[1]=facetBaseCoord[1];
						closestFacetBaseCoord[2]=facetBaseCoord[2];

					}
				}
				//DEBUG
				if(printInfo){
					cout <<"closestFacetBaseCoord to endPt1 ("<< closestFacetBaseCoord[0]<<" "
						<<closestFacetBaseCoord[1]<<" "
						<<closestFacetBaseCoord[2]<<") ";
					cout <<"numclose "<< numClose <<endl;
					cout <<"closestDir "<< closestFacetDir<<endl;
				}
				bool debugFlagA=false, debugFlagB=false, debugFlagC=false;
				if(numClose!=0 && endPt1CloseToAFacet)
				{
					bool endPt2CloseToAFacet = true;

					//endPt1 is close to a grid facet and not close to any edge.
					find_closest_facet(DIM3, endPt2Coord, closestFacetDir, facetBaseCoord);
					if(is_epsilon_close(endPt2Coord, facetBaseCoord, closestFacetDir, epsilon))
					{
						debugFlagA = true;
						if(compareCoords<GRID_COORD_TYPE>(facetBaseCoord, closestFacetBaseCoord))
						{
							debugFlagB=true;
							const int d1 = (closestFacetDir+1)%DIM3;
							const int d2 = (closestFacetDir+2)%DIM3;
							//bool endPt2CloseToEdge = false;
							find_closest_facet(DIM3, endPt2Coord, d1, facetBaseCoord);
							if(is_epsilon_close(endPt2Coord, facetBaseCoord, d1, epsilon))
							{
								//endPt2CloseToEdge=true;
								endPt2CloseToAFacet = false;
							}
							if( endPt2CloseToAFacet )
							{
								find_closest_facet(DIM3, endPt2Coord, d2, facetBaseCoord);
								if(is_epsilon_close(endPt2Coord, facetBaseCoord, d2, epsilon))
								{
									//endPt2CloseToEdge=true;
									endPt2CloseToAFacet = false;
								}
							}

							if(endPt2CloseToAFacet)
							{
								debugFlagC=true;
								//endPt2  is close to grid facet in the direction  closestFacetDir
								//and not close to any edge
								if (is_permitted_collapse (iso_vlist, endPt2,  &(closestFacetBaseCoord[0]),
									closestFacetDir,  FACET)
									&&
									is_permitted_collapse (iso_vlist, endPt1,  &(closestFacetBaseCoord[0]),
									closestFacetDir,  FACET)
									)
								{
									if (printInfo){
										print_collapse_info("\ncollapse across facets [permitted]", 
											endPt2, endPt1, &(closestFacetBaseCoord [0]));
										cout <<"endPT2Coord "<< endPt2Coord[0]<<" "<<endPt2Coord[1]<<" "<<endPt2Coord[2]<<endl;
										cout <<"endPT1Coord "<< endPt1Coord[0]<<" "<<endPt1Coord[1]<<" "<<endPt1Coord[2]<<endl;

									}
									collapseToFacet(endPt1, endPt1Coord, endPt2Coord, 
										closestFacetBaseCoord, closestFacetDir,
										vertex_coord, iso_vlist,  printInfo);
									dualiso_info.col_info.permitted_facet_restriction++;
									update_collapse_edges(collapse_map,endPt1, endPt2);
								}
								else
								{
									dualiso_info.col_info.not_permitted_facet_restriction++;
									if (printInfo){
										print_collapse_info("collapse across facets [NOT permitted]", 
											endPt2, endPt1, &(closestFacetBaseCoord [0]));
									}
								}
							}
						}
					}
				}//
				//debug
				if (printInfo)
				{
					if (debugFlagA)
						cout <<"endpt2 is epsilon close to  facetBase  in dir of closest facet" << endl;
					if(debugFlagB)
						cout <<" the are the same vertex"<<endl;
					if(debugFlagC)
						cout <<"endpt2 is not close to an edge or vertex "<<endl;
				}
			}
		}
	}
}


//subfunction for collapse across facet version C
//Is EndPt epsilon close to a grid edge. 
bool isEpsilonCloseToGridEdge
	( 
	const VERTEX_INDEX endPt,
	const COORD_TYPE* endPtCoord,
	const float epsilon,
	//returns
	IJK::ARRAY<GRID_COORD_TYPE> &closestEdgeBase2EndPt,
	DIRECTION_TYPE & closestEdge2EndPt1Dir
	)
{
	IJK::ARRAY<GRID_COORD_TYPE> tempEdgeBaseCoord(DIM3, 0);
	float tempDistance=0.0;
	float closestEdgeDist2Endpt=0.0;
	for (int d = 0; d < DIM3; d++)
	{
		find_closest_edge_and_distance
			(DIM3, endPtCoord, d, tempEdgeBaseCoord, tempDistance);
		if(d==0)
		{
			closestEdge2EndPt1Dir = 0;
			closestEdgeDist2Endpt = tempDistance;
			closestEdgeBase2EndPt[0] = tempEdgeBaseCoord[0];
			closestEdgeBase2EndPt[1] = tempEdgeBaseCoord[1];
			closestEdgeBase2EndPt[2] = tempEdgeBaseCoord[2];
		}
		if (tempDistance < closestEdgeDist2Endpt )
		{
			closestEdge2EndPt1Dir = d;
			closestEdgeDist2Endpt = tempDistance;
			closestEdgeBase2EndPt[0] = tempEdgeBaseCoord[0];
			closestEdgeBase2EndPt[1] = tempEdgeBaseCoord[1];
			closestEdgeBase2EndPt[2] = tempEdgeBaseCoord[2];
		}
	}
	if (closestEdgeDist2Endpt <= epsilon)
	{	
		return true;
	}
	else
	{
		return false;
	}
}



//Collapse across facets
void collapse_across_facets(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
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
				bool vert_simple1 = vert_simple(isodual_table, iso_vlist[endPt1].table_index);
				bool vert_simple2 = vert_simple(isodual_table, iso_vlist[endPt2].table_index);
				if(vert_simple1 == vert_simple2)
				{
					if (iso_vlist[endPt2].cube_index < iso_vlist[endPt1].cube_index)
					{
						swap(endPt1, endPt2);
					}
				}
				else if (vert_simple2){
					swap(endPt1, endPt2);
				}
				const COORD_TYPE * endPt1_coord = & (vertex_coord[DIM3*endPt1]);
				const COORD_TYPE * endPt2_coord = & (vertex_coord[DIM3*endPt2]);
				//DEBUG  GREEN
				if(print_info){
					cout <<"\nquad "<< q<<"  endpts "<< endPt1 <<" "<<endPt2 <<endl;
					cout <<"endpt coords ("<<endPt1<<") "<< endPt1_coord[0] <<" "<< endPt1_coord[1]<<" "<<endPt1_coord[2]<<endl;
					cout <<" ("<<endPt2<<") "<< endPt2_coord[0]<<" "<<endPt2_coord[1]<<" "<<endPt2_coord[2]<<endl;
				}
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
								closest_facet,  FACET)
								&& (iso_vlist[endPt2].flag_fixed == false))
							{ 

								if (print_info){
									print_collapse_info("collapse across facets [permitted]", 
										endPt1, endPt2, &(facet_base_coord[0]));
									cout <<"endpt " << endPt1 <<" is close to a facet."<<endl;
									cout <<"facet direction "<< closest_facet <<endl;
								}

								dualiso_info.col_info.permitted_facet_restriction++;
								update_collapse_edges(collapse_map,endPt1, endPt2);
							}
							else
							{
								dualiso_info.col_info.not_permitted_facet_restriction++;
								if (print_info)
								{
									print_collapse_info("collapse across facets [NOT permitted]", 
										endPt1, endPt2, &(facet_base_coord[0]));

									cout <<"is Endpt "<< endPt2 <<"  fixed ? "
										<< iso_vlist[endPt2].flag_fixed <<endl;
									cout <<"endpt " << endPt1 <<" is close to a facet."<<endl;
									cout <<"facet direction "<< closest_facet <<endl;
								}
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
									closest_facet, FACET)
									&& (iso_vlist[endPt1].flag_fixed == false))
								{ 
									if (print_info){
										print_collapse_info("collapse across facets [permitted]", 
											endPt2, endPt1, &(facet_base_coord[0]));
										cout <<"facet direction "<< closest_facet <<endl;
										cout <<"endpt " << endPt2 <<" is close to a facet."<<endl;
									}
									dualiso_info.col_info.permitted_facet_restriction++;
									update_collapse_edges(collapse_map,endPt2, endPt1);
								}
								else{
									dualiso_info.col_info.not_permitted_facet_restriction++;
									if (print_info)
									{
										print_collapse_info("collapse across facets [NOT permitted]", 
											endPt2, endPt1, &(facet_base_coord[0]));

										cout <<"is Endpt "<< endPt1 <<"  fixed ? "
											<< iso_vlist[endPt1].flag_fixed <<endl;
										cout <<"endpt " << endPt2 <<" is close to a facet."<<endl;
										cout <<"facet direction "<< closest_facet <<endl;
									}
								}
							}
						}
					}
				}
			}//if, only if the endPts are not already contracted		
		} // for each edge end 
	}
}

//sub function for collapse across edgesB
bool isClose2Edge(
	const VERTEX_INDEX endPt,
	const COORD_TYPE* endPtCoord,
	const DIRECTION_TYPE d, 
	const COORD_TYPE epsilon,
	const bool debug
	)
{
	int numCloseFacetsOtherThan_d=0;
	DIRECTION_TYPE d1 = (d+1)%DIM3;
	DIRECTION_TYPE d2 = (d+2)%DIM3;
	IJK::ARRAY<GRID_COORD_TYPE> facetBaseCoord(DIM3,0);
	find_closest_facet(DIM3, endPtCoord, d1, facetBaseCoord);
	if(is_epsilon_close(endPtCoord, facetBaseCoord, d1, epsilon))
	{
		numCloseFacetsOtherThan_d++;
	}
	find_closest_facet(DIM3, endPtCoord, d2, facetBaseCoord);
	if(is_epsilon_close(endPtCoord, facetBaseCoord, d2, epsilon))
	{
		numCloseFacetsOtherThan_d++;
	}
	if (numCloseFacetsOtherThan_d==1)
	{
		return true;
	}
	else
		return false;
}
//Sub function for collapse across edgesB
//Returns the facetBaseCoord for the common facet and the facetDir.
bool findFacetCommonToEndPts(
	const VERTEX_INDEX endPt1,
	const VERTEX_INDEX endPt2,
	const std::vector<COORD_TYPE> & vertex_coord,
	const std::vector<DUAL_ISOVERT> & iso_vlist,
	const bool debug,
	//returns
	IJK::ARRAY<GRID_COORD_TYPE> &facetBaseCoord,
	DIRECTION_TYPE & facetDir
	)
{
	const COORD_TYPE * cube0Coord = &(iso_vlist[endPt1].cube_coord[0]);
	const COORD_TYPE * cube1Coord = &(iso_vlist[endPt2].cube_coord[0]);
	//find facetDir
	facetDir=-1;
	int count=0;
	for (int d = 0; d < DIM3; d++)
	{
		if (!AreSame<COORD_TYPE>(cube0Coord[d],cube1Coord[d]))
		{
			facetDir=d;
			count++;
		}
	}

	if (count==1)
	{
		DIRECTION_TYPE d1 = (facetDir+1)%DIM3;
		DIRECTION_TYPE d2 = (facetDir+2)%DIM3;
		facetBaseCoord[d1]= cube0Coord[d1];
		facetBaseCoord[d2]= cube0Coord[d2];
		facetBaseCoord[facetDir] = max(cube0Coord[facetDir], cube1Coord[facetDir]);
		return true;
	}
	else
		return false;
}

//Collapse across facets version C
void collapse_across_facetsC
	(	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	std::vector<DUAL_ISOVERT> & iso_vlist,
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	const std::vector<DIRECTION_TYPE> & orth_dir,
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	const float epsilon,
	const bool printInfo,
	DUALISO_INFO & dualiso_info)
{
	const int numQuads=quad_vert.size()/VERT_PER_QUAD;
	const int dimension=scalar_grid.Dimension();
	int v1,v2;
	IJK::ARRAY<GRID_COORD_TYPE> facetBaseCoord(dimension,0);
	IJK::ARRAY<GRID_COORD_TYPE> closestFacetBaseCoord(dimension,0);
	//foreach quad
	for (QUAD_INDEX q = 0; q < numQuads; q++){
		//foreach edge
		for (v1 = VERT_PER_QUAD-1, v2 = 0; v2 < VERT_PER_QUAD; v1 = v2++)
		{
			QUAD_INDEX endPt1=quad_vert[q*VERT_PER_QUAD+v1];
			QUAD_INDEX endPt2=quad_vert[q*VERT_PER_QUAD+v2];
			if(printInfo){
				cout <<"\nCollapseAcrossFacets quad "<< q<<"  endpts "<< endPt1 <<" "<<endPt2 <<endl;
			}
			endPt1=find_vertex(collapse_map, endPt1);
			endPt2=find_vertex(collapse_map, endPt2);

			if(endPt1!=endPt2){
				const COORD_TYPE* endPt1Coord=&(vertex_coord[DIM3*endPt1]);
				const COORD_TYPE* endPt2Coord=&(vertex_coord[DIM3*endPt2]);
				IJK::ARRAY<GRID_COORD_TYPE> facetBaseCoord(DIM3,0);
				DIRECTION_TYPE facetDir=-1;
				bool hasCommonFacet = findFacetCommonToEndPts
					(endPt1, endPt2,  vertex_coord, iso_vlist, printInfo, facetBaseCoord, facetDir);
				if(printInfo){
					cout <<"Updated endpts "<< endPt1 <<" "<<endPt2 <<endl;
					cout <<"endpt coords ("<<endPt1<<") "<< endPt1Coord[0] <<" "<< endPt1Coord[1]<<" "<<endPt1Coord[2]<<endl;
					cout <<" ("<<endPt2<<") "<< endPt2Coord[0]<<" "<<endPt2Coord[1]<<" "<<endPt2Coord[2]<<endl;
					cout <<"facetBaseCoord "<< facetBaseCoord[0]<<" "<<facetBaseCoord[1]<<" "<<facetBaseCoord[2]<<endl;
				}
				if (hasCommonFacet){
					bool endPt1EpsilonClose2F =
						is_epsilon_close(endPt1Coord, facetBaseCoord, facetDir, epsilon);
					bool endPt2EpsilonClose2F =
						is_epsilon_close(endPt2Coord, facetBaseCoord, facetDir, epsilon);
					bool endPt1IsPermittedAcrossF=
						is_permitted_collapse (iso_vlist, endPt1,  &(facetBaseCoord[0]),
						facetDir,  FACET);
					bool endPt2IsPermittedAcrossF=
						is_permitted_collapse (iso_vlist, endPt2,  &(facetBaseCoord[0]),
						facetDir,  FACET);
					if(endPt1EpsilonClose2F && endPt2EpsilonClose2F 
						&& endPt1IsPermittedAcrossF && endPt2IsPermittedAcrossF)
					{
						IJK::ARRAY<GRID_COORD_TYPE> closestEdgeBase2EndPt1(DIM3, 0);
						DIRECTION_TYPE closestEdge2EndPt1Dir=-1;
						bool isEndPt1EpsilonCloseToGridEdge = isEpsilonCloseToGridEdge 
							( endPt1, endPt1Coord, epsilon, closestEdgeBase2EndPt1, closestEdge2EndPt1Dir);

						IJK::ARRAY<GRID_COORD_TYPE> closestEdgeBase2EndPt2(DIM3, 0);
						DIRECTION_TYPE closestEdge2EndPt2Dir=-1;
						bool isEndPt2EpsilonCloseToGridEdge = isEpsilonCloseToGridEdge 
							( endPt2, endPt2Coord, epsilon, closestEdgeBase2EndPt2, closestEdge2EndPt2Dir);

						if ((isEndPt1EpsilonCloseToGridEdge == false)
							&&(isEndPt2EpsilonCloseToGridEdge == false))
						{
							collapseToFacet
								( endPt1, endPt1Coord, endPt2Coord, facetBaseCoord,
								facetDir, vertex_coord, iso_vlist, printInfo);
							dualiso_info.col_info.permitted_facet_restriction++;
							update_collapse_edges(collapse_map,endPt1, endPt2);
							if(printInfo)
							{
								cout <<"endpt "<< endPt1 <<"  and endpt "<< endPt2 <<" are not close to any edge\n";
								cout <<isEndPt1EpsilonCloseToGridEdge<<" "<<isEndPt2EpsilonCloseToGridEdge<<endl;
								cout <<"**endpt "<< endPt2 <<" collapses to "<< endPt1 << endl;
							}
						}
						else if(isEndPt1EpsilonCloseToGridEdge==false)
						{
							update_collapse_edges(collapse_map,endPt1, endPt2);
							dualiso_info.col_info.permitted_edge_collapse++;
							iso_vlist[endPt1].flag_fixed=true;
							if(printInfo)
							{
								cout <<"endpt "<< endPt1 <<"  and endpt "<< endPt2 <<". endpt1 is not close to any edge\n";
								cout <<"**endpt "<< endPt2 <<" collapses to "<< endPt1 << endl;
							}
						}
						else if(isEndPt2EpsilonCloseToGridEdge==false)
						{
							update_collapse_edges(collapse_map,endPt2, endPt1);
							dualiso_info.col_info.permitted_edge_collapse++;
							iso_vlist[endPt2].flag_fixed=true;
							if(printInfo)
							{
								cout <<"endpt "<< endPt1 <<"  and endpt "<< endPt2 <<". endpt2 is not close to any edge\n";
								cout <<"**endpt "<< endPt1 <<" collapses to "<< endPt2 << endl;
							}
						}
					}
				}
			}
		}
	}
}

//Collapse across edges version B
void collapse_across_edgesB(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	std::vector<DUAL_ISOVERT> & iso_vlist,
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	const float epsilon,
	const bool printInfo,
	DUALISO_INFO & dualiso_info)
{
	const int numQuads = quad_vert.size()/VERT_PER_QUAD;
	int v1, v2;
	for (int q = 0; q < numQuads; q++)
	{
		for (int v1 = VERT_PER_QUAD-1, v2 = 0; v2 < VERT_PER_QUAD; v1 = v2++)
		{
			int endPt1 = quad_vert[VERT_PER_QUAD*q+v1];
			int endPt2 = quad_vert[VERT_PER_QUAD*q+v2];

			//update the endpts
			endPt1=find_vertex(collapse_map, endPt1);
			endPt2=find_vertex(collapse_map, endPt2);
			const COORD_TYPE* endPt1Coord=&(vertex_coord[DIM3*endPt1]);
			const COORD_TYPE* endPt2Coord=&(vertex_coord[DIM3*endPt2]);

			if ( (endPt1!=endPt2) && (iso_vlist[endPt2].flag_fixed == false)
				&& (iso_vlist[endPt1].flag_fixed == false))
			{
				IJK::ARRAY<GRID_COORD_TYPE> facetBaseCoord(DIM3,0);
				DIRECTION_TYPE facetDir=-1;
				bool hasCommonFacet = findFacetCommonToEndPts
					(endPt1, endPt2,  vertex_coord, iso_vlist, printInfo, facetBaseCoord, facetDir);
				if(hasCommonFacet){
					bool endPt1EpsilonClose2F =
						is_epsilon_close(endPt1Coord, facetBaseCoord, facetDir, epsilon);
					bool endPt2EpsilonClose2F =
						is_epsilon_close(endPt2Coord, facetBaseCoord, facetDir, epsilon);
					bool endPt1IsPermittedAcrossF=
						is_permitted_collapse (iso_vlist, endPt1,  &(facetBaseCoord[0]),
						facetDir,  FACET);
					bool endPt2IsPermittedAcrossF=
						is_permitted_collapse (iso_vlist, endPt2,  &(facetBaseCoord[0]),
						facetDir,  FACET);
					if(endPt1EpsilonClose2F && endPt2EpsilonClose2F 
						&& endPt1IsPermittedAcrossF && endPt2IsPermittedAcrossF)
					{
						//find Grid edge closest to endPt1
						IJK::ARRAY<GRID_COORD_TYPE> tempEdgeBaseCoord(DIM3, 0);
						IJK::ARRAY<GRID_COORD_TYPE> closestEdgeBase2EndPt1(DIM3, 0);
						float tempDistance=0.0;
						DIRECTION_TYPE closestEdge2EndPt1Dir=-1;
						float closestEdgeDist2Endpt1=0.0;
						for (int d = 0; d < DIM3; d++)
						{
							find_closest_edge_and_distance
								(DIM3, endPt1Coord, d, tempEdgeBaseCoord, tempDistance);
							if(d==0)
							{
								closestEdge2EndPt1Dir = 0;
								closestEdgeDist2Endpt1 = tempDistance;
								closestEdgeBase2EndPt1[0] = tempEdgeBaseCoord[0];
								closestEdgeBase2EndPt1[1] = tempEdgeBaseCoord[1];
								closestEdgeBase2EndPt1[2] = tempEdgeBaseCoord[2];
							}
							if (tempDistance < closestEdgeDist2Endpt1 )
							{
								closestEdge2EndPt1Dir = d;
								closestEdgeDist2Endpt1 = tempDistance;
								closestEdgeBase2EndPt1[0] = tempEdgeBaseCoord[0];
								closestEdgeBase2EndPt1[1] = tempEdgeBaseCoord[1];
								closestEdgeBase2EndPt1[2] = tempEdgeBaseCoord[2];
							}
						}
						if (closestEdgeDist2Endpt1 <= epsilon)
						{
							find_closest_edge_and_distance
								(DIM3, endPt2Coord, closestEdge2EndPt1Dir,  tempEdgeBaseCoord, tempDistance);
							if ( tempDistance <= epsilon)
							{
								if (compareCoords<GRID_COORD_TYPE>( tempEdgeBaseCoord, closestEdgeBase2EndPt1))
								{
									collapseToFacet
										(endPt1, endPt1Coord, endPt2Coord, closestEdgeBase2EndPt1,
										closestEdge2EndPt1Dir, vertex_coord, iso_vlist, printInfo);
									update_collapse_edges(collapse_map,endPt1, endPt2);
									dualiso_info.col_info.permitted_edge_collapse++;
								}
							}
						}
						/*
						/// DEPRECATED ** NOT WORKING **
						bool endPt1Close2Edge=
						isClose2Edge(endPt1, endPt1Coord, facetDir, epsilon, printInfo);
						bool endPt2Close2Edge=
						isClose2Edge(endPt2, endPt2Coord, facetDir, epsilon, printInfo);

						if (endPt1Close2Edge && !endPt2Close2Edge)
						{
						if (printInfo){
						print_collapse_info("\ncollapse across edges [permitted]", 
						endPt2, endPt1, &(facetBaseCoord [0]));
						cout <<"endPT2Coord "<< endPt2Coord[0]<<" "<<endPt2Coord[1]<<" "<<endPt2Coord[2]<<endl;
						cout <<"endPT1Coord "<< endPt1Coord[0]<<" "<<endPt1Coord[1]<<" "<<endPt1Coord[2]<<endl;
						cout <<"is Endpt "<< endPt2 <<"  fixed ? "
						<< iso_vlist[endPt2].flag_fixed <<endl;
						cout <<"is Endpt "<< endPt1 <<"  fixed ? "
						<< iso_vlist[endPt1].flag_fixed <<endl;
						}
						collapseToFacet(endPt1, endPt1Coord, endPt2Coord, facetBaseCoord, facetDir, 
						vertex_coord,  iso_vlist, printInfo);
						dualiso_info.col_info.permitted_edge_collapse++;
						update_collapse_edges(collapse_map,endPt1, endPt2);
						}
						if (endPt2Close2Edge && !endPt1Close2Edge)
						{
						if (printInfo){
						print_collapse_info("\ncollapse across edges [permitted]", 
						endPt2, endPt1, &(facetBaseCoord [0]));
						cout <<"endPT2Coord "<< endPt2Coord[0]<<" "<<endPt2Coord[1]<<" "<<endPt2Coord[2]<<endl;
						cout <<"endPT1Coord "<< endPt1Coord[0]<<" "<<endPt1Coord[1]<<" "<<endPt1Coord[2]<<endl;
						cout <<"is Endpt "<< endPt2 <<"  fixed ? "
						<< iso_vlist[endPt2].flag_fixed <<endl;
						cout <<"is Endpt "<< endPt1 <<"  fixed ? "
						<< iso_vlist[endPt1].flag_fixed <<endl;
						}
						collapseToFacet(endPt1, endPt1Coord, endPt2Coord, facetBaseCoord, facetDir, 
						vertex_coord, iso_vlist, printInfo);
						update_collapse_edges(collapse_map,endPt1, endPt2);
						dualiso_info.col_info.permitted_edge_collapse++;
						}
						*/

					}
				}
			}
		}
	}
}


//Collapse across edges
void collapse_across_edges(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
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
				bool vert_simple1 = vert_simple(isodual_table, iso_vlist[endPt1].table_index);
				bool vert_simple2 = vert_simple(isodual_table, iso_vlist[endPt2].table_index);
				if(vert_simple1 == vert_simple2){
					if (iso_vlist[endPt2].cube_index < iso_vlist[endPt1].cube_index){
						swap(endPt1, endPt2);
					}
				}
				else if (vert_simple2){
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
								closest_edge_dir, EDGE)
								&& (iso_vlist[endPt2].flag_fixed == false))
							{ 
								if (print_info){
									print_collapse_info("collapse across edges [permitted]", 
										endPt1, endPt2, &(edge_base_coord[0]));
									cout <<"endPt1 "<< endPt1<<" "<< iso_vlist[endPt1].flag_fixed <<endl;
									cout <<"endPt2 "<< endPt2<<" "<< iso_vlist[endPt2].flag_fixed <<endl;
									cout <<"endpt2 "<< endPt2<<"("<< endPt2_coord[0]<<","<< endPt2_coord[1]<<","<< endPt2_coord[2]<<"\n";
									cout <<"endpt1 "<< endPt1<<"("<< endPt1_coord[0]<<","<< endPt1_coord[1]<<","<< endPt1_coord[2]<<"\n";
									cout <<"edge base "<< edge_base_coord[0]<<","
										<< edge_base_coord[1]<<"," << edge_base_coord[2] << endl;
									cout <<"edge dir "<< closest_edge_dir <<"\n"<<endl;
								}

								dualiso_info.col_info.permitted_edge_collapse++;
								update_collapse_edges(collapse_map,endPt1, endPt2);
							}
							else
							{
								dualiso_info.col_info.not_permitted_edge_restriction++;
								if (print_info){
									print_collapse_info("collapse across edges [NOT permitted]", 
										endPt1, endPt2, &(edge_base_coord[0]));
									cout <<"endPt1 "<< endPt1<<" "<< iso_vlist[endPt1].flag_fixed <<endl;
									cout <<"endPt2 "<< endPt2<<" "<< iso_vlist[endPt2].flag_fixed <<endl;

									cout <<scalar_grid.ComputeVertexIndex(&edge_base_coord[0]);
									cout <<" "<<scalar_grid.ComputeVertexIndex(&closest_edge_base_coord[0]) <<endl;
									int v=0;
									v=scalar_grid.ComputeVertexIndex(&edge_base_coord[0]);
									cout <<"edgebase vertex index "<<v<<endl;

									cout <<"endpt2 "<< endPt2<<"("<< endPt2_coord[0]<<","<< endPt2_coord[1]<<","<< endPt2_coord[2]<<"\n";
									cout <<"endpt1 "<< endPt1<<"("<< endPt1_coord[0]<<","<< endPt1_coord[1]<<","<< endPt1_coord[2]<<"\n";
									cout <<"edge base "<< edge_base_coord[0]<<","
										<< edge_base_coord[1]<<"," << edge_base_coord[2] << endl;
									cout <<"edge dir "<< closest_edge_dir <<"\n"<<endl;
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
									closest_edge_dir, EDGE) 
									&& (iso_vlist[endPt1].flag_fixed == false))
								{ 
									if (print_info){
										print_collapse_info("collapse across edges [permitted]", 
											endPt2, endPt1, &(edge_base_coord[0]));
										cout <<"endPt1 "<< endPt1<<" "<< iso_vlist[endPt1].flag_fixed <<endl;
										cout <<"endPt2 "<< endPt2<<" "<< iso_vlist[endPt2].flag_fixed <<endl;
										int v=0;
										v=scalar_grid.ComputeVertexIndex(&edge_base_coord[0]);
										cout <<"edgebase vertex index "<<v<<endl;

										cout <<"endpt1 "<< endPt1<<"("<< endPt1_coord[0]<<","<< endPt1_coord[1]<<","<< endPt1_coord[2]<<"\n";

										cout <<"edge base "<< edge_base_coord[0]<<","
											<< edge_base_coord[1]<<"," << edge_base_coord[2] << endl;
										cout <<"edge dir "<< closest_edge_dir <<endl;
									}

									dualiso_info.col_info.permitted_edge_collapse++;
									update_collapse_edges(collapse_map,endPt2, endPt1);
								}
								else
								{
									dualiso_info.col_info.not_permitted_edge_restriction++;
									if (print_info){
										print_collapse_info("collapse across edges [ NOT permitted]", 
											endPt2, endPt1, &(edge_base_coord[0]));
										cout <<"else case"<<endl;
										cout <<"endPt1 "<< endPt1<<" "<< iso_vlist[endPt1].flag_fixed <<endl;
										cout <<"endPt2 "<< endPt2<<" "<< iso_vlist[endPt2].flag_fixed <<endl;
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
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
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
					bool vert_simple1 = vert_simple(isodual_table, iso_vlist[endPt1].table_index);
					bool vert_simple2 = vert_simple(isodual_table, iso_vlist[endPt2].table_index);
					if(vert_simple1 == vert_simple2)
					{
						if (iso_vlist[endPt2].cube_index < iso_vlist[endPt1].cube_index)
						{
							swap(endPt1, endPt2);
						}
					}
					else if (vert_simple2){
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
						else
						{
							if(print_info)
							{
								cout <<" closest distance  endPt2 "<<closest_distance<<" greater than "<< epsilon << endl;
							}
						}
					}
					else
					{
						if(print_info)
						{
							cout <<" closest distance endPt1 "<<closest_distance<<" greater than "<< epsilon << endl;
						}
					}
				}
			}//for each edge
}

//collapse across vertices
void collapse_across_verticesB(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	std::vector<DUAL_ISOVERT> & iso_vlist,
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	const float epsilon,
	const bool print_info,
	DUALISO_INFO & dualiso_info)
{
	const int numQuads = quad_vert.size()/ VERT_PER_QUAD;
	int v1, v2;
	for (int q = 0; q < numQuads; q++)
	{
		for (v1 = VERT_PER_QUAD-1, v2=0; v2 < VERT_PER_QUAD; v1 = v2++)
		{
			int endPt1 = quad_vert[q*4+v1];
			int endPt2 = quad_vert[q*4+v2];
			endPt1 = find_vertex(collapse_map, endPt1);
			endPt2 = find_vertex(collapse_map, endPt2);

			if ( (endPt1!=endPt2) && (iso_vlist[endPt2].flag_fixed == false)
				&& (iso_vlist[endPt1].flag_fixed == false)
				&& (iso_vlist[endPt1].ver_degree == 3)
				&& (iso_vlist[endPt2].ver_degree == 3)
				&& (iso_vlist[endPt1].sep_vert == iso_vlist[endPt2].sep_vert)
				)
			{
				const COORD_TYPE* endPt1Coord=&(vertex_coord[DIM3*endPt1]);
				const COORD_TYPE* endPt2Coord=&(vertex_coord[DIM3*endPt2]);
				IJK::ARRAY<GRID_COORD_TYPE> facetBaseCoord(DIM3,0);
				DIRECTION_TYPE facetDir=-1;
				bool hasCommonFacet = findFacetCommonToEndPts
					(endPt1, endPt2,  vertex_coord, iso_vlist, print_info, facetBaseCoord, facetDir);
				if(hasCommonFacet)
				{
					bool endPt1EpsilonClose2F =
						is_epsilon_close(endPt1Coord, facetBaseCoord, facetDir, epsilon);
					bool endPt2EpsilonClose2F =
						is_epsilon_close(endPt2Coord, facetBaseCoord, facetDir, epsilon);
					bool endPt1IsPermittedAcrossF=
						is_permitted_collapse (iso_vlist, endPt1,  &(facetBaseCoord[0]),
						facetDir,  FACET);
					bool endPt2IsPermittedAcrossF=
						is_permitted_collapse (iso_vlist, endPt2,  &(facetBaseCoord[0]),
						facetDir,  FACET);
					if(endPt1EpsilonClose2F && endPt2EpsilonClose2F 
						&& endPt1IsPermittedAcrossF && endPt2IsPermittedAcrossF)
					{
						if (print_info){
							print_collapse_info("\ncollapse across vertices [permitted]", 
								endPt2, endPt1, &(facetBaseCoord [0]));
							cout <<"endPT2Coord "<< endPt2Coord[0]<<" "<<endPt2Coord[1]<<" "<<endPt2Coord[2]<<endl;
							cout <<"endPT1Coord "<< endPt1Coord[0]<<" "<<endPt1Coord[1]<<" "<<endPt1Coord[2]<<endl;
							cout <<"is Endpt "<< endPt2 <<"  fixed ? "
								<< iso_vlist[endPt2].flag_fixed <<endl;
							cout <<"is Endpt "<< endPt1 <<"  fixed ? "
								<< iso_vlist[endPt1].flag_fixed <<endl;
						}
						collapseToFacet(endPt1, endPt1Coord, endPt2Coord, facetBaseCoord, facetDir, 
							vertex_coord, iso_vlist, print_info);
						update_collapse_edges(collapse_map,endPt1, endPt2);
						dualiso_info.col_info.permitted_vertex_restriction++;
					}
				}
			}
		}
	}
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
	std::vector<DUAL_ISOVERT> & iso_vlist,
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
			update_collapse_edges(collapse_map, v, v2);
		}
	}
	iso_vlist[v].flag_fixed = true;
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
	std::vector<DUAL_ISOVERT> & iso_vlist,
	IJK::ARRAY<GRID_COORD_TYPE> &edgeBase,
	int & edgeDir,
	bool &capQuad,
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	DUALISO_INFO & dualiso_info,
	bool printInfo
	)
{
	//if (printInfo)
	//{
	//	cout <<"d1 "<< d1 <<" sepvert "<< iso_vlist[d1].sep_vert<< " degree "<<iso_vlist[d1].ver_degree <<endl;
	//	cout <<"d2 "<< d2 <<" sepvert "<< iso_vlist[d2].sep_vert<< " degree "<<iso_vlist[d2].ver_degree <<endl;
	//}
	if (iso_vlist[d1].sep_vert == iso_vlist[d2].sep_vert
		&& iso_vlist[d1].ver_degree == 3 && iso_vlist[d2].ver_degree == 3)
	{ 

		if(printInfo)
			cout <<"capQuad"<<endl;

		//CAP QUAD 
		float closestDistance=0;
		IJK::ARRAY<GRID_COORD_TYPE> d2EdgeBase(DIM3,0);
		IJK::ARRAY<GRID_COORD_TYPE> d1EdgeBase(DIM3,0);

		find_closest_edge_and_distance(DIM3, &(d2Coord[0]), edgeDir, d2EdgeBase, closestDistance);
		/*if ( printInfo)
		{
		cout <<"close distance "<< closestDistance <<endl;
		cout <<scalar_grid.ComputeVertexIndex(&d2EdgeBase[0])
		<<" == "<< scalar_grid.ComputeVertexIndex(&edgeBase[0])<<endl;
		cout <<"edgebase "<< edgeBase[0]<<" "<<edgeBase[1]<<" "<<edgeBase[2]<<endl;
		cout <<"d2edgeBase "<< d2EdgeBase[0]<<" "<<d2EdgeBase[1]<<" "<<d2EdgeBase[2]<<endl;
		cout <<"d2coord " << d2Coord[0]<<" "<<d2Coord[1]<<" "<<d2Coord[2]<<endl;

		}*/
		if(closestDistance < epsilon)
		{
			if (scalar_grid.ComputeVertexIndex(&d2EdgeBase[0])
				== scalar_grid.ComputeVertexIndex(&edgeBase[0]))
			{
				if (is_permitted_collapse(iso_vlist, d2,  &(d2EdgeBase[0]),
					edgeDir, EDGE))
				{
					if ( printInfo)
					{
						cout <<"*** CapQuad **** "
							<<"quad index "<<quad_vert[VERT_PER_QUAD*q] <<" coord "
							<<vertex_coord[DIM3*quad_vert[VERT_PER_QUAD*q]]
						<<" "<<vertex_coord[DIM3*quad_vert[VERT_PER_QUAD*q]+1]
						<<" "<<vertex_coord[DIM3*quad_vert[VERT_PER_QUAD*q]+2]

						<<"\nquad index "<<quad_vert[VERT_PER_QUAD*q+1] <<" coord "
							<<vertex_coord[DIM3*quad_vert[VERT_PER_QUAD*q+1]]
						<<" "<<vertex_coord[DIM3*quad_vert[VERT_PER_QUAD*q+1]+1]
						<<" "<<vertex_coord[DIM3*quad_vert[VERT_PER_QUAD*q+1]+2]

						<<"\nquad index "<<quad_vert[VERT_PER_QUAD*q+2] <<" coord "
							<<vertex_coord[DIM3*quad_vert[VERT_PER_QUAD*q+2]]
						<<" "<<vertex_coord[DIM3*quad_vert[VERT_PER_QUAD*q+2]+1]
						<<" "<<vertex_coord[DIM3*quad_vert[VERT_PER_QUAD*q+2]+2]



						<<"\nquad index "<<quad_vert[VERT_PER_QUAD*q+3] <<" coord "
							<<vertex_coord[DIM3*quad_vert[VERT_PER_QUAD*q+3]]
						<<" "<<vertex_coord[DIM3*quad_vert[VERT_PER_QUAD*q+3]+1]
						<<" "<<vertex_coord[DIM3*quad_vert[VERT_PER_QUAD*q+3]+2]<< endl;

						cout <<"Edgebase "
							<<edgeBase[0]<<" "<<edgeBase[1]<<" "<<edgeBase[2]<<endl;
						cout <<"Edge dir "<< edgeDir <<endl;

					}
					collapseToEdge(d2, q, quad_vert, vertex_coord, iso_vlist
						,collapse_map, edgeBase, edgeDir, dualiso_info);

					dualiso_info.cp_info.moved2Edge++;
				}
			}
		}
		else
		{
			find_closest_edge_and_distance(DIM3, &(d1Coord[0]), edgeDir, d1EdgeBase, closestDistance);
			/*if ( printInfo)
			{
			cout <<"close distance "<< closestDistance <<endl;
			cout <<scalar_grid.ComputeVertexIndex(&d2EdgeBase[0])
			<<" == "<< scalar_grid.ComputeVertexIndex(&edgeBase[0])<<endl;

			}*/
			if (closestDistance <   epsilon)
			{
				if (scalar_grid.ComputeVertexIndex(&d1EdgeBase[0])
					== scalar_grid.ComputeVertexIndex(&edgeBase[0]))
				{
					if (is_permitted_collapse(iso_vlist, d1,  &(d1EdgeBase[0]),
						edgeDir, EDGE))
					{
						if ( printInfo)
						{
							cout <<"*** CapQuad *** "
								<<quad_vert[VERT_PER_QUAD*q] <<" "
								<<quad_vert[VERT_PER_QUAD*q+1] <<" "
								<<quad_vert[VERT_PER_QUAD*q+2] <<" "
								<<quad_vert[VERT_PER_QUAD*q+3] << endl;
							cout <<"edgebase "
								<<edgeBase[0]<<" "<<edgeBase[1]<<" "<<edgeBase[2]<<endl;
							cout <<"dir "<< edgeDir <<endl;
						}
						collapseToEdge(d1, q, quad_vert, vertex_coord, iso_vlist
							,collapse_map, edgeBase, edgeDir, dualiso_info);

						dualiso_info.cp_info.moved2Edge++;
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


		if (printInfo)
		{
			cout <<"q "<< q<<" quad index          "<< A <<" "<<B<<" "<<C<<" "<<D <<endl; 
		}
		//update vertices
		A = find_vertex(collapse_map,A);
		B = find_vertex(collapse_map,B);
		C = find_vertex(collapse_map,C);
		D = find_vertex(collapse_map,D);
		bool flagCapQuad = false;

		//if (printInfo)
		//{
		//	cout <<"q "<< q<<" quad indexUpdated to "<< A <<" "<<B<<" "<<C<<" "<<D <<endl; 
		//}

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

		//if (printInfo)
		//{
		//	cout <<"edgebase "<< edgeBase[0]<<" "<<edgeBase[1]<<" "<<edgeBase[2]<<" edgeDir "<< edgeDir <<endl;

		//}
		capQuadComputation(A,B, ACoord, BCoord, q, quad_vert, vertex_coord,
			scalar_grid, epsilon, notDegreeThree, iso_vlist, 
			edgeBase, edgeDir, flagCapQuad, collapse_map, dualiso_info, printInfo);

		IJK::ARRAY<COORD_TYPE> CCoord(DIM3,0);
		IJK::ARRAY<COORD_TYPE> DCoord(DIM3,0);

		for (int d=0;d<DIM3;d++)
		{
			DCoord[d]=vertex_coord[DIM3*C+d];
			CCoord[d]=vertex_coord[DIM3*D+d];
		}
		capQuadComputation(C,D, CCoord,DCoord, q, quad_vert, vertex_coord,
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
// RETURNS:
// *The collpase map.
void QCOLLAPSE::dual_collapse(
	const DUALISO_DATA & dualiso_data,
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	std::vector<VERTEX_INDEX> & quad_vert,
	IJKDUALTABLE::ISODUAL_CUBE_TABLE &isodual_table,
	std::vector<DUAL_ISOVERT> & iso_vlist, 
	std::vector<COORD_TYPE> & vertex_coord,
	const std::vector<DIRECTION_TYPE> &orth_dir,
	const float epsilon,
	IJK::ARRAY<VERTEX_INDEX> &collapse_map,
	DUALISO_INFO & dualiso_info
	)
{
	//Reordering QuadVert 
	IJK::reorder_quad_vertices(quad_vert);
	clock_t t0,t1,t2,t3,t4;
	t0 = clock();

	//CYCLIC ORDER
	if (dualiso_data.flag_cap_col)
	{
		collapse_caps(scalar_grid, quad_vert, vertex_coord, orth_dir, 
			epsilon, iso_vlist, collapse_map, dualiso_info, dualiso_data.flag_collapse_debug);
	}

	if (dualiso_data.flag_use_collapse_B)
	{
		t1=clock();
		collapse_across_facetsB(scalar_grid, isodual_table, iso_vlist, quad_vert,
			vertex_coord, orth_dir, collapse_map, epsilon,
			dualiso_data.flag_collapse_debug, dualiso_info);
		t2=clock();
		collapse_across_edgesB(scalar_grid, isodual_table, iso_vlist, quad_vert,
			vertex_coord, collapse_map, epsilon, dualiso_data.flag_collapse_debug, dualiso_info);
		t3=clock();

		collapse_across_verticesB(scalar_grid, isodual_table, iso_vlist, quad_vert, 
			vertex_coord, collapse_map, epsilon, dualiso_data.flag_collapse_debug, dualiso_info);
		t4=clock();
	}
	else if (dualiso_data.flag_use_collapse_C)
	{
		t1=clock();
		collapse_across_facetsC(scalar_grid, isodual_table, iso_vlist, quad_vert,
			vertex_coord, orth_dir, collapse_map, epsilon,
			dualiso_data.flag_collapse_debug, dualiso_info);
		t2=clock();
		collapse_across_edgesB(scalar_grid, isodual_table, iso_vlist, quad_vert,
			vertex_coord, collapse_map, epsilon, dualiso_data.flag_collapse_debug, dualiso_info);
		t3=clock();

		collapse_across_verticesB(scalar_grid, isodual_table, iso_vlist, quad_vert, 
			vertex_coord, collapse_map, epsilon, dualiso_data.flag_collapse_debug, dualiso_info);
		t4=clock();
	}
	else
	{
		t1=clock();
		collapse_across_facets(scalar_grid, isodual_table, iso_vlist, quad_vert,
			vertex_coord, collapse_map, epsilon, dualiso_data.flag_collapse_debug, dualiso_info);

		t2=clock();

		collapse_across_edges(scalar_grid, isodual_table, iso_vlist, quad_vert,
			vertex_coord, collapse_map, epsilon, dualiso_data.flag_collapse_debug, dualiso_info);
		t3=clock();

		collapse_across_vertices(scalar_grid, isodual_table, iso_vlist, quad_vert, 
			vertex_coord, collapse_map, epsilon, dualiso_data.flag_collapse_debug, dualiso_info);
		t4=clock();
	}
	update_quads(collapse_map, quad_vert);
	//reorder the quads back to the original.
	IJK::reorder_quad_vertices(quad_vert);

	IJK::clock2seconds(t1-t0, dualiso_info.time.collapse_caps);
	IJK::clock2seconds(t2-t1, dualiso_info.time.collapse_across_facets);
	IJK::clock2seconds(t3-t2, dualiso_info.time.collapse_across_edges);
	IJK::clock2seconds(t4-t3, dualiso_info.time.collapse_across_vertices);

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


