#include "qdual_remove_degen.h"
#include "qdual_collapse.h"


#include <cmath>

using namespace std;
using namespace IJK;
using namespace NamedConstants;
using namespace QCOLLAPSE;
//Returns true if degen quads exists
//Returns the quadvert without the degen QUADS.
//Returns the track_quad_indices. 
bool QTRIANGULATE::triangulate_degenerate_quads (
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<VERTEX_INDEX> & tri_vert,
	const std::vector<COORD_TYPE> & vertex_coord,
	unordered_map<QUAD_INDEX, QUAD_INDEX> &track_quad_indices)
{

	std::vector<VERTEX_INDEX> non_degen_quad_vert;
	remove_degenerate_quads(quad_vert, non_degen_quad_vert,
		tri_vert, vertex_coord, track_quad_indices);
	//reset the quad vert
	quad_vert.clear();
	quad_vert = non_degen_quad_vert;
	if (tri_vert.size() == 0 )
		return false;
	else
		return true;
}


bool sharesOppositeVerticesComputations
	(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const QUAD_INDEX q, //qth quad
	QUAD_INDEX & q2,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
	std::vector<VERTEX_INDEX> & quad_vert,
	const std::vector<VERTEX_INDEX> & origQuadVert,
	std::unordered_map<VERTEX_INDEX,VERTEX_INDEX> & diagonalMap,
	const std::vector<DIRECTION_TYPE> &orth_dir,
	vector< VERTEX_INDEX> & v,
	unordered_map<QUAD_INDEX, QUAD_INDEX > & track_quad_indices,
	IJK::ARRAY<VERTEX_INDEX> & collapse_map
	)
{
	if (diagonalMap.empty())
		return false;

	v.clear();
	unordered_map<QUAD_INDEX, QUAD_INDEX>::const_iterator it = 
		track_quad_indices.find(q);
	QUAD_INDEX mappedQ = it->second;

	VERTEX_INDEX B = iso_vlist [origQuadVert[VERT_PER_QUAD*mappedQ+2]].cube_index;
	VERTEX_INDEX tempEdgeDir = (int)orth_dir[q];
	VERTEX_INDEX ie = DIM3*B + tempEdgeDir;

	VERTEX_INDEX edgeDir = ie%DIM3;
	VERTEX_INDEX iend = (ie-edgeDir)/DIM3;

	IJK::ARRAY<GRID_COORD_TYPE> edgeBase(DIM3,0);
	scalar_grid.ComputeCoord(iend, &(edgeBase[0]));

	if (edgeBase[edgeDir]!=0)
	{
		VERTEX_INDEX iprev = scalar_grid.PrevVertex(iend, edgeDir );
		VERTEX_INDEX ie2 = iprev*DIM3 + edgeDir;
		std::unordered_map<VERTEX_INDEX, VERTEX_INDEX>::iterator ind 
			= diagonalMap.find(ie2);
		if (ind  != diagonalMap.end())
		{
			q2 = diagonalMap.at(ie2);
			if ((
				(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD])) == 
				find_vertex(collapse_map, (origQuadVert[q2*VERT_PER_QUAD]))
				)&&( 
				find_vertex(collapse_map, quad_vert[q*VERT_PER_QUAD+2])
				== find_vertex(collapse_map, origQuadVert[q2*VERT_PER_QUAD+2])
				))
			{
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD]));
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+2]));
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+1]));
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+3]));
				return true;		
			}
			if ((
				(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+1])) == 
				find_vertex(collapse_map, (origQuadVert[q2*VERT_PER_QUAD+3]))
				)&&( 
				find_vertex(collapse_map, quad_vert[q*VERT_PER_QUAD+3])
				== find_vertex(collapse_map, origQuadVert[q2*VERT_PER_QUAD+1])
				))
			{
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+1]));
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+3]));
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD]));
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+2]));
				return true;		
			}
		}
	}
	const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize(); 
	if (edgeBase[edgeDir] != (axis_size[edgeDir]-1))
	{
		VERTEX_INDEX inext = scalar_grid.NextVertex(iend, edgeDir);
		VERTEX_INDEX ie2 = inext*DIM3 + edgeDir;

		std::unordered_map<VERTEX_INDEX, VERTEX_INDEX>::iterator ind 
			= diagonalMap.find(ie2);
		if (ind != diagonalMap.end())
		{
			q2 = diagonalMap.at(ie2);
			if ((
				(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD])) == 
				find_vertex(collapse_map, (origQuadVert[q2*VERT_PER_QUAD]))
				)&&( 
				find_vertex(collapse_map, quad_vert[q*VERT_PER_QUAD+2])
				== find_vertex(collapse_map, origQuadVert[q2*VERT_PER_QUAD+2])
				))
			{
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD]));
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+2]));
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+1]));
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+3]));
				return true;		
			}
			if ((
				(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+1])) == 
				find_vertex(collapse_map, (origQuadVert[q2*VERT_PER_QUAD+1]))
				)&&( 
				find_vertex(collapse_map, quad_vert[q*VERT_PER_QUAD+3])
				== find_vertex(collapse_map, origQuadVert[q2*VERT_PER_QUAD+3])
				))
			{
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+1]));
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+3]));
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD]));
				v.push_back(find_vertex( collapse_map, quad_vert[q*VERT_PER_QUAD+2]));
				return true;		
			}
		}
	}
	return false;
}

//Setup the hashmap to find quads which share diagonals
void QTRIANGULATE::hashQuadsDual2GridEdge(
	std::unordered_map<VERTEX_INDEX,VERTEX_INDEX> & diagonalMap,
	std::vector<VERTEX_INDEX> & quad_vert,
	const std::vector<DIRECTION_TYPE> &orth_dir,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	const std::vector<COORD_TYPE> & vertex_coord)
{
	const int numQuads = quad_vert.size() / VERT_PER_QUAD;
	IJK::reorder_quad_vertices(quad_vert);
	for (int q = 0; q < numQuads; q++)
	{
		VERTEX_INDEX B = iso_vlist [quad_vert[VERT_PER_QUAD*q+2]].cube_index;
		int edgeDir = (int)orth_dir[q];
		int index = DIM3*B + edgeDir;
		diagonalMap.insert(make_pair(index, q));
		//debug
		/*
		cout <<"q "<< q << endl;
		cout <<"B(cube index) "<< B<<" ";
		VERTEX_INDEX temp = quad_vert[VERT_PER_QUAD*q+2];
		cout <<" issovertex index [" << temp <<"]  coord [";
		cout <<vertex_coord[DIM3*temp]<<" ";
		cout <<vertex_coord[DIM3*temp+1]<<" ";
		cout <<vertex_coord[DIM3*temp+2]<<"] "<<endl;
		cout <<"quad vert  "<< quad_vert[VERT_PER_QUAD*q]<<" : " <<iso_vlist [quad_vert[VERT_PER_QUAD*q]].cube_coord[0]<<" "
		<<iso_vlist [quad_vert[VERT_PER_QUAD*q]].cube_coord[1]<<" "
		<<iso_vlist [quad_vert[VERT_PER_QUAD*q]].cube_coord[2]<<
		" cubeindex "<< iso_vlist [quad_vert[VERT_PER_QUAD*q]].cube_index <<endl;

		cout <<"quad vert  "<< quad_vert[VERT_PER_QUAD*q+1]<<" : " <<iso_vlist [quad_vert[VERT_PER_QUAD*q+1]].cube_coord[0]<<" "
		<<iso_vlist [quad_vert[VERT_PER_QUAD*q+1]].cube_coord[1]<<" "
		<<iso_vlist [quad_vert[VERT_PER_QUAD*q+1]].cube_coord[2]<<
		" cubeindex "<< iso_vlist [quad_vert[VERT_PER_QUAD*q+1]].cube_index <<endl;


		cout <<"quad vert  "<< quad_vert[VERT_PER_QUAD*q+2]<<" : " <<iso_vlist [quad_vert[VERT_PER_QUAD*q+2]].cube_coord[0]<<" "
		<<iso_vlist [quad_vert[VERT_PER_QUAD*q+2]].cube_coord[1]<<" "
		<<iso_vlist [quad_vert[VERT_PER_QUAD*q+2]].cube_coord[2]<<
		" cubeindex "<< iso_vlist [quad_vert[VERT_PER_QUAD*q+2]].cube_index <<endl;


		cout <<"quad vert "<< quad_vert[VERT_PER_QUAD*q+3]<<" : " <<iso_vlist [quad_vert[VERT_PER_QUAD*q+3]].cube_coord[0]<<" "
		<<iso_vlist [quad_vert[VERT_PER_QUAD*q+3]].cube_coord[1]<<" "
		<<iso_vlist [quad_vert[VERT_PER_QUAD*q+3]].cube_coord[2]<<
		" cubeindex "<< iso_vlist [quad_vert[VERT_PER_QUAD*q+3]].cube_index <<endl;


		cout <<"other vertices "<<endl;
		cout <<quad_vert[VERT_PER_QUAD*q] <<" "<<quad_vert[VERT_PER_QUAD*q+1]<<" "
		<<quad_vert[VERT_PER_QUAD*q+3]<<endl;
		cout <<"edgeDir "<< edgeDir <<endl;
		cout <<"index "<< index <<endl;
		*/

	}
	IJK::reorder_quad_vertices(quad_vert);

}

// Remove degenerate quads
void remove_degen (
	const int q, vector<VERTEX_INDEX> & quad_vert,
	vector<VERTEX_INDEX> & non_degen_vert,
	int & num_non_degen)
{
	num_non_degen=0;
	int i0 = q*4;
	if (quad_vert[i0]==quad_vert[i0+2])
		return;
	if (quad_vert[i0+1]==quad_vert[i0+3])
		return;
	int v1, v2;
	const int count=4;
	for ( v1 = count-1, v2 = 0; v2 < count; v1 = v2++ ) 
	{
		int endPt1 = quad_vert[q*4+v1];
		int endPt2 = quad_vert[q*4+v2];
		if (endPt1!=endPt2)
		{
			non_degen_vert.push_back(endPt2);
			num_non_degen++;
		}
	}
	if (num_non_degen <=2)
		num_non_degen=0;
}

// Remove degenerate quads
void QTRIANGULATE::remove_degenerate_quads(
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<VERTEX_INDEX> & non_degen_quad_vert, // only non degenerate quads
	std::vector<VERTEX_INDEX> & tri_vert,
	const std::vector<COORD_TYPE> & vertex_coord,
	unordered_map<QUAD_INDEX, QUAD_INDEX > &track_quad_indices
	)
{

	//set up the vertex collapse map
	int num_vertex = vertex_coord.size();
	//Reordering QuadVert 
	IJK::reorder_quad_vertices(quad_vert);
	const int num_quads = quad_vert.size()/VERT_PER_QUAD;

	QUAD_INDEX q2 = 0; // number of non-degen quad
	for (int q=0; q < num_quads; q++)
	{
		vector<VERTEX_INDEX> non_degen_verts;
		int num_non_degen = 0;
		remove_degen(q, quad_vert, non_degen_verts, num_non_degen);

		if ( num_non_degen == 0 )
			continue;
		else if (num_non_degen == VERT_PER_QUAD )
		{
			for (int d=0;d<VERT_PER_QUAD;d++)
			{
				non_degen_quad_vert.push_back(quad_vert[4*q+d]);
			}
			track_quad_indices.insert(make_pair(q2,q));
			q2++;
		}
		else // triangle
		{
			for (int k = 0; k < num_non_degen; k++)
			{
				tri_vert.push_back(non_degen_verts[k]);
			}
		}
	}
	//reorder the quads back to the original.
	IJK::reorder_quad_vertices(quad_vert);
	IJK::reorder_quad_vertices(non_degen_quad_vert);
}



// Compute degree of each vertex
// Only for non degenerate poly
// param 1 : num vertex in the poly.
// param 2 : non degenerate polys
void QTRIANGULATE::compute_degree_per_vertex(
	const int vert_per_poly,
	std::vector<VERTEX_INDEX> & poly_vert, // only non degenerate quads
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist		
	)
{
	const int num_poly = poly_vert.size()/vert_per_poly;
	for (int q=0; q<num_poly; q++)
	{
		for (int d=0;d<vert_per_poly;d++){
			VERTEX_INDEX v = poly_vert[vert_per_poly*q+d];
			iso_vlist[v].ver_degree++;
		}
	}
}

void QTRIANGULATE::reset_degree_per_vertex(
	const int vert_per_poly,
	std::vector<VERTEX_INDEX> & poly_vert, // only non degenerate quads
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist		
	)
{
	const int num_poly = poly_vert.size()/vert_per_poly;
	for (int q=0; q<num_poly; q++)
	{
		for (int d=0;d<vert_per_poly;d++){
			VERTEX_INDEX v = poly_vert[vert_per_poly*q+d];
			iso_vlist[v].ver_degree=0;
		}
	}
}

//Compute the length of the side of a triangle 
//given the vertex index of the endpts and the vertexCoord
float length_of_side (	const VERTEX_INDEX A, const VERTEX_INDEX B,
					  const std::vector<COORD_TYPE> & vertex_coord)
{

	return 	sqrt( (vertex_coord[3*A]-vertex_coord[3*B])*(vertex_coord[3*A]-vertex_coord[3*B])
		+ (vertex_coord[3*A+1]-vertex_coord[3*B+1])*(vertex_coord[3*A+1]-vertex_coord[3*B+1]) 
		+ (vertex_coord[3*A+2]-vertex_coord[3*B+2])*(vertex_coord[3*A+2]-vertex_coord[3*B+2]) );
}


//compute the cosine of the angle subtended by the FIRST PARAM.
//cosine rule
//c2=a2+b2-2abcosC
float compute_cosine_rule_angle(const float A, const float B, const float C)
{
	return ((B*B + C*C - A*A)/ (2.0*B*C));
}
//compute the cosine of the min angle in the triangle

void compute_cos_min_angle (
	const VERTEX_INDEX A, const VERTEX_INDEX B, const VERTEX_INDEX C,
	const std::vector<COORD_TYPE> & vertex_coord,
	float &cos_angle)
{
	float AB = length_of_side (A,B,vertex_coord);
	float BC = length_of_side (B, C, vertex_coord);
	float CA = length_of_side (C, A, vertex_coord);


	float cos_angle_C = compute_cosine_rule_angle(AB, BC, CA);
	float cos_angle_B = compute_cosine_rule_angle(BC, CA, AB);
	float cos_angle_A = compute_cosine_rule_angle(CA, BC, AB);

	if ((cos_angle_A > cos_angle_B) && (cos_angle_A > cos_angle_C))
		cos_angle = cos_angle_A;
	else if ((cos_angle_B > cos_angle_A)&&(cos_angle_B > cos_angle_C))
		cos_angle = cos_angle_B;
	else
		cos_angle = cos_angle_C;
}

// update collapses in trivert with collapse map
void update_tris(
	const IJK::ARRAY<VERTEX_INDEX> & collapse_map,
	std::vector<VERTEX_INDEX> & tri_vert
	)
{
	const int tri_vert_size = tri_vert.size();
	for (int v=0; v<tri_vert_size; v++)
	{
		int k = tri_vert[v];
		tri_vert[v]=collapse_map[k];
	}
}

// push triangle ABC into trivert
void push_triangle(
	const int A, 
	const int B, 
	const int C,
	std::vector<VERTEX_INDEX> & tri_vert)
{
	tri_vert.push_back(C);
	tri_vert.push_back(B);
	tri_vert.push_back(A);
}

//Remove the triangle duplciates
void remove_tri_degenerates(
	std::vector<VERTEX_INDEX> & tri_vert,
	std::vector<VERTEX_INDEX> & tri_vert_non_dupli)
{
	const int size = tri_vert.size()/DIM3;

	for (int i = 0; i < size; i++)
	{
		bool flag_degen = false;
		VERTEX_INDEX v1,v2=0;
		for (  v1 =  DIM3-1,  v2 = 0; v2 <DIM3; v1 = v2++ ) 
		{
			if (tri_vert[DIM3*i+v1] == tri_vert[DIM3*i+v2])
			{
				flag_degen = true;
				break;
			}
		}
		if (!flag_degen)
		{
			tri_vert_non_dupli.push_back(tri_vert[DIM3*i]);
			tri_vert_non_dupli.push_back(tri_vert[DIM3*i+1]);
			tri_vert_non_dupli.push_back(tri_vert[DIM3*i+2]);
		}
	}
}

//Check if the vertex is a boundary 
//param 1 index into isovlist
//returns flag_boundary true if in the boundary
void  QTRIANGULATE::isBoundaryIsoVertex(
	const int vertex, //index into isovlist
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist,
	IJK::BOOL_GRID<DUALISO_GRID> &boundary_grid,
	QDUAL_TABLE & qdual_table,
	bool & flag_boundary
	)
{
	flag_boundary = false;
	const VERTEX_INDEX v = iso_vlist[vertex].cube_index;
	if (!boundary_grid.Scalar(v))
	{
		flag_boundary = false;
		return;
	}
	for (int i = 0; i < NUM_CUBE_FACETS; i++)
	{
		bool flagFacetIsBoundary = false;
		int orthDir = i%DIM3;
		if (i<DIM3)
		{
			flagFacetIsBoundary = boundary_grid.IsCubeFacetOnGridBoundary(v, orthDir, false);
		}
		else
		{
			flagFacetIsBoundary = boundary_grid.IsCubeFacetOnGridBoundary(v, orthDir, true);
		}
		if (flagFacetIsBoundary)
		{
			int it = iso_vlist[vertex].table_index;
			unsigned char ip = iso_vlist[vertex].patch_index;
			QDUAL_TABLE::DIR_BITS w_edge_flag = qdual_table.connectDir(it,ip);
			if ((w_edge_flag & (1<<i))!=0)
			{
				flag_boundary = true;
			}
		}
	}
}


// Triangulate quads based on their angles
// Param 1: set of NON_DEGENERATE quads
void QTRIANGULATE::triangulate_quad_angle_based(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	std::vector<VERTEX_INDEX> & quadVert, // only non degenerate quads
	const std::vector<VERTEX_INDEX> & origQuadVert,
	std::vector<VERTEX_INDEX> & tri_vert,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
	const std::vector<COORD_TYPE> & vertex_coord,
	IJK::BOOL_GRID<DUALISO_GRID> &boundary_grid,
	QDUAL_TABLE & qdual_table,
	std::unordered_map<VERTEX_INDEX,VERTEX_INDEX> & diagonalMap,
	const std::vector<DIRECTION_TYPE> &orth_dir,
	unordered_map< QUAD_INDEX, QUAD_INDEX > & track_quad_indices,
	IJK::ARRAY<VERTEX_INDEX> & origCollapse_map,
	const bool printInfo
	)
{
	const int num_quad = quadVert.size()/VERT_PER_QUAD;
	IJK::reorder_quad_vertices(quadVert);
	compute_degree_per_vertex(4, quadVert, iso_vlist);
	compute_degree_per_vertex(3, tri_vert, iso_vlist);

	//set up the vertex collapse map
	int num_vertex = vertex_coord.size();
	IJK::ARRAY<VERTEX_INDEX> collapse_map(num_vertex,0);
	//setup collapse_map.
	QCOLLAPSE::setup_collapse_map(collapse_map, num_vertex);

	vector <QUAD_INDEX> quadIndices;
	for (int q=0; q < num_quad; q++)
	{
		int w1=0;// has the index NOT the vertex
		bool flag_deg2 = false;
		bool flag_boundary = false;
		int vertex;
		int num_deg2=0;

		for (int j=0; j<VERT_PER_QUAD; j++)
		{
			vertex = QCOLLAPSE::find_vertex(collapse_map,quadVert[q*VERT_PER_QUAD+j]);
			if(printInfo){
				cout <<"q "<< q<<" j "<<j<< " " <<quadVert[q*VERT_PER_QUAD+j];
				cout <<" mappped to "<< vertex << " degree "<< iso_vlist[vertex].ver_degree <<endl;
				cout <<"cube_coord "<<iso_vlist[vertex].cube_coord[0]<<" "<<iso_vlist[vertex].cube_coord[1]
				<<" "<<iso_vlist[vertex].cube_coord[2]<<endl;
			}

			if (iso_vlist[vertex].ver_degree == 2)
			{
				w1 = j;
				num_deg2++;
				isBoundaryIsoVertex(vertex, iso_vlist, boundary_grid,
					qdual_table, flag_boundary);
				flag_deg2 = true;
				//break;
				if(printInfo){
					cout <<"vertex "<<vertex <<" boundary "<< flag_boundary <<endl;
				}
			}
		}

		bool flag_non_degen_quad = false;

		int num_non_degen = 0;
		VERTEX_INDEX temp_tri[3]={0};
		VERTEX_INDEX v1,v2=0;
		for (  v1 =  VERT_PER_QUAD-1,  v2 = 0; v2 <VERT_PER_QUAD; v1 = v2++ ) 
		{
			int endPt1 = QCOLLAPSE::find_vertex(collapse_map,quadVert[q*4+v1]);
			int endPt2 = QCOLLAPSE::find_vertex(collapse_map,quadVert[q*4+v2]);
			if (endPt1!=endPt2)
			{
				temp_tri[num_non_degen]=endPt2;
				num_non_degen++;
			}
		}
		if (num_non_degen <=2)
		{
			num_non_degen=0;
		}

		if ( num_non_degen == 0 )
			continue;
		else if (num_non_degen == 4 )
		{
			// do computation
			flag_non_degen_quad = true;
		}
		else // triangle
		{
			for (int d=0; d< num_non_degen; d++)
			{
				tri_vert.push_back(temp_tri[d]);
			}
		}
		//
		if (flag_non_degen_quad)
		{
			vertex = quadVert[q*VERT_PER_QUAD+ (w1)%VERT_PER_QUAD];	
			int B = quadVert[q*VERT_PER_QUAD+ (w1+2)%VERT_PER_QUAD];
			int C = quadVert[q*VERT_PER_QUAD+ (w1+3)%VERT_PER_QUAD];
			int D = quadVert[q*VERT_PER_QUAD+ (w1+1)%VERT_PER_QUAD];


			vertex = QCOLLAPSE::find_vertex(collapse_map,vertex);
			B = QCOLLAPSE::find_vertex(collapse_map,B);
			C = QCOLLAPSE::find_vertex(collapse_map,C);
			D = QCOLLAPSE::find_vertex(collapse_map, D);
			//find the 4 angles.

			float cos_min_angle_abc; // min angle of triangle ABC.
			compute_cos_min_angle(vertex, B, C, vertex_coord, cos_min_angle_abc);
			float cos_min_angle_abd; // min angle of triangle ABD.
			compute_cos_min_angle(vertex, B, D, vertex_coord, cos_min_angle_abd);

			float cos_min_angle_dbc; // min angle of triangle DBC.
			compute_cos_min_angle(D, B, C, vertex_coord, cos_min_angle_dbc);
			float cos_min_angle_dca; // min angle of triangle DCA.
			compute_cos_min_angle(vertex, D, C, vertex_coord, cos_min_angle_dca);

			// min of the cosines of the first set of angles ABC and ABD 
			float cos_min_1;
			if (cos_min_angle_abc > cos_min_angle_abd )
				cos_min_1 = cos_min_angle_abc;
			else
				cos_min_1 = cos_min_angle_abd;

			//min of the cosine of the second set of angles DBC and DCA
			float cos_min_2;
			if (cos_min_angle_dbc > cos_min_angle_dca)
				cos_min_2 = cos_min_angle_dbc;
			else
				cos_min_2 = cos_min_angle_dca;


			//degree 2 and not boundary   
			if (flag_deg2 && !flag_boundary)
			{
				if(printInfo){
					cout <<"flag_deg2 && !flag_boundary "<< endl;
				}
				if (cos_min_1 < cos_min_2)
				{
					push_triangle(C, B, vertex, tri_vert);
					push_triangle(vertex, B, D, tri_vert);
				}
				else
				{
					if (length_of_side(vertex, C, vertex_coord) <
						length_of_side(vertex, D, vertex_coord))
					{
						// MAP vertex to C.
						collapse_map[vertex] = C;
						if(printInfo){
							cout <<vertex<<" mapped2 to " << C <<endl;
						}
					}
					else
					{
						// MAP vertex to D
						collapse_map[vertex] = D;
						if(printInfo){
							cout <<vertex<<" mapped2 to " << D <<endl;
						}
					}

					push_triangle(D, C, B, tri_vert);
				}

			}
			//shares opposite vertices

			else
			{
				QUAD_INDEX commonQuad;
				if ( sharesOppositeVerticesComputations
					(scalar_grid, q, commonQuad, iso_vlist, quadVert,
					origQuadVert, diagonalMap,
					orth_dir, quadIndices, track_quad_indices, origCollapse_map) )
				{
					push_triangle( quadIndices[0], quadIndices[2], quadIndices[3], tri_vert);
					push_triangle(quadIndices[3], quadIndices[2], quadIndices[1], tri_vert);
					if(printInfo){
						cout <<"diagonal "<<endl;
					}
				}
				else
				{ // degree 2
					if(printInfo){
						cout <<"generic angle based"<<endl;
					}
					if (cos_min_1 < cos_min_2)
					{
						push_triangle(vertex, C, B, tri_vert);
						push_triangle(vertex, B, D, tri_vert);
						if(printInfo){
							cout <<vertex<<" "<<C<<" "<<B<<endl;
							cout <<vertex<<" "<<B<<" "<<D<<endl;
						}
					}
					else
					{
						push_triangle(D, C, B, tri_vert);
						push_triangle(vertex, C, D, tri_vert);
						if(printInfo){
							cout <<D<<" "<<C<<" "<<B<<endl;
							cout <<vertex<<" "<<C<<" "<<D<<endl;
						}
					}
				}
			}
		}//if	
	}

	update_tris(collapse_map, tri_vert);
	vector<VERTEX_INDEX> temp_tri_vert;
	remove_tri_degenerates (tri_vert, temp_tri_vert);
	tri_vert.clear();
	tri_vert = temp_tri_vert;
	//reorder back to original 
	IJK::reorder_quad_vertices(quadVert);

}
