#include "qdual_remove_degen.h"
#include "qdual_collapse.h"


#include <cmath>

using namespace std;
using namespace IJK;
using namespace NamedConstants;

bool QTRIANGULATE::triangulate_non_degen_quads (
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<VERTEX_INDEX> & tri_vert,
	const std::vector<COORD_TYPE> & vertex_coord)
{

	std::vector<VERTEX_INDEX> non_degen_quad_vert;
	remove_degenerate_quads(quad_vert, non_degen_quad_vert, tri_vert, vertex_coord);
    quad_vert.clear();
	quad_vert = non_degen_quad_vert;
	if (tri_vert.size() == 0 )
		return false;
	else
		return true;
}

//Returns true if the quads param1 and param2 share diagonal vertices
//NOTE: This returns a vector, in the following format
//the first two are the vertex indices of the  quad which are in common and the last two are the other ones.
bool doQuadsShareDiagVertices
	(
	const VERTEX_INDEX q,
	const VERTEX_INDEX q2,
	std::vector<VERTEX_INDEX> & quad_vert,
    vector < VERTEX_INDEX > v 
	)
{
	//if (
	//	(
	//	(quad_vert[q*VERT_PER_QUAD] == quad_vert[q2*VERT_PER_QUAD])
	//	&& (quad_vert[q*VERT_PER_QUAD+2] == quad_vert[q2*VERT_PER_QUAD+2])
	//	)
	//	|| 
	//	(
	//	(quad_vert[q*VERT_PER_QUAD+1] == quad_vert[q2*VERT_PER_QUAD+1])
	//	&& (quad_vert[q*VERT_PER_QUAD+3]==quad_vert[q2*VERT_PER_QUAD+3])
	//	) 
	//	)
	//{
	//	cout <<" They share vertices"<<endl;
	//	return true;
	//}
	//else
 //       return false;

    v.clear();
	if ((quad_vert[q*VERT_PER_QUAD] == quad_vert[q2*VERT_PER_QUAD])
		&& (quad_vert[q*VERT_PER_QUAD+2] == quad_vert[q2*VERT_PER_QUAD+2]
	))
	{
		v.push_back(0);
		v.push_back(2);
		v.push_back(1);
		v.push_back(3);
		return true;
	}

	if(
		(quad_vert[q*VERT_PER_QUAD+1] == quad_vert[q2*VERT_PER_QUAD+1])
		&& (quad_vert[q*VERT_PER_QUAD+3]==quad_vert[q2*VERT_PER_QUAD+3])
		)
	{
		v.push_back(1);
		v.push_back(3);
		v.push_back(0);
		v.push_back(2);
		return true;
	}
    return false;
}

// PreCondition :  quadVert must be reordered before
//NOTE: This returns a vector, in the following format
//the first two are the vertex index of the  quad in common and the last two are the other ones.
bool sharesOppositeVerticesComputations
	(
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	const VERTEX_INDEX q, //qth quad
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
	std::vector<VERTEX_INDEX> & quad_vert,
	std::unordered_map<VERTEX_INDEX,VERTEX_INDEX> & diagonalMap,
	const std::vector<DIRECTION_TYPE> &orth_dir,
	vector< VERTEX_INDEX> & v
	)
{
	if (diagonalMap.empty())
		return false; 

	//VERTEX_INDEX ie = (int)orth_dir[q];
	VERTEX_INDEX B = iso_vlist [quad_vert[VERT_PER_QUAD*q+2]].cube_index;
	VERTEX_INDEX tempEdgeDir = (int)orth_dir[q];
	VERTEX_INDEX ie = DIM3*B + tempEdgeDir;
	
	VERTEX_INDEX edgeDir = ie%DIM3;
	VERTEX_INDEX iend = (ie-edgeDir)/DIM3;
	IJK::ARRAY<GRID_COORD_TYPE> edgeBase(DIM3,0);
	scalar_grid.ComputeCoord(iend, &(edgeBase[0]));
	cout<< "\n checking q "<< q 
		<<" index "<< ie <<endl;
	cout <<"edgebase "<< edgeBase[0]<<" "<< edgeBase[1]<<" "<<edgeBase[2]<<endl;
	cout <<"edgeDir "<< edgeDir	<<" axix size "<<(scalar_grid.AxisSize(edgeDir)-1) <<endl;
	
	if (edgeBase[edgeDir]!=0)
	{
		VERTEX_INDEX iprev = scalar_grid.PrevVertex(iend, edgeDir );
		VERTEX_INDEX ie2 = iprev*DIM3 + edgeDir;
		std::unordered_map<VERTEX_INDEX, VERTEX_INDEX>::iterator ind 
			= diagonalMap.find(ie2);
		cout <<"compare with index1 "<< ie2 <<endl;
		if (ind  != diagonalMap.end())
		{
			VERTEX_INDEX q2 = diagonalMap.at(ie2);

			cout<<" which is quad " << q2 << endl;
			if (doQuadsShareDiagVertices(q, q2, quad_vert,v))
			{
				//debug 
				cout <<"quad q "<< q << " "
					<< "ie "<< ie <<endl;
				cout <<"iend "<<iend <<" ";
				cout <<"coord "<< edgeBase[0]<<" "<<edgeBase[1]<<" "<<edgeBase[2]<<endl;
				cout <<"ie2 "<< ie2 <<endl;
				cout <<"the second quad is "<< q2<<endl;
				return true;
			}
		}
	}
	const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
	cout <<"edgedir "<<edgeDir<<endl; 
	if (edgeBase[edgeDir] != (axis_size[edgeDir]-1))
	{
		VERTEX_INDEX inext = scalar_grid.NextVertex(iend, edgeDir);
		VERTEX_INDEX ie2 = inext*DIM3 + edgeDir;

		std::unordered_map<VERTEX_INDEX, VERTEX_INDEX>::iterator ind 
			= diagonalMap.find(ie2);
		cout <<"compare with index2 "<< ie2 <<endl;
		if (ind  == diagonalMap.end())
			return false;
		else
		{
			VERTEX_INDEX q2 = diagonalMap.at(ie2);
			cout <<" which is quad " << q2 << endl;
			if (doQuadsShareDiagVertices(q, q2, quad_vert,v))
			{
				//debug 
				cout <<"quad q "<< q << " "
					<< "ie "<< ie <<endl;
				cout <<"iend "<<iend <<" ";
				cout <<"coord "<< edgeBase[0]<<" "<<edgeBase[1]<<" "<<edgeBase[2]<<endl;
				cout <<"ie2 "<< ie2 <<endl;
				cout <<"the second quad is "<< q2<<endl;
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

	cout <<"after reordering"<<endl;
	for (int i=0;i<8;i++)
	{
		cout <<" "<< quad_vert[VERT_PER_QUAD*30+i];
	}
	cout <<"\n";
	for (int q = 0; q < numQuads; q++)
	{
		VERTEX_INDEX B = iso_vlist [quad_vert[VERT_PER_QUAD*q+2]].cube_index;
		int edgeDir = (int)orth_dir[q];
		int index = DIM3*B + edgeDir;
		diagonalMap.insert(make_pair(index, q));
	//debug
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

	}
    IJK::reorder_quad_vertices(quad_vert);

}
// Triangulate all quads
void QTRIANGULATE::triangulate_quads (
    const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<VERTEX_INDEX> & tri_vert,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
	const std::vector<COORD_TYPE> & vertex_coord,
	QDUAL_TABLE & qdual_table,
	IJK::BOOL_GRID<DUALISO_GRID> &boundary_grid,
	const std::vector<DIRECTION_TYPE> &orth_dir,
	std::unordered_map<VERTEX_INDEX,VERTEX_INDEX>  &diagonalMap)
{


	bool has_non_degen_quad = triangulate_non_degen_quads(quad_vert, tri_vert, vertex_coord);
	// OBSOLETE
	//std::unordered_map<VERTEX_INDEX,VERTEX_INDEX>  diagonalMap;
	//cout <<"before quad angle based"<<endl;
	//for (int i=0;i<8;i++)
	//{
	//	cout <<" "<< quad_vert[VERT_PER_QUAD*30+i];
	//}
	//cout <<"\n";
	//hashQuadsDual2GridEdge(diagonalMap, quad_vert, orth_dir, iso_vlist, vertex_coord);

	// only non degenerate quads and triangles remain
	triangulate_quad_angle_based(scalar_grid, quad_vert, tri_vert, iso_vlist, 
		vertex_coord, boundary_grid, qdual_table, diagonalMap, orth_dir);
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
	const std::vector<COORD_TYPE> & vertex_coord
	)
{

	//set up the vertex collapse map
	int num_vertex = vertex_coord.size();
	//Reordering QuadVert 
	IJK::reorder_quad_vertices(quad_vert);
	const int num_quads = quad_vert.size()/VERT_PER_QUAD;

	for (int q=0;q<num_quads;q++)
	{
		vector<VERTEX_INDEX> non_degen_verts;
		int num_non_degen = 0;
		remove_degen(q, quad_vert, non_degen_verts, num_non_degen);

		if ( num_non_degen == 0 )
			continue;
		else if (num_non_degen == 4 )
		{
			for (int d=0;d<4;d++)
				non_degen_quad_vert.push_back(quad_vert[4*q+d]);
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
void remove_tri_duplicates(
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
		boundary_grid.IsCubeFacetOnGridBoundary(v, i, flagFacetIsBoundary);
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
	std::vector<VERTEX_INDEX> & non_degen_quad_vert, // only non degenerate quads
	std::vector<VERTEX_INDEX> & tri_vert,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
	const std::vector<COORD_TYPE> & vertex_coord,
	IJK::BOOL_GRID<DUALISO_GRID> &boundary_grid,
	QDUAL_TABLE & qdual_table,
    std::unordered_map<VERTEX_INDEX,VERTEX_INDEX> & diagonalMap,
    const std::vector<DIRECTION_TYPE> &orth_dir
	)
{
	const int num_quad = non_degen_quad_vert.size()/VERT_PER_QUAD;
	IJK::reorder_quad_vertices(non_degen_quad_vert);
	compute_degree_per_vertex(4, non_degen_quad_vert, iso_vlist);
	compute_degree_per_vertex(3, tri_vert, iso_vlist);

	//set up the vertex collapse map
	int num_vertex = vertex_coord.size();
	IJK::ARRAY<VERTEX_INDEX> collapse_map(num_vertex,0);
	//setup collapse_map.
	QCOLLAPSE::setup_collapse_map(collapse_map, num_vertex);


    vector <VERTEX_INDEX> quadIndices;
	for (int q=0; q < num_quad; q++)
	{
		int w1=0;// has the index NOT the vertex
		bool flag_deg2 = false;
		bool flag_boundary = false;
		int vertex;
		int num_deg2=0;
		
		for (int j=0; j<VERT_PER_QUAD; j++)
		{
			vertex = QCOLLAPSE::find_vertex(collapse_map,non_degen_quad_vert[q*VERT_PER_QUAD+j]);
			if (iso_vlist[vertex].ver_degree == 2)
			{
				w1 = j;
				num_deg2++;
				isBoundaryIsoVertex(vertex, iso_vlist, boundary_grid,
					qdual_table, flag_boundary);

				flag_deg2 = true;
				//break;
			}
		}

		bool flag_non_degen_quad = false;


		int num_non_degen = 0;
		VERTEX_INDEX temp_tri[3]={0};
		VERTEX_INDEX v1,v2=0;
		for (  v1 =  VERT_PER_QUAD-1,  v2 = 0; v2 <VERT_PER_QUAD; v1 = v2++ ) 
		{
			int endPt1 = QCOLLAPSE::find_vertex(collapse_map,non_degen_quad_vert[q*4+v1]);
			int endPt2 = QCOLLAPSE::find_vertex(collapse_map,non_degen_quad_vert[q*4+v2]);
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
			vertex = non_degen_quad_vert[q*VERT_PER_QUAD+ (w1)%VERT_PER_QUAD];	
			int B = non_degen_quad_vert[q*VERT_PER_QUAD+ (w1+2)%VERT_PER_QUAD];
			int C = non_degen_quad_vert[q*VERT_PER_QUAD+ (w1+3)%VERT_PER_QUAD];
			int D = non_degen_quad_vert[q*VERT_PER_QUAD+ (w1+1)%VERT_PER_QUAD];
			
			
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
					}
					else
					{
						// MAP vertex to D
						collapse_map[vertex] = D;
					}

					push_triangle(D, C, B, tri_vert);
				}

			}
            //shares opposite vertices
			
			else if ( sharesOppositeVerticesComputations
				(scalar_grid, q, iso_vlist, non_degen_quad_vert, diagonalMap,
                orth_dir, quadIndices) )
			{
				cout <<"test"<<endl;
                cout <<"shares vertices"<<endl;
				int a = QCOLLAPSE::find_vertex(collapse_map, 
					non_degen_quad_vert[q*VERT_PER_QUAD+quadIndices[0]]);
				int b = QCOLLAPSE::find_vertex(collapse_map, 
					non_degen_quad_vert[q*VERT_PER_QUAD+quadIndices[1]]);

				int c = QCOLLAPSE::find_vertex(collapse_map, 
					non_degen_quad_vert[q*VERT_PER_QUAD+quadIndices[2]]);
				int d = QCOLLAPSE::find_vertex(collapse_map, 
					non_degen_quad_vert[q*VERT_PER_QUAD+quadIndices[3]]);

				push_triangle( a, c, d, tri_vert);
				push_triangle(d, c, b, tri_vert);
			}
			else
			{ // degree 2
				if (cos_min_1 < cos_min_2)
				{
					push_triangle(vertex, C, B, tri_vert);
					push_triangle(vertex, B, D, tri_vert);
				}
				else
				{
					push_triangle(D, C, B, tri_vert);
					push_triangle(vertex, C, D, tri_vert);
				}
			}
		}//if	
	}
	
	update_tris(collapse_map, tri_vert);
	vector<VERTEX_INDEX> temp_tri_vert;
	remove_tri_duplicates (tri_vert, temp_tri_vert);
	tri_vert.clear();
	tri_vert = temp_tri_vert;
	//reorder back to original 
	IJK::reorder_quad_vertices(non_degen_quad_vert);

}
