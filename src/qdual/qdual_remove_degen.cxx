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
	quad_vert = non_degen_quad_vert;
	if (tri_vert.size() == 0 )
		return false;
	else
		return true;
}

// Triangulate all quads
void QTRIANGULATE::triangulate_quads (
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<VERTEX_INDEX> & tri_vert,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
	const std::vector<COORD_TYPE> & vertex_coord)
{
	bool has_non_degen_quad = triangulate_non_degen_quads(quad_vert, tri_vert, vertex_coord);
	// only non degenerate quads and triangles remain
	triangulate_quad_angle_based(quad_vert, tri_vert, iso_vlist, vertex_coord);
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
			for (int d=0;d<num_non_degen;d++)
				tri_vert.push_back(non_degen_verts[d]);
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
void compute_degree_per_vertex(
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
	tri_vert.push_back(A);
	tri_vert.push_back(B);
	tri_vert.push_back(C);
}


// Triangulate quads based on their angles
// Param 1: set of NON_DEGENERATE quads
void QTRIANGULATE::triangulate_quad_angle_based(
	std::vector<VERTEX_INDEX> & non_degen_quad_vert, // only non degenerate quads
	std::vector<VERTEX_INDEX> & tri_vert,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
	const std::vector<COORD_TYPE> & vertex_coord
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

	for (int q=0; q < num_quad; q++)
	{
		int w1=0;// has the index NOT the vertex
		bool flag_deg2 = false;
		int vertex;
		for (int j=0; j<VERT_PER_QUAD; j++)
		{
			vertex = QCOLLAPSE::find_vertex(collapse_map,non_degen_quad_vert[q*VERT_PER_QUAD+j]);
			if (iso_vlist[vertex].ver_degree == 2)
			{
				w1 = j;
				flag_deg2 = true;
				//break;
			}
		}

		//
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
			num_non_degen=0;

		if ( num_non_degen == 0 )
			continue;
		else if (num_non_degen == 4 )
		{
			// do computation
			flag_non_degen_quad = true;
		}
		else // triangle
		{
			for (int d=0;d<num_non_degen;d++)
				tri_vert.push_back(temp_tri[d]);
		}
		//
		if (flag_non_degen_quad)
		{
			vertex = non_degen_quad_vert[q*VERT_PER_QUAD+ (w1)%VERT_PER_QUAD];
			vertex = QCOLLAPSE::find_vertex(collapse_map,vertex);
			int B = non_degen_quad_vert[q*VERT_PER_QUAD+ (w1+2)%VERT_PER_QUAD];
			B = QCOLLAPSE::find_vertex(collapse_map,B);
			int C = non_degen_quad_vert[q*VERT_PER_QUAD+ (w1+3)%VERT_PER_QUAD];
			C = QCOLLAPSE::find_vertex(collapse_map,C);
			int D = non_degen_quad_vert[q*VERT_PER_QUAD+ (w1+1)%VERT_PER_QUAD];
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

			if (flag_deg2)
			{	
				if (cos_min_1 < cos_min_2)
				{
					push_triangle(vertex, B, C, tri_vert);
					push_triangle(vertex, B, D, tri_vert);
					//tri_vert.push_back(vertex);
					//tri_vert.push_back(B);
					//tri_vert.push_back(C);
					//tri_vert.push_back(vertex);
					//tri_vert.push_back(D);
					//tri_vert.push_back(B);

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
					//tri_vert.push_back(B);
					//tri_vert.push_back(C);
					//tri_vert.push_back(D);
					push_triangle(B, C, D, tri_vert);
				}

			}
			else
			{
				if (cos_min_1 < cos_min_2)
				{
					//tri_vert.push_back(vertex);
					//tri_vert.push_back(B);
					//tri_vert.push_back(C);

					//tri_vert.push_back(vertex);
					//tri_vert.push_back(D);
					//tri_vert.push_back(B);
					push_triangle(vertex, B, C, tri_vert);
					push_triangle(vertex, B, D, tri_vert);

				}
				else
				{
					/*tri_vert.push_back(D);
					tri_vert.push_back(B);
					tri_vert.push_back(C);

					tri_vert.push_back(D);
					tri_vert.push_back(C);
					tri_vert.push_back(vertex);*/
					push_triangle(D, B, C, tri_vert);
					push_triangle(D, C, vertex, tri_vert);
				}
			}
		}//if	
	}
	update_tris(collapse_map, tri_vert);
	//reorder back to original 
	IJK::reorder_quad_vertices(non_degen_quad_vert);

}
