#include "qdual_remove_degen.h"
using namespace std;
using namespace IJK;

void QTRIANGULATE::triangulate_non_degen_quads (
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<VERTEX_INDEX> & tri_vert,
	const std::vector<COORD_TYPE> & vertex_coord)
{
	cout <<"triangulate only the non degen quads"<<endl;

	std::vector<VERTEX_INDEX> non_degen_quad_vert;
	remove_degenerate_quads(quad_vert, non_degen_quad_vert, tri_vert, vertex_coord);
	quad_vert = non_degen_quad_vert;

}

// Triangulate all quads
void QTRIANGULATE::triangulate_quads (
	std::vector<VERTEX_INDEX> & quad_vert,
	std::vector<VERTEX_INDEX> & tri_vert,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
	const std::vector<COORD_TYPE> & vertex_coord)
{
	triangulate_non_degen_quads(quad_vert, tri_vert, vertex_coord);
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
	const int num_quads = quad_vert.size()/4;

	for (int q=0;q<num_quads;q++)
	{
		vector<VERTEX_INDEX> non_degen_verts;
		int num_non_degen =0;
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
	/*
	int v1,v2;
	const int count = 4;
	for (int q=0;q<num_quads;q++)
	{
	bool degenerate = false;
	for ( v1 = count-1, v2 = 0; v2 < count; v1 = v2++ ) {
	int endPt1 = quad_vert[q*4+v1];
	int endPt2 = quad_vert[q*4+v2];
	if (endPt1==endPt2)
	{
	degenerate = true;
	break;
	}
	}
	if (!degenerate)
	{
	for (int i=0; i<4; i++)
	non_degen_quad_vert.push_back(quad_vert[q*4+i]);
	}
	}
	*/
	//reorder the quads back to the original.
	IJK::reorder_quad_vertices(quad_vert);
	IJK::reorder_quad_vertices(non_degen_quad_vert);
}



// Compute degree of each vertex
// Only for non degenerate poly
void compute_degree_per_vertex(
	const int vert_per_poly,
	std::vector<VERTEX_INDEX> & poly_vert, // only non degenerate quads
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist		
	)
{

	cout <<"compute_degree_per_vertex"<<endl;

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

// Triangulate quads based on their angles
void QTRIANGULATE::triangulate_quad_angle_based(
	std::vector<VERTEX_INDEX> & non_degen_quad_vert, // only non degenerate quads
	std::vector<VERTEX_INDEX> & tri_vert,
	std::vector<QDUAL::DUAL_ISOVERT> & iso_vlist, 
	const std::vector<COORD_TYPE> & vertex_coord
	)
{
	cout <<"Triangulate "<< endl;
	const int NUM_VERTEX_PER_POLY = 4;
	const int num_quad = non_degen_quad_vert.size()/NUM_VERTEX_PER_POLY;
	IJK::reorder_quad_vertices(non_degen_quad_vert);


	//set vertex_degrees in isovlist
	for (int j=0;j<iso_vlist.size();j++)
	{
		iso_vlist[j].ver_degree=0;
	}
	compute_degree_per_vertex(4, non_degen_quad_vert, iso_vlist);
	compute_degree_per_vertex(3, tri_vert, iso_vlist);

	for (int q=0; q < num_quad; q++)
	{
		int w1=0;// has the index NOT the vertex
		bool flag_deg2 = false;
		int vertex;
		for (int j=0; j<NUM_VERTEX_PER_POLY; j++)
		{
			vertex = non_degen_quad_vert[q*NUM_VERTEX_PER_POLY+j];
			cout <<"j="<<j<<" "<<vertex<<" ";
			if (iso_vlist[vertex].ver_degree == 2)
			{
				w1 = j;
				flag_deg2 = true;
				cout <<"-deg2 ";
				//break;
			}
		}




		cout <<"\n w1="<<w1<<", deg="<< iso_vlist[ non_degen_quad_vert[q*NUM_VERTEX_PER_POLY+w1]].ver_degree;
		vertex = non_degen_quad_vert[q*NUM_VERTEX_PER_POLY+ (w1)%NUM_VERTEX_PER_POLY];
		int B = non_degen_quad_vert[q*NUM_VERTEX_PER_POLY+ (w1+2)%NUM_VERTEX_PER_POLY];
		int C = non_degen_quad_vert[q*NUM_VERTEX_PER_POLY+ (w1+3)%NUM_VERTEX_PER_POLY];
		int D = non_degen_quad_vert[q*NUM_VERTEX_PER_POLY+ (w1+1)%NUM_VERTEX_PER_POLY];
		cout <<" A="<< vertex << " B="<< B << " C="<< C<<" D="<<D<<endl;


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
				tri_vert.push_back(vertex);
				tri_vert.push_back(B);
				tri_vert.push_back(C);
                cout <<"tris  "<<vertex<<","<<B<<","<<C;
				tri_vert.push_back(vertex);
				tri_vert.push_back(D);
				tri_vert.push_back(B);
                cout <<" "<<vertex<<","<<B<<","<<D<<endl;
			}
			else
			{
                cout <<"tris_  "<<B<<","<<C<<","<<D<<endl;
				if (length_of_side(vertex, C, vertex_coord) <
					length_of_side(vertex, D, vertex_coord))
					///***************** PROBLEM
				{
					tri_vert.push_back(B);
					tri_vert.push_back(C);
					tri_vert.push_back(D);
				}
				else
				{
                    tri_vert.push_back(B);
					tri_vert.push_back(C);
					tri_vert.push_back(D);

				}
			}

		}
		else
		{
			if (cos_min_1 < cos_min_2)
			{
				tri_vert.push_back(vertex);
				tri_vert.push_back(B);
				tri_vert.push_back(C);

				tri_vert.push_back(vertex);
				tri_vert.push_back(D);
				tri_vert.push_back(B);

			}
			else
			{
				tri_vert.push_back(D);
				tri_vert.push_back(B);
				tri_vert.push_back(C);
				
				tri_vert.push_back(D);
                tri_vert.push_back(C);
                tri_vert.push_back(vertex);
			}
		}
	}

	//reorder back to original 
	IJK::reorder_quad_vertices(non_degen_quad_vert);

}