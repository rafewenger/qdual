
#include "qdual_table.h"

using namespace std;
using namespace NamedConstants;



// Constructor
QDUAL_TABLE::QDUAL_TABLE(const int dimension)
{
	Create(dimension);
}

void QDUAL_TABLE::CreateConnectDir()
{
	int d1coef[NUM_CUBE_FACET_VERT] = {0,1,0,1}; 
	int d2coef [NUM_CUBE_FACET_VERT] = {0,0,1,1};
	const int connect_dir_size = num_table_entries * MAX_NUM_ISOVERT_PER_CUBE;
	
	connect_dir = new DIR_BITS [num_table_entries * MAX_NUM_ISOVERT_PER_CUBE]();
	for (int k=0; k < connect_dir_size-1; k++)
    { connect_dir[k]=0;	}

	for (int it=0; it <= num_table_entries-1; it++)
	{
		for (int ie=0; ie <= NUM_CUBE_EDGES-1; ie++)
		{
			if (IsBipolar(it,ie))
			{
				int edgeDir = ie/NUM_CUBE_FACET_VERT;
				int d1 = (edgeDir + 1 ) % DIM3;
				int d2 = (edgeDir + 2 ) % DIM3;
				int j=ie%NUM_CUBE_FACET_VERT;

				int signedDir1 = (d1 + d1coef[j]*DIM3)%(2*DIM3);
				int signedDir2 = (d2 + d2coef[j]*DIM3)%(2*DIM3);


				int ipatch =IncidentIsoVertex(it,ie);
				int k = IsoVIndex(it,ipatch);

				connect_dir[k]=connect_dir[k]|(1<<signedDir1)|(1<<signedDir2);

			} // is bipolar 
		}
	}
}

void QDUAL_TABLE::Create(const int dimension)
{
	ISODUAL_CUBE_TABLE::Create(dimension);
	CreateConnectDir();
	
	CreateV1();
	CreateV2();
}

void QDUAL_TABLE::CreateV1() 
{
	std::vector<unsigned char> T(64, 255);
	const int V1_size =  num_table_entries * MAX_NUM_ISOVERT_PER_CUBE;

	array_V1 = new CORNER_INDEX [V1_size]();
	for (int k=0; k < V1_size; k++)
	{ array_V1[k]=255;}

	for (int i_corner=0; i_corner < NUM_CUBE_VERT; i_corner++)
	{
		TABLE_INDEX it = (1<<i_corner);
		int k = IsoVIndex(it, 0);
		int edgeFlag = connect_dir[k];
		T[edgeFlag] = i_corner; 
	}

	for (int it = 0; it < this->NumTableEntries(); it++) {
    for (int j = 0; j < this->NumIsoVertices(it); j++) {
      int k = it*MAX_NUM_ISOVERT_PER_CUBE+j;
      int edgeFlag = connect_dir[k];
      if (T[edgeFlag] != 255) {
        // Quads incident on the isosurface vertex separate one cube corner 
        // from other cube corners
        unsigned char icorner = T[edgeFlag];
        array_V1[k]=icorner;
      }
    }
  }
}

void QDUAL_TABLE::CreateV2() 
{
	std::vector<QDUAL::GRID_EDGE> T(64);
	const int V2_size =  num_table_entries * MAX_NUM_ISOVERT_PER_CUBE;
  IJK::CUBE_FACE_INFO<int, int, QDUAL::VERTEX_INDEX> cube(DIM3);

	array_V2 = new QDUAL::GRID_EDGE [V2_size]();
	for (int k=0; k < V2_size; k++)
	{
		array_V2[k].endpoint0 = 255;
    array_V2[k].direction = 0;
	}

  for (int i = 0; i < T.size(); i++) {
    T[i].endpoint0 = 255;
    T[i].direction = 0;
  }

  for (int edge_dir = 0; edge_dir < DIM3; edge_dir++) {
    for (int j = 0; j < cube.NumFacetVertices(); j++) {

      QDUAL::VERTEX_INDEX endpoint0 = cube.FacetVertex(edge_dir, j);
      QDUAL::VERTEX_INDEX endpoint1 = cube.VertexNeighbor(endpoint0, edge_dir);

      TABLE_INDEX it0 = (1<<endpoint0);
      TABLE_INDEX it1 = (1<<endpoint1);
      TABLE_INDEX it = (it0 | it1);
      int k = IsoVIndex(it, 0);
      int edgeFlag = connect_dir[k];
      T[edgeFlag].endpoint0 = endpoint0;
      T[edgeFlag].direction = edge_dir;
    }
  }

	for (int it = 0; it < this->NumTableEntries(); it++) {
    for (int j = 0; j < this->NumIsoVertices(it); j++) {
      int k = it*MAX_NUM_ISOVERT_PER_CUBE+j;
      int edgeFlag = connect_dir[k];
      if (T[edgeFlag].endpoint0 != 255) {
        // Quads incident on the isosurface vertex separate one cube edge
        // from other cube edges
        array_V2[k] = T[edgeFlag];
      }
    }
  }
}

/// Destructor
QDUAL_TABLE::~QDUAL_TABLE()
{
  delete[] connect_dir;
  delete[] array_V1;
  delete[] array_V2;

  connect_dir = NULL;
  array_V1 = NULL;
  array_V2 = NULL;
}

