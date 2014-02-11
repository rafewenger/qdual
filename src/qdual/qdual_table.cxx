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
	cout <<"In function CreateConnrctDir\n";
	cout <<"table entries " << num_table_entries <<endl;
	int d1coef[NUM_CUBE_FACET_VERT] = {0,1,0,1}; 
	int d2coef [NUM_CUBE_FACET_VERT] = {0,0,1,1};
	const int connect_dir_size = num_table_entries * MAX_NUM_ISOVERT_PER_CUBE;
	cout <<"conenct dir size "<<num_table_entries * MAX_NUM_ISOVERT_PER_CUBE<< endl;
	connect_dir = new DIR_BITS [num_table_entries * MAX_NUM_ISOVERT_PER_CUBE]();
	for (int k=0; k < connect_dir_size-1; k++)
	{
		connect_dir[k]=0;
	}
	for (int it=0; it <= num_table_entries-1; it++)
	{
		cout <<" table entry "<<it <<",\n ";
		for (int ie=0; ie <= NUM_CUBE_EDGES-1; ie++)
		{
			if (IsBipolar(it,ie))
			{
				int edgeDir = ie/NUM_CUBE_FACET_VERT;
				cout <<"edg="<< ie<<"\n ";
				cout <<" dir="<< edgeDir <<",";
				int d1 = (edgeDir + 1 ) % DIM3;
				int d2 = (edgeDir + 2 ) % DIM3;
				int j=ie%NUM_CUBE_FACET_VERT;
				cout <<" j="<<j<<", ";
				int signedDir1 = (d1 + d1coef[j]*DIM3)%(2*DIM3);
				int signedDir2 = (d2 + d2coef[j]*DIM3)%(2*DIM3);
				cout <<"d1="<<d1<<",d2="<<d2<<", ";
				cout <<"signed_ d1="<<signedDir1<<",d2="<<signedDir2;

				int ipatch =IncidentIsoVertex(it,ie);
				int k = IsoVIndex(it,ipatch);
				cout <<" ipatch= "<<ipatch <<", k="<<k;

				cout <<"\n";
				cout << "connect_dir " << int (connect_dir[k]); 
				connect_dir[k]=connect_dir[k]|(1<<signedDir1)|(1<<signedDir2);
				cout <<" "<< int (connect_dir[k]) <<endl;
			} // is bipolar 
		}
		cout <<" *****************\n";
	}
}

void QDUAL_TABLE::Create(const int dimension)
{
	ISODUAL_CUBE_TABLE::Create(dimension);
	CreateConnectDir();
	cout <<"createV1"<<endl;
	CreateV1();
}

void QDUAL_TABLE::CreateV1() 
{
	std::vector<unsigned char> T(64, 255);
	const int V1_size =  num_table_entries * MAX_NUM_ISOVERT_PER_CUBE;
	for (int i_corner=0; i_corner < NUM_CUBE_VERT-1; i_corner++)
	{
		TABLE_INDEX it = (1<<i_corner);
		cout <<"table index "<< it << endl;
		int k = IsoVIndex(it, 0);
		cout <<"k "<<  k << endl;
		int edgeFlag = connect_dir[k];
		cout <<"edge_flag "<< edgeFlag<<endl;
		T[edgeFlag]= i_corner; 
		cout <<"T["<<edgeFlag <<"]="<<int(T[edgeFlag])<<endl;
	}

	array_V1 = new CORNER_INDEX [V1_size]();

	for (int k = 0; k < V1_size; k++)
	{
		int edgeFlag = connect_dir[k];
		if(T[edgeFlag]!=255)
			// Quads incident on the isosurface vertex separate one cube corner 
				// from other cube corners
		{
			int icorner = T[edgeFlag];
			array_V1[k]=icorner;
		}
	}
}