#ifndef _QDUALTABLE_
#define _QDUALTABLE_
#include "qdual_types.h"
#include "ijkdualtable.h"



class QDUAL_TABLE:public IJKDUALTABLE::ISODUAL_CUBE_TABLE{
public:
    typedef unsigned char DIR_BITS;
    typedef unsigned char CORNER_INDEX;
    
   // static const int MAX_NUM_ISOVERT_PER_CUBE = 4;
protected:
    /// connect_dir [] = Array  indicating connection directions.
    ///                 size is (MAX_NUM_ISOVERT_PER_CUBE*NumTableEntries()).
    /// Directions 0,1,2 are along the X,Y,Z axis in negative direction. 
    /// Directions 3,4,5 are along the X,Y,Z axis in positive direction.
    /// it - table entry, k - isosurface vertex
    /// j'th bit of connect_dir [IsoVIndex(it,k)] is 1
    ///     if k'th isosurface vertex w in the table entry "it" incident on 
    ///     an isosurface edge pointing away from w in direction j.
    DIR_BITS * connect_dir;

    /// V1[] = Array storing  corner at V^1_w (if V^1_w is  not empty).
    ///     Size is  (MAX_NUM_ISOVERT_PER_CUBE*NumTableEntries())
    /// V1[IsoVindex(it,k)] = Cube Corner (0,1,...,7) in V^1_w 
    ///     where w is the k'th isosurface vertex in table entry "it".
    /// V1[isoVindex(it,k)] is defined only if V^1_w is not empty,
    ///   i.e degree (w) = 3.
    CORNER_INDEX * array_V1;
    
    /// V2[] = Array storing  edge at V^2_w (if V^2_w is  not empty).
    ///     Size is  (MAX_NUM_ISOVERT_PER_CUBE*NumTableEntries())
    /// V2[IsoVindex(it,k)] = Grid edge in V^2_w 
    ///     where w is the k'th isosurface vertex in table entry "it".
    /// V2[isoVindex(it,k)] is defined only if V^2_w is not empty, 
    ///   i.e degree (w) = 4 and w is not flat.
    QDUAL::GRID_EDGE * array_V2;

    void CreateConnectDir();
    void CreateV1();
    void CreateV2();

    /// Create table.
    /// Calls ISODUAL_CUBE_TABLE::Create(), CreateConnectDir() 
    ///   CreateV1() and CreateV2() (in order).
    void Create(const int dimension);

public:
    //Constructors
	QDUAL_TABLE(){};
    QDUAL_TABLE(const int);
    //Destructors
	~QDUAL_TABLE();
	
	IJKDUALTABLE::TABLE_INDEX  IsoVIndex (int it, int k)  const
	{return (it*NamedConstants::MAX_NUM_ISOVERT_PER_CUBE+k);};
    
	DIR_BITS connectDir (int it, int k) const
	{
		return (connect_dir[IsoVIndex(it,k)]);
	}

  CORNER_INDEX V1(int it, int k) const
	{
		return (array_V1[IsoVIndex(it,k)]);
	}

  QDUAL::GRID_EDGE & V2(int it, int k) const
  {
    return(array_V2[IsoVIndex(it,k)]);
  }
};


#endif
