/// \file qdual.cxx
/// Quality dual contouring isosurface generation
/// Version 0.0.1

/*
  QDUAL: Quality Dual Isosurface Generation
  Copyright (C) 2014 Arindam Bhattacharya, Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <unordered_map>
#include "ijkisopoly.txx"
#include "ijkmesh.txx"
#include "ijkmesh_geom.txx"
#include "ijktime.txx"

#include "ijkdualtable.h"
#include "ijkdualtable_ambig.h"

#include "qdual.h"
#include "qdual_datastruct.h"
#include "qdual_extract.h"
#include "qdual_position.h"

#include "qdual_collapse.h"
#include "qdual_remove_degen.h"
#include "qdual_restrictions.h"
#include "qdual_table.h"

using namespace IJK;
using namespace QDUAL;
using namespace QCOLLAPSE; 
using namespace QTRIANGULATE;
using namespace NamedConstants;
using namespace std;

// Helper function for quality dual contouring
// Help function for quality dual contouring.
void setUpIsoVlist
(const DUALISO_DATA & dualiso_data, QDUAL_TABLE &qdual_table,
 IJK::BOOL_GRID<DUALISO_GRID> &cube_boundary_grid,
 std::vector<QDUAL::DUAL_ISOVERT> &iso_vlist);

// **************************************************
// DUAL CONTOURING (HYPERCUBES)
// **************************************************

/// Dual Contouring Algorithm.
void QDUAL::quality_dual_contouring
(const DUALISO_DATA & dualiso_data, const SCALAR_TYPE isovalue, 
 DUAL_ISOSURFACE & dual_isosurface, DUALISO_INFO & dualiso_info)
{
	const int dimension = dualiso_data.ScalarGrid().Dimension();
	const AXIS_SIZE_TYPE * axis_size = dualiso_data.ScalarGrid().AxisSize();
	const bool flag_select_split = dualiso_data.SelectSplitFlag();
	const bool flag_separate_neg = dualiso_data.SeparateNegFlag();
	const VERTEX_POSITION_METHOD vpos_method = 
		dualiso_data.VertexPositionMethod();

	float merge_time = 0.0;
	PROCEDURE_ERROR error("dual_contouring");

	clock_t t_start = clock();

	if (!dualiso_data.Check(error)) { throw error; };

	dual_isosurface.Clear();
	dualiso_info.time.Clear();

	ISO_MERGE_DATA merge_data(dimension, axis_size);

	// List of isosurface vertices. 
	std::vector<DUAL_ISOVERT> iso_vlist;

	bool flag_always_separate_opposite(true);
	IJKDUALTABLE::ISODUAL_CUBE_TABLE 
		isodual_table(dimension, flag_separate_neg, 
                  flag_always_separate_opposite);

	dual_contouring
		(dualiso_data.ScalarGrid(), isovalue, dualiso_data,
     dual_isosurface.isopoly_vert, dual_isosurface.dual_edge, 
     dual_isosurface.vertex_coord, iso_vlist,
     isodual_table, merge_data, dualiso_info);

  dual_isosurface.orth_dir.clear();
  for (int i = 0; i < dual_isosurface.dual_edge.size(); i++) {
    DIRECTION_TYPE edge_dir = 
      DIRECTION_TYPE(dual_isosurface.dual_edge[i]%DIM3);
    dual_isosurface.orth_dir.push_back(edge_dir);
  }

	std::unordered_map<VERTEX_INDEX,VERTEX_INDEX>  diagonalMap;

	if (!dualiso_data.flag_NO_collapse) {
    clock_t t0,t1,t3,t2,t4,t5;

    // set up restriction conditions
    QDUAL_TABLE qdual_table(DIM3);

    IJK::BOOL_GRID<DUALISO_GRID> vertex_boundary_grid;
    vertex_boundary_grid.SetSize(dualiso_data.ScalarGrid());
    vertex_boundary_grid.SetAll(false);
    compute_boundary_grid(vertex_boundary_grid);


    IJK::BOOL_GRID<DUALISO_GRID> cube_boundary_grid;
    cube_boundary_grid.SetSize(dualiso_data.ScalarGrid());
    cube_boundary_grid.SetAll(false);
    flag_boundary_cubes(cube_boundary_grid);

    ////set variables in iso_vlist
    setUpIsoVlist( dualiso_data, qdual_table, cube_boundary_grid, iso_vlist);

    //store the original quads 
    vector <QUAD_INDEX> origQuadVert = dual_isosurface.isopoly_vert;
    IJK::reorder_quad_vertices(origQuadVert);

    DUALISO_INDEX_GRID first_isov;
    first_isov.SetSize(dualiso_data.ScalarGrid());
    first_isov.SetAll(-1);

    // set up first_isov
    for (int i = iso_vlist.size()-1; i>=0; i--)
      {
        first_isov.Set(iso_vlist[i].cube_index, i);
      }

    t0=clock();
    // Setup sep  vert
    compute_sep_vert(dualiso_data.ScalarGrid(), iso_vlist, qdual_table);
    compute_sep_edges(dualiso_data.ScalarGrid(), iso_vlist, qdual_table);

    vector<VERTEX_INDEX> isolatedList; // keep track of the isolated vertices.
    set_restrictions 
      (dualiso_data,  dualiso_data.ScalarGrid(), isovalue,  
       dual_isosurface.isopoly_vert, iso_vlist, isodual_table, 
       first_isov, qdual_table, dual_isosurface.vertex_coord, 
       dualiso_info, vertex_boundary_grid, isolatedList);
    t1=clock();
    IJK::clock2seconds(t1-t0, dualiso_info.time.set_restrictions);

    if (dualiso_data.flag_collapse_debug)
      {
        using namespace std;
        GRID_COORD_TYPE * c = new GRID_COORD_TYPE [3];
        cout <<"Num isosurface loops "
             << dualiso_info.rs_info.restriction_BList_size << endl;
        for (int l=0; l< dualiso_info.rs_info.restriction_BList_size; l++)
          {
            int v1 = dualiso_info.rs_info.restricted_edges_info[l].first;
            int v2 = dualiso_data.ScalarGrid().NextVertex
              (v1, dualiso_info.rs_info.restricted_edges_info[l].second);

            dualiso_data.ScalarGrid().ComputeCoord(v1,c);
            cout <<c[0]<<" "<<c[1]<<" "<<c[2] <<" - ";
            dualiso_data.ScalarGrid().ComputeCoord(v2,c);
            cout <<c[0]<<" "<<c[1]<<" "<<c[2]<<endl;
          }
        cout << "Num isosurface boxes "
             << dualiso_info.rs_info.restriction_CList_size << endl;
        for (int l=0; l < dualiso_info.rs_info.restriction_CList_size; l++)
          {
            int v = dualiso_info.rs_info.restricted_vertex_info[l];

            dualiso_data.ScalarGrid().ComputeCoord(v,c);
            cout <<c[0]<<" "<<c[1]<<" "<<c[2]<<endl;
          }
        delete[] c;
      }

    // Set up the vertex collapse map.
    const int numVertex = dual_isosurface.vertex_coord.size();
    IJK::ARRAY<VERTEX_INDEX> collapse_map(numVertex,0);
    // Set up collapse_map.
    setup_collapse_map(collapse_map, numVertex);

    t2=clock();

    // Collapse function calls.
    dual_collapse
      (dualiso_data, dualiso_data.ScalarGrid(), dual_isosurface.isopoly_vert, 
       isodual_table, iso_vlist, dual_isosurface.vertex_coord, 
       dual_isosurface.orth_dir, dualiso_data.qdual_epsilon, 
       collapse_map, dualiso_info);
    t3=clock();
    IJK::clock2seconds(t3-t2, dualiso_info.time.dual_collapse);
		
    // Delete isolated vertices.
    if (dualiso_data.flag_delete_isolate) {
      delIsolated
        (dual_isosurface.isopoly_vert, dual_isosurface.dual_edge,
         isolatedList, dualiso_data.ScalarGrid(),
         iso_vlist, first_isov, isodual_table, 
         dualiso_data.flag_collapse_debug);
    }

    // Keep track of the quad indices when collapsing,
    //   used by function triangulate_degenerate_quads
    //   and for keeping track of quads that share diagonal vertices.
    unordered_map<QUAD_INDEX, QUAD_INDEX> track_quad_indices;
    if (dualiso_data.use_quad_tri_mesh) {

      dual_isosurface.flag_has_degen_quads = 
        triangulate_degenerate_quads 
        (dual_isosurface.isopoly_vert, dual_isosurface.dual_edge,  
         dual_isosurface.tri_vert, dual_isosurface.vertex_coord, 
         track_quad_indices);
		
			
      if (dual_isosurface.flag_has_degen_quads)
        IJK::reorder_quad_vertices(dual_isosurface.isopoly_vert);
    }
    else if (dualiso_data.use_triangle_mesh) {

      t4=clock();
      dual_isosurface.flag_has_degen_quads =  
        triangulate_degenerate_quads 
        (dual_isosurface.isopoly_vert, dual_isosurface.dual_edge,  
         dual_isosurface.tri_vert, dual_isosurface.vertex_coord, 
         track_quad_indices);

      triangulate_quad_angle_based
        (dualiso_data.ScalarGrid(), dual_isosurface.isopoly_vert, 
         dual_isosurface.dual_edge, origQuadVert, 
         dual_isosurface.tri_vert, iso_vlist, 
         dual_isosurface.vertex_coord, cube_boundary_grid, qdual_table, 
         dual_isosurface.orth_dir, track_quad_indices, collapse_map,
         dualiso_data.flag_collapse_debug);
			
      t5=clock();
      IJK::clock2seconds(t5-t4, dualiso_info.time.triangulate);
    }
		
  }

	clock_t t_end = clock();
	clock2seconds(t_end-t_start, dualiso_info.time.total);
}


// **************************************************
// DUAL CONTOURING (HYPERCUBES)
// **************************************************

// Extract isosurface using Dual Contouring algorithm
// Returns list of isosurface simplex vertices
//   and list of isosurface vertex coordinates
void QDUAL::dual_contouring
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, 
 const DUALISO_DATA_FLAGS & dualiso_data_flags,
 std::vector<VERTEX_INDEX> & quad_vert,
 std::vector<VERTEX_INDEX> & dual_edge,
 std::vector<COORD_TYPE> & vertex_coord,
 std::vector<DUAL_ISOVERT> &iso_vlist,
 IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 MERGE_DATA & merge_data, 
 DUALISO_INFO & dualiso_info)
{
	const int dimension = scalar_grid.Dimension();
	const VERTEX_POSITION_METHOD vertex_position_method
		= dualiso_data_flags.vertex_position_method;
	const bool flag_select_split = dualiso_data_flags.flag_select_split;
	const bool flag_separate_neg = dualiso_data_flags.flag_separate_neg;
	const COORD_TYPE center_offset = 0.1;
	PROCEDURE_ERROR error("dual_contouring");
	clock_t t0, t1, t2, t3;


	t0 = clock();

	quad_vert.clear();
	vertex_coord.clear();
	dualiso_info.time.Clear();

	// List of cubes containing isosurface quadrilateral vertices.
	// quad_cube[iq*4+j] = Index of cube containing j'th vertex of quad iq.
	std::vector<ISO_VERTEX_INDEX> quad_cube;

	// Position on quadrilateral of quadrilateral vertices.
	// facet_vertex[iq*4+j] = Position on quad iquad of j'th vertex of quad iq.
	std::vector<FACET_VERTEX_INDEX> facet_vertex;

	// For each bipolar edge e,
	//   For each vertex iv of quad iq dual to e,
	//     Add cube containing iv to quad_cube[].
	//     Add position of iv on quad iq to facet_vertex[].

	extract_dual_isoquad_and_dual_edge
		(scalar_grid, isovalue, quad_cube, facet_vertex, dual_edge, dualiso_info);
	t1 = clock();

	// List of cubes containing isosurface vertices.
	std::vector<ISO_VERTEX_INDEX> cube_list;

	// List of cubes containing isosurface quadrilateral vertices.
	// quad_cube2[iq*4+j] = Index in cube_list[] of cube 
	//     containing j'th vertex of quad iq.
	// quad_cube[k] and quad_cube2[k] represent the same cube but
	//     quad_cube[k] is the cube index in scalar_grid
	//     while quad_cube2[k] is the location in cube_list[].
	std::vector<ISO_VERTEX_INDEX> quad_cube2;

	// Merge identical cubes in list quad_cube[].
	// Return list of cubes cube_list[] with no duplicate entries.
	// quad_cube2[k] = Location in cube_list[] of quad_cube[k].
	merge_identical(quad_cube, cube_list, quad_cube2, merge_data);
	t2 = clock();



	VERTEX_INDEX num_split;
	IJKDUALTABLE::ISODUAL_CUBE_TABLE_AMBIG_INFO ambig_info(dimension);
	int num_non_manifold_split;
	int num_1_2_change;

	if (flag_select_split) {
		split_dual_isovert_ambig
			(scalar_grid, isodual_table, ambig_info, isovalue, 
       cube_list, quad_cube2, facet_vertex,
       iso_vlist, quad_vert, num_split, 
       num_non_manifold_split, num_1_2_change);

		dualiso_info.multi_isov.num_non_manifold_split = num_non_manifold_split;
		dualiso_info.multi_isov.num_1_2_change = num_1_2_change;
	}
	else {
		// Split isosurface vertices using isodual_table and ambig_info
		//   to ensure that the isosurface is a manifold.
		// Initially, each cube contains exactly one isosurface vertex.
		// iso_vlist[] = Resulting list of isosurface vertices.
		// quad_vert[iq*4+j] = Index in iso_vlist[] of j'th vertex of quad iq.
		split_dual_isovert_manifold
			(scalar_grid, isodual_table, ambig_info, isovalue, 
       cube_list, quad_cube2, facet_vertex,
       iso_vlist, quad_vert, num_split, num_non_manifold_split);

		dualiso_info.multi_isov.num_non_manifold_split = num_non_manifold_split;
	}

	dualiso_info.multi_isov.num_cubes_multi_isov = num_split;
	dualiso_info.multi_isov.num_cubes_single_isov =
		cube_list.size() - num_split;

	if (vertex_position_method == CUBE_CENTER) {
		position_dual_isovertices_near_cube_center_multi
			(scalar_grid, isodual_table, isovalue, iso_vlist, center_offset, 
       vertex_coord);
	}  

  else if (vertex_position_method == RANDOM_POS) {
		// set up restriction conditions
		QDUAL_TABLE qdual_table(DIM3);

    RANDOM_SEED_TYPE seed = dualiso_data_flags.random_seed;
    int num_intervals = dualiso_data_flags.random_num_intervals;
    bool flag_V1w_close = dualiso_data_flags.flag_V1w_close;
    position_dual_isovertices_random
			(scalar_grid, qdual_table, iso_vlist, seed, flag_V1w_close, 
       num_intervals, vertex_coord);
  }	
  else {

		// Compute the coordinates of each isosurface vertex iv in iso_vlist[]
		//   as the centroid of the intersection points of the quadrilaterals
		//   incident on vertex iv and the cube edges.
		position_dual_isovertices_centroid_multi
			(scalar_grid, isodual_table, isovalue, iso_vlist, vertex_coord);
	}


	t3 = clock();

	// store times
	clock2seconds(t1-t0, dualiso_info.time.extract);
	clock2seconds(t2-t1, dualiso_info.time.merge);
	clock2seconds(t3-t2, dualiso_info.time.position);
	clock2seconds(t3-t0, dualiso_info.time.total);
}


// **************************************************
// CONVERT QUADRILATERALS TO TRIANGLES
// **************************************************

/// Convert quadrilaterals (embedded in 3D) to triangles.
/// Assumes input isosurface is quadrilaterals embedded in 3D.
void QDUAL::convert_quad_to_tri
(const DUALISO_DATA & dualiso_data, DUAL_ISOSURFACE & dual_isosurface)
{
	const int NUM_QUAD_VERTICES = 4;
	const int DIM3 = 3;
	const int numv_per_iso_poly = dual_isosurface.NumVerticesPerIsoPoly();
	const COORD_TYPE max_smal_magnitude = dualiso_data.MaxSmallMagnitude();
	VERTEX_INDEX_ARRAY quad_vert(dual_isosurface.isopoly_vert);
	IJK::PROCEDURE_ERROR error("convert_quad_to_tri");

	if (numv_per_iso_poly != NUM_QUAD_VERTICES) {
		error.AddMessage("Number of vertices per isosurface polygon is ",
                     numv_per_iso_poly, ".");
		error.AddMessage("Number of vertices per isosurface polygon must be ",
                     NUM_QUAD_VERTICES, ".");
		throw error;
	}

	IJK::reorder_quad_vertices(quad_vert);

	if (dualiso_data.QuadTriangulationMethod() == SPLIT_MAX_ANGLE) {
		IJK::triangulate_quad_split_max_angle
			(DIM3, dual_isosurface.vertex_coord, quad_vert,
       dualiso_data.MaxSmallMagnitude(), dual_isosurface.tri_vert);
	}
	else {
		IJK::triangulate_quad(quad_vert, dual_isosurface.tri_vert);
	}
}

// Help function for quality dual contouring.
void setUpIsoVlist(
                   const DUALISO_DATA & dualiso_data,
                   QDUAL_TABLE &qdual_table,
                   IJK::BOOL_GRID<DUALISO_GRID> &cube_boundary_grid,
                   //returns
                   std::vector<QDUAL::DUAL_ISOVERT> &iso_vlist
                   )
{
	const int iso_vlistSize = iso_vlist.size();
	//set variables in iso_vlist
	IJK::ARRAY<COORD_TYPE> a (DIM3,0);
	bool flag_boundary = false;

	for (int j=0;j<iso_vlist.size();j++) {
    iso_vlist[j].ver_degree=0;
    iso_vlist[j].flag_restrictionC = false;
    iso_vlist[j].flag_isolated = false;
    iso_vlist[j].flag_fixed = false;

    (dualiso_data.ScalarGrid()).ComputeCoord
      (iso_vlist[j].cube_index, &(a[0]));

    //Set iso_vlist[j].cube_coord
    for (int d=0; d<DIM3; d++)
      iso_vlist[j].cube_coord.push_back(a[d]);
    flag_boundary = false;
    //this ones needs the cube boundaries
    isBoundaryIsoVertex(j, iso_vlist, cube_boundary_grid,
                        qdual_table, flag_boundary);

    if (flag_boundary)
      iso_vlist[j].restricted_facets=255;
    else
      iso_vlist[j].restricted_facets=0;
  }
}
