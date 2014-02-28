/// \file qdual.cxx
/// Quality dual contouring isosurface generation
/// Version 0.1.0

/*
QDUAL: Quality Dual Contouring Isosurface Generation
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


#include <iostream>

#include "qdualIO.h"
#include "qdual.h"

using namespace IJK;
using namespace QDUAL;

using namespace std;

// local subroutines
void memory_exhaustion();
void construct_isosurface
	(const IO_INFO & io_info, const DUALISO_DATA & dualiso_data,
	DUALISO_TIME & dualiso_time, IO_TIME & io_time);


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
time_t start_time;
	time(&start_time);

	DUALISO_TIME dualiso_time;
	IO_TIME io_time = {0.0, 0.0};
	IO_INFO io_info;
	IJK::ERROR error;
	try {

		std::set_new_handler(memory_exhaustion);

		parse_command_line(argc, argv, io_info);

		DUALISO_SCALAR_GRID full_scalar_grid;
		NRRD_HEADER nrrd_header;
		read_nrrd_file
			(io_info.input_filename, full_scalar_grid,  nrrd_header, io_time);

		if (!check_input(io_info, full_scalar_grid, error)) 
		{ throw(error); };

		nrrd_header.GetSpacing(io_info.grid_spacing);

		// set DUAL datastructures and flags
		DUALISO_DATA dualiso_data;

		// Note: dualiso_data.SetScalarGrid must be called before set_mesh_data.
		dualiso_data.SetScalarGrid
			(full_scalar_grid, io_info.flag_subsample, io_info.subsample_resolution,
			io_info.flag_supersample, io_info.supersample_resolution);
		set_dualiso_data(io_info, dualiso_data, dualiso_time);
		report_num_cubes(full_scalar_grid, io_info, dualiso_data);

		construct_isosurface(io_info, dualiso_data, dualiso_time, io_time);

		if (io_info.report_time_flag) {

			time_t end_time;
			time(&end_time);
			double total_elapsed_time = difftime(end_time, start_time);

			cout << endl;
			report_time(io_info, io_time, dualiso_time, total_elapsed_time);
		};

	} 
	catch (ERROR & error) {
		if (error.NumMessages() == 0) {
			cerr << "Unknown error." << endl;
		}
		else { error.Print(cerr); }
		cerr << "Exiting." << endl;
		exit(20);
	}
	catch (...) {
		cerr << "Unknown error." << endl;
		exit(50);
	};

}

void construct_isosurface
	(const IO_INFO & io_info, const DUALISO_DATA & dualiso_data,
	DUALISO_TIME & dualiso_time, IO_TIME & io_time)
{
	const int dimension = dualiso_data.ScalarGrid().Dimension();
	const int num_cube_vertices = dualiso_data.ScalarGrid().NumCubeVertices();
	const int num_facet_vertices = dualiso_data.ScalarGrid().NumFacetVertices();
	const int num_cubes = dualiso_data.ScalarGrid().ComputeNumCubes();

	io_time.write_time = 0;
	for (unsigned int i = 0; i < io_info.isovalue.size(); i++) {

		const SCALAR_TYPE isovalue = io_info.isovalue[i];

		DUAL_ISOSURFACE dual_isosurface(num_facet_vertices);
		DUALISO_INFO dualiso_info(dimension);
		dualiso_info.grid.num_cubes = num_cubes;

		quality_dual_contouring(dualiso_data, isovalue, dual_isosurface, dualiso_info);
		dualiso_time.Add(dualiso_info.time);

		OUTPUT_INFO output_info;
		set_output_info(io_info, i, output_info);

		int grow_factor = 1;
		int shrink_factor = 1;
		if (io_info.flag_subsample) 
		{ grow_factor = io_info.subsample_resolution; }
		if (io_info.flag_supersample) 
		{ shrink_factor = io_info.supersample_resolution; }

		rescale_vertex_coord(grow_factor, shrink_factor, io_info.grid_spacing,
			dual_isosurface.vertex_coord);
		
		if (dimension == 3 && dualiso_data.UseTriangleMesh() && dualiso_data.flag_NO_collapse) {
			convert_quad_to_tri(dualiso_data, dual_isosurface);
		}
		
		output_dual_isosurface
			(output_info, dualiso_data, dual_isosurface, dualiso_info, io_time);
	}
}

void memory_exhaustion()
{
	cerr << "Error: Out of memory.  Terminating program." << endl;
	exit(10);
}

