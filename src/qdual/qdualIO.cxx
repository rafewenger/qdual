// \file ijkdualIO.cxx
/// IO classes and routines for qdual.

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

#include <assert.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "ijkcommand_line.txx"
#include "ijkIO.txx"
#include "ijkstring.txx"
#include "ijkmesh.txx"

#include "qdualIO.h"

using namespace IJK;
using namespace QDUAL;

using namespace std;

// **************************************************
// PARSE COMMAND LINE
// **************************************************

// local namespace
namespace {

	typedef enum
	{SUBSAMPLE_PARAM, SUPERSAMPLE_PARAM, POSITION_PARAM, 
	SINGLE_ISOV_PARAM,  MULTI_ISOV_PARAM, 
	SPLIT_NON_MANIFOLD_PARAM, SELECT_SPLIT_PARAM,
	SEP_NEG_PARAM, SEP_POS_PARAM, EPSILON, MOVE_VERTEX, DEL_ISOLATE,
	TRIMESH_PARAM, UNIFORM_TRIMESH_PARAM,COLLAPSE_DEBUG,COLLAPSE_INFO,
	NO_COLLAPSE, NO_CAPCOL, QTMESH,
	NO_RESTR_AB,NO_RESTR_B, NO_RESTR_C,
  SEED_PARAM,
	HELP_PARAM, OFF_PARAM, IV_PARAM, 
	OUTPUT_FILENAME_PARAM, STDOUT_PARAM, 
	NOWRITE_PARAM, SILENT_PARAM, TIME_PARAM, UNKNOWN_PARAM} PARAMETER;
	const char * parameter_string[] = 
	{"-subsample", "-supersample", "-position", 
	"-single_isov", "-multi_isov", "-split_non_manifold", "-select_split",
	"-sep_neg", "-sep_pos","-epsilon","-move_vertex","-del_isolate",
	"-trimesh", "-uniform_trimesh","-collapse_debug", "-collapse_info",
	"-no_collapse","-no_capcol","-qt_mesh",
	"-no_res_AB","-no_res_B","-no_res_C",
  "-seed",
	"-help", "-off", "-iv", "-o", "-stdout",
	"-nowrite", "-s", "-time", "-unknown"};

	PARAMETER get_parameter_token(char * s)
		// convert string s into parameter token
	{
		for (int i = 0; i < int(UNKNOWN_PARAM); i++)
			if (strcmp(parameter_string[i], s) == 0)
				return(PARAMETER(i));
		return(UNKNOWN_PARAM);
	}

	INTERPOLATION_TYPE get_interpolation_type(char * s)
		// convert string s into parameter token
	{
		INTERPOLATION_TYPE type = LINEAR_INTERPOLATION;

		if (strcmp(s, "linear") == 0) 
		{ type = LINEAR_INTERPOLATION; }
		else if (strcmp(s, "multilinear") == 0) 
		{ type = MULTILINEAR_INTERPOLATION; }
		else {
			cerr << "Error in input parameter -interpolate.  Illegal interpolation type: " 
				<< s << "." << endl;
			exit(1030);
		}

		return(type);
	}

	VERTEX_POSITION_METHOD get_vertex_position_method(char * s)
		// convert string s into parameter token
	{
		VERTEX_POSITION_METHOD method = CENTROID_EDGE_ISO;
		string str = s;

		if (str == "cube_center") {
			method = CUBE_CENTER;
		}
		else if (str == "centroid") {
			method = CENTROID_EDGE_ISO;
		}
    else if (str == "random") {
			method = RANDOM_POS;
    }
		else {
			cerr << "Error in input parameter -position.  Illegal position method: " 
				<< str << "." << endl;
			exit(1030);
		}

		return(method);
	}

}

void QDUAL::parse_command_line(int argc, char **argv, IO_INFO & io_info)
	// parse command line
	// control parameters, followed by one or more isovalues, 
	// followed by input file name
{
	IJK::ERROR error;

	if (argc == 1) { usage_error(); };

	int iarg = 1;
	while (iarg < argc && argv[iarg][0] == '-') {
		PARAMETER param = get_parameter_token(argv[iarg]);
		if (param == UNKNOWN_PARAM) break;

		switch(param) {

		case SUBSAMPLE_PARAM:
			io_info.subsample_resolution = get_arg_int(iarg, argc, argv, error);
			io_info.flag_subsample = true;
			iarg++;
			break;

		case SUPERSAMPLE_PARAM:
			io_info.supersample_resolution = get_arg_int(iarg, argc, argv, error);
			io_info.flag_supersample = true;
			iarg++;
			break;

		case POSITION_PARAM:
			iarg++;
			if (iarg >= argc) usage_error();
			io_info.vertex_position_method = 
				get_vertex_position_method(argv[iarg]);
			break;

		case SELECT_SPLIT_PARAM:
			io_info.flag_select_split = true;
			io_info.allow_multiple_iso_vertices = true;
			break;

		case SEP_NEG_PARAM:
			io_info.allow_multiple_iso_vertices = true;
			io_info.flag_separate_neg = true;
			break;

		case SEP_POS_PARAM:
			io_info.allow_multiple_iso_vertices = true;
			io_info.flag_separate_neg = false;
			break;

		case TRIMESH_PARAM:
			io_info.use_triangle_mesh = true;
			io_info.quad_tri_method = SPLIT_MAX_ANGLE;
			break;

		case QTMESH:
			io_info.use_quad_tri_mesh = true;
			break;

		case EPSILON:
			iarg++;
			if (iarg >= argc) usage_error();
			io_info.qdual_epsilon = atof(argv[iarg]);
			break;

		case MOVE_VERTEX:
			io_info.flag_move_vertices = true;
			break;

		case DEL_ISOLATE:
			io_info.flag_delete_isolate = true;
			break;


		case NO_CAPCOL:
			io_info.flag_cap_col = false;
			break;

    case SEED_PARAM:
      iarg++;
			if (iarg >= argc) usage_error();
			io_info.random_seed = atoi(argv[iarg]);
			break;

		case UNIFORM_TRIMESH_PARAM:
			io_info.use_triangle_mesh = true;
			io_info.quad_tri_method = UNIFORM_TRI;
			break;

		case OFF_PARAM:
			io_info.output_format = OFF;
			break;

		case IV_PARAM:
			io_info.output_format = IV;
			break;

		case OUTPUT_FILENAME_PARAM:
			iarg++;
			if (iarg >= argc) usage_error();
			io_info.output_filename = argv[iarg];
			break;

		case STDOUT_PARAM:
			io_info.use_stdout = true;
			break;

		case NOWRITE_PARAM:
			io_info.nowrite_flag = true;
			break;

		case SILENT_PARAM:
			io_info.flag_silent = true;
			break;

		case TIME_PARAM:
			io_info.report_time_flag = true;
			break;

		case HELP_PARAM:
			help();
			break;

		case COLLAPSE_DEBUG:
			io_info.flag_collapse_debug = true;
			break;

		case COLLAPSE_INFO:
			io_info.flag_collapse_info = true;
			break;

		case NO_COLLAPSE:
			io_info.flag_NO_collapse = true;
			break;

		case NO_RESTR_AB:
			io_info.flag_no_restriction_AB = true;
			break;
		case NO_RESTR_B:
			io_info.flag_no_restriciton_B = true;
			break;
		case NO_RESTR_C:
			io_info.flag_no_restriciton_C = true;
			break;

		};

		iarg++;
	};

	// remaining parameters should be list of isovalues followed
	// by input file name

	// check for more parameter tokens
	for (int j = iarg; j+1 < argc; j++) {
		if ((get_parameter_token(argv[j]) != UNKNOWN_PARAM) ||
			(argv[j][0] == '-' && !is_type<float>(argv[j]))) {
				// argv[iarg] is not an isovalue
				cerr << "Error. Illegal parameter: " << argv[iarg] << endl;
				cerr << endl;
				usage_error();
		}
	}

	if (iarg+2 > argc) {
		cerr << "Error.  Missing input isovalue or input file name." << endl;
		cerr << endl;
		usage_error();
	};

	// store isovalues
	for (int j = iarg; j+1 < argc; j++) {
		SCALAR_TYPE value;
		if (!IJK::string2val(argv[j], value)) {
			cerr << "Error. \"" << argv[j] << "\" is not a valid input isovalue." 
				<< endl;
			usage_error();
		}

		io_info.isovalue_string.push_back(argv[j]);
		io_info.isovalue.push_back(value);
	}

	io_info.input_filename = argv[argc-1];

	if (io_info.flag_subsample && io_info.subsample_resolution <= 1) {
		cerr << "Error.  Subsample resolution must be an integer greater than 1."
			<< endl;
		exit(230);
	};

	if (io_info.output_filename != NULL && io_info.use_stdout) {
		cerr << "Error.  Can't use both -o and -stdout parameters."
			<< endl;
		exit(230);
	};

	if (io_info.flag_subsample && io_info.flag_supersample) {
		cerr << "Error.  Can't use both -subsample and -supersample parameters."
			<< endl;
		exit(555);
	}
}

// Check input information/flags.
bool QDUAL::check_input
	(const IO_INFO & io_info, 
	const DUALISO_SCALAR_GRID_BASE & scalar_grid,
	IJK::ERROR & error)
{
	// Construct isosurface
	if (io_info.isovalue.size() > 1 && io_info.use_stdout) {
		error.AddMessage
			("Error.  Cannot use stdout for more than one isovalue.");
		return(false);
	}

	if (io_info.isovalue.size() > 1 && io_info.output_filename != NULL) {
		error.AddMessage
			("Error.  Cannot specify output file for more than one isovalue.");
		return(false);
	}

	return(true);
}

// **************************************************
// READ NEARLY RAW RASTER DATA (nrrd) FILE
// **************************************************

void QDUAL::read_nrrd_file
	(const char * input_filename, DUALISO_SCALAR_GRID & scalar_grid, 
	NRRD_HEADER & nrrd_header, IO_TIME & io_time)
{
	ELAPSED_TIME wall_time;
	GRID_NRRD_IN<int, AXIS_SIZE_TYPE> nrrd_in;
	IJK::PROCEDURE_ERROR error("read_nrrd_file");

	nrrd_in.ReadScalarGrid(input_filename, scalar_grid, nrrd_header, error);
	if (nrrd_in.ReadFailed()) { throw error; }

	io_time.read_nrrd_time = wall_time.getElapsed();
}

// **************************************************
// PATH_DELIMITER
// **************************************************

namespace {

#ifdef _WIN32
	const char PATH_DELIMITER = '\\';
#else
	const char PATH_DELIMITER = '/';
#endif
}


// **************************************************
// OUTPUT ISOSURFACE
// **************************************************

void QDUAL::output_dual_isosurface
	(const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
	const DUAL_ISOSURFACE & dual_isosurface,
	const DUALISO_INFO & dualiso_info, IO_TIME & io_time)
{
	if (output_info.use_triangle_mesh) {
		output_dual_tri_isosurface
			(output_info, dualiso_data, dual_isosurface.vertex_coord, 
			dual_isosurface.tri_vert, dualiso_info, io_time);
	}
	else if (output_info.use_quad_tri_mesh && dual_isosurface.flag_has_degen_quads)
	{
		output_dual_tri_quad_isosurface
			(output_info, dualiso_data, dual_isosurface.vertex_coord, 
			dual_isosurface.tri_vert, dual_isosurface.isopoly_vert,
			dualiso_info, io_time);
	}
	else {
		output_dual_isosurface
			(output_info, dualiso_data, dual_isosurface.vertex_coord, 
			dual_isosurface.isopoly_vert, dualiso_info, io_time);
	}
}

void QDUAL::output_dual_isosurface
	(const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
	const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
	const DUALISO_INFO & dualiso_info, IO_TIME & io_time)
{
	if (!output_info.use_stdout && !output_info.flag_silent) {
		report_iso_info(output_info, dualiso_data, 
			vertex_coord, slist, dualiso_info);
	}

	if (!output_info.nowrite_flag) 
	{ write_dual_mesh(output_info, vertex_coord, slist, io_time); }
}

/// Output isosurface of triangles.
void QDUAL::output_dual_tri_isosurface
	(const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
	const std::vector<COORD_TYPE> & vertex_coord,
	const std::vector<VERTEX_INDEX> & tri_vert,
	const DUALISO_INFO & dualiso_info, IO_TIME & io_time)
{
	if (!output_info.use_stdout && !output_info.flag_silent) {
		report_iso_info(output_info, dualiso_data,
			vertex_coord, tri_vert, dualiso_info);
	}

	if (!output_info.nowrite_flag) {
		write_dual_tri_mesh(output_info, vertex_coord, tri_vert, io_time);
	}
}

/// Output isosurface of quadrilaterals and triangles.
void QDUAL::output_dual_tri_quad_isosurface
	(const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
	const std::vector<COORD_TYPE> & vertex_coord,
	const std::vector<VERTEX_INDEX> & tri_vert,
	const std::vector<VERTEX_INDEX> & quad_vert,
	const DUALISO_INFO & dualiso_info, IO_TIME & io_time)
{
	if (!output_info.use_stdout && !output_info.flag_silent) {
		report_iso_info(output_info, dualiso_data,
			vertex_coord, tri_vert, quad_vert, dualiso_info);
	}

	if (!output_info.nowrite_flag) {
		write_dual_tri_quad_mesh
			(output_info, vertex_coord, tri_vert, quad_vert, io_time);
	}
}

void QDUAL::output_dual_isosurface_color
	(const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
	const DUAL_ISOSURFACE & dual_isosurface,
	const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
	const DUALISO_INFO & dualiso_info, IO_TIME & io_time)
{
	output_dual_isosurface_color
		(output_info, dualiso_data, 
		dual_isosurface.vertex_coord, dual_isosurface.isopoly_vert,
		front_color, back_color, dualiso_info, io_time);
}

void QDUAL::output_dual_isosurface_color
	(const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
	const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
	const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
	const DUALISO_INFO & dualiso_info, IO_TIME & io_time)
{
	if (!output_info.use_stdout && !output_info.flag_silent) {
		report_iso_info(output_info, dualiso_data, vertex_coord, slist, dualiso_info);
	}

	if (!output_info.nowrite_flag) {
		write_dual_mesh_color
			(output_info, vertex_coord, slist, front_color, back_color, io_time);
	}
}

void QDUAL::output_dual_isosurface_color_alternating
	(const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
	const DUAL_ISOSURFACE & dual_isosurface,
	const DUALISO_INFO & dualiso_info, IO_TIME & io_time)
{
	const VERTEX_INDEX num_poly = dual_isosurface.NumIsoPoly();

	IJK::ARRAY<COLOR_TYPE> front_color(4*num_poly);
	IJK::ARRAY<COLOR_TYPE> back_color(4*num_poly);
	set_color_alternating
		(dualiso_data.ScalarGrid(), dual_isosurface.cube_containing_isopoly,
		front_color.Ptr());
	set_color_alternating
		(dualiso_data.ScalarGrid(), dual_isosurface.cube_containing_isopoly,
		back_color.Ptr());

	output_dual_isosurface_color
		(output_info, dualiso_data, dual_isosurface.vertex_coord,
		dual_isosurface.isopoly_vert, 
		front_color.PtrConst(), back_color.PtrConst(),
		dualiso_info, io_time);
}

// **************************************************
// WRITE_DUAL_MESH
// **************************************************

void QDUAL::write_dual_mesh
	(const OUTPUT_INFO & output_info,
	const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & plist)
{
	const int dimension = output_info.dimension;
	const int numv_per_simplex = output_info.num_vertices_per_isopoly;
	const bool use_stdout = output_info.use_stdout;
	ofstream output_file;
	ERROR error_mcube("write_dual_mesh");

	// Output vertices in counter-clockwise order around quadrilateral.
	const bool flag_reorder_quad_vertices = true;

	string ofilename = output_info.output_filename;

	switch (output_info.output_format) {

	case OFF:
		if (!use_stdout) {
			output_file.open(ofilename.c_str(), ios::out);
			if (dimension == 3) {
				ijkoutQuadOFF(output_file, dimension, vertex_coord, plist, 
					flag_reorder_quad_vertices);
			}
			else {
				ijkoutOFF(output_file, dimension, numv_per_simplex,
					vertex_coord, plist);
			}
			output_file.close();
		}
		else {
			if (dimension == 3) {
				ijkoutQuadOFF(dimension, vertex_coord, plist,
					flag_reorder_quad_vertices);
			}
			else {
				ijkoutOFF(dimension, numv_per_simplex, vertex_coord, plist);
			}
		};
		break;

	case IV:
		if (dimension == 3) {
			if (!use_stdout) {
				output_file.open(ofilename.c_str(), ios::out);
				ijkoutIV(output_file, dimension, vertex_coord, plist);
				output_file.close();
			}
			else {
				ijkoutIV(dimension, vertex_coord, plist);
			}
		}
		else throw error_mcube("Illegal dimension. OpenInventor format is only for dimension 3.");
		break;

	default:
		throw error_mcube("Illegal output format.");
		break;
	}

	if (!use_stdout && !output_info.flag_silent)
		cout << "Wrote output to file: " << ofilename << endl;
}

void QDUAL::write_dual_mesh
	(const OUTPUT_INFO & output_info,
	const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & plist,
	IO_TIME & io_time)
{
	ELAPSED_TIME wall_time;

	write_dual_mesh(output_info, vertex_coord, plist);

	io_time.write_time += wall_time.getElapsed();
}

void QDUAL::write_dual_mesh_color
	(const OUTPUT_INFO & output_info,
	const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & plist,
	const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
{
	const int dimension = output_info.dimension;
	const int numv_per_simplex = output_info.num_vertices_per_isopoly;
	const bool use_stdout = output_info.use_stdout;

	ofstream output_file;
	ERROR error_mcube("write_dual_mesh_color");

	string ofilename = output_info.output_filename;

	switch (output_info.output_format) {

	case OFF:
		if (!use_stdout) {
			output_file.open(ofilename.c_str(), ios::out);
			ijkoutColorFacesOFF(output_file, dimension, numv_per_simplex,
				vertex_coord, plist, front_color, back_color);
			output_file.close();
		}
		else {
			ijkoutColorFacesOFF(std::cout, dimension, numv_per_simplex,
				vertex_coord, plist, front_color, back_color);
		};
		break;

	case IV:
		if (dimension == 3) {
			if (!use_stdout) {
				output_file.open(ofilename.c_str(), ios::out);
				ijkoutIV(output_file, dimension, vertex_coord, plist);
				output_file.close();
			}
			else {
				ijkoutOFF(dimension, vertex_coord, plist);
			}
		}
		else throw error_mcube("Illegal dimension. OpenInventor format is only for dimension 3.");
		break;

	default:
		throw error_mcube("Illegal output format.");
		break;
	}

	if (!use_stdout && !output_info.flag_silent)
		cout << "Wrote output to file: " << ofilename << endl;
}

void QDUAL::write_dual_mesh_color
	(const OUTPUT_INFO & output_info,
	const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & plist,
	const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
	IO_TIME & io_time)
{
	ELAPSED_TIME wall_time;

	write_dual_mesh_color(output_info, vertex_coord, plist, 
		front_color, back_color);

	io_time.write_time += wall_time.getElapsed();
}

/// Write dual isosurface triangular mesh.
/// @param output_info Output information.
/// @param vertex_coord List of vertex coordinates.
/// @param tri_vert[] List of triangle vertices.
///        tri_vert[3*i+k] is k'th vertex of triangle i.
void QDUAL::write_dual_tri_mesh
	(const OUTPUT_INFO & output_info,
	const std::vector<COORD_TYPE> & vertex_coord,
	const std::vector<VERTEX_INDEX> & tri_vert)
{
	const int NUMV_PER_QUAD = 4;
	const int NUMV_PER_TRI = 3;
	const int dimension = output_info.dimension;
	const bool use_stdout = output_info.use_stdout;
	ofstream output_file;
	PROCEDURE_ERROR error("write_dual_tri_mesh");

	if (dimension != 3) {
		error.AddMessage("Programming error.  Illegal dimension ", dimension, ".");
		error.AddMessage("   Routine only allowed for dimension 3.");
		throw error;
	}

	string ofilename = output_info.output_filename;

	switch (output_info.output_format) {

	case OFF:
		if (!use_stdout) {
			output_file.open(ofilename.c_str(), ios::out);
			ijkoutOFF(output_file, dimension, NUMV_PER_TRI,
				vertex_coord, tri_vert);
			output_file.close();
		}
		else {
			ijkoutOFF(dimension, NUMV_PER_TRI, vertex_coord, tri_vert);
		};
		break;

	case IV:
		ijkoutIV(dimension, vertex_coord, tri_vert);
		break;

	default:
		throw error("Illegal output format.");
		break;
	}

	if (!use_stdout && !output_info.flag_silent)
		cout << "Wrote output to file: " << ofilename << endl;
}

void QDUAL::write_dual_tri_mesh
	(const OUTPUT_INFO & output_info,
	const vector<COORD_TYPE> & vertex_coord,
	const vector<VERTEX_INDEX> & tri_vert,
	IO_TIME & io_time)
{
	ELAPSED_TIME wall_time;

	write_dual_tri_mesh(output_info, vertex_coord, tri_vert);

	io_time.write_time += wall_time.getElapsed();
}

// **************************************************
// WRITE DUAL QUAD + TRI MESH 
// **************************************************

/* OBSOLETE
void QDUAL::write_dual_quad_tri_mesh
( const OUTPUT_INFO & output_info,
const DUALISO_DATA & dualiso_data,
const DUAL_ISOSURFACE & dual_isosurface
)
{
const int dimension = output_info.dimension;
const int numv_per_simplex = output_info.num_vertices_per_isopoly;
const bool use_stdout = output_info.use_stdout;
const int numv = dual_isosurface.vertex_coord.size()/dimension;
const int num_tri = dual_isosurface.tri_vert.size()/dimension;
const int num_quad = dual_isosurface.isopoly_vert.size()/4;
const int tri = 3;
const int quad = 4;

report_iso_info(output_info, dualiso_data, 
vertex_coord, slist, dualiso_info);

ofstream output_file;
ERROR error_mcube("write_dual_mesh");
string ofilename = output_info.output_filename;
output_file.open(ofilename.c_str(), ios::out);
ijkoutOFF(output_file, dimension, &(dual_isosurface.vertex_coord[0]), numv, 
&(dual_isosurface.tri_vert[0]), tri, num_tri,
&(dual_isosurface.isopoly_vert[0]), quad, num_quad);
output_file.close();
}
*/


void QDUAL::write_dual_tri_quad_mesh
	(const OUTPUT_INFO & output_info,
	const std::vector<COORD_TYPE> & vertex_coord, 
	const std::vector<VERTEX_INDEX> & tri_vert,
	const std::vector<VERTEX_INDEX> & quad_vert,
	IO_TIME & io_time)
{
	const int NUM_VERT_PER_TRI(3);
	const int NUM_VERT_PER_QUAD(4);
	const int dimension = output_info.dimension;
	//IJK::reorder_quad_vertices(quad_vert);
	ofstream output_file;
	ERROR error_mcube("write_dual_mesh");
	string ofilename = output_info.output_filename;
	output_file.open(ofilename.c_str(), ios::out);
	ijkoutOFF(output_file, dimension, vertex_coord,
		tri_vert, NUM_VERT_PER_TRI, quad_vert, NUM_VERT_PER_QUAD);
	output_file.close();
}

// **************************************************
// RESCALE ROUTINES
// **************************************************

namespace {

	void grow_coord(const int scale, vector<COORD_TYPE> & vertex_coord)
	{
		for (unsigned int i = 0; i < vertex_coord.size(); i++) {
			vertex_coord[i] = scale * vertex_coord[i];
		};
	}

	void shrink_coord(const int scale, vector<COORD_TYPE> & vertex_coord)
	{
		for (unsigned int i = 0; i < vertex_coord.size(); i++) {
			vertex_coord[i] = vertex_coord[i]/scale;
		};
	}

	bool unit_spacing(const std::vector<COORD_TYPE> & spacing)
		// return true if spacing not defined or spacing along all axes equals 1.0
	{
		for (unsigned int d = 0; d < spacing.size(); d++) {
			if (!AIR_EXISTS(spacing[d])) { return(true); }
			else if (spacing[d] != 1.0) { return(false); };
		}

		return(true);
	}

	void rescale_coord(const std::vector<COORD_TYPE> & grid_spacing,
		std::vector<COORD_TYPE> & vertex_coord)
	{
		const int dimension = grid_spacing.size();

		if (unit_spacing(grid_spacing)) { return; }

		const VERTEX_INDEX numv = vertex_coord.size()/dimension;
		for (int iv = 0; iv < numv; iv++) {
			for (int d = 0; d < dimension; d++) {
				vertex_coord[iv*dimension+d] *= grid_spacing[d];
			}
		};
	}

}

/// Rescale subsampled/supersampled vertex coordinates.
/// Also rescale to reflect grid spacing.
void QDUAL::rescale_vertex_coord
	(const OUTPUT_INFO & output_info, vector<COORD_TYPE> & vertex_coord)
{
	const int grow_factor = output_info.grow_factor;
	const int shrink_factor = output_info.shrink_factor;
	PROCEDURE_ERROR error("rescale_vertex_coord");

	if (grow_factor <= 0) {
		error.AddMessage("Illegal grow factor ", grow_factor, ".");
		error.AddMessage("  Grow factor must be a positive integer");
	}

	if (shrink_factor <= 0) {
		error.AddMessage("Illegal shrink factor ", shrink_factor, ".");
		error.AddMessage("  Shrink factor must be a positive integer");
	}

	if (output_info.dimension != output_info.grid_spacing.size()) {
		error.AddMessage("Size of grid spacing array does not equal volume dimension.");
		error.AddMessage("  Grid spacing array has ", 
			output_info.grid_spacing.size(), " elements.");
		error.AddMessage("  Volume dimension = ", output_info.dimension, ".");
	}

	if (output_info.grow_factor != 1) 
	{ grow_coord(output_info.grow_factor, vertex_coord); };

	if (output_info.shrink_factor != 1) 
	{ shrink_coord(output_info.shrink_factor, vertex_coord); };

	rescale_coord(output_info.grid_spacing, vertex_coord);
}

/// Rescale subsampled/supersampled vertex coordinates.
/// Also rescale to reflect grid spacing.
void QDUAL::rescale_vertex_coord
	(const int grow_factor, const int shrink_factor,
	const COORD_ARRAY & grid_spacing, COORD_ARRAY & vertex_coord)
{
	PROCEDURE_ERROR error("rescale_vertex_coord");

	if (grow_factor <= 0) {
		error.AddMessage("Illegal grow factor ", grow_factor, ".");
		error.AddMessage("  Grow factor must be a positive integer");
	}

	if (shrink_factor <= 0) {
		error.AddMessage("Illegal shrink factor ", shrink_factor, ".");
		error.AddMessage("  Shrink factor must be a positive integer");
	}

	if (vertex_coord.size() == 0) { return; };

	if (grid_spacing.size() < 1) {
		error.AddMessage("Illegal size ", grid_spacing.size(), 
			" of array grid spacing.");
		error.AddMessage("Size must equal vertex dimension.");
		throw error;
	}

	if (grow_factor != 1) 
	{ grow_coord(grow_factor, vertex_coord); };

	if (shrink_factor != 1) 
	{ shrink_coord(shrink_factor, vertex_coord); };

	rescale_coord(grid_spacing, vertex_coord);
}

// **************************************************
// REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
// **************************************************

void QDUAL::report_num_cubes
	(const DUALISO_GRID & full_scalar_grid, const IO_INFO & io_info, 
	const DUALISO_DATA & dualiso_data)
{
	const int num_grid_cubes = full_scalar_grid.ComputeNumCubes();
	const int num_cubes_in_dualiso_data = 
		dualiso_data.ScalarGrid().ComputeNumCubes();

	if (!io_info.use_stdout && !io_info.flag_silent) {

		if (io_info.flag_subsample) {
			// subsampled grid
			cout << num_grid_cubes << " grid cubes.  "
				<< num_cubes_in_dualiso_data << " subsampled grid cubes." << endl;
		}
		else if (io_info.flag_supersample) {
			// supersample grid
			cout << num_grid_cubes << " grid cubes.  "
				<< num_cubes_in_dualiso_data << " supersampled grid cubes." << endl;
		}
		else {
			// use full_scalar_grid
			cout << num_grid_cubes << " grid cubes." << endl;
		}
	}

}

void QDUAL::report_iso_info
	(const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
	const vector<COORD_TYPE> & vertex_coord, 
	const vector<VERTEX_INDEX> & plist, 
	const DUALISO_INFO & dualiso_info)
{
	const int dimension = output_info.dimension;
	const int numv_per_simplex = output_info.num_vertices_per_isopoly;

	const char * indent4 = "    ";
	string grid_element_name = "cubes";
	if (dimension == 2) { grid_element_name = "squares"; };

	VERTEX_INDEX numv = (vertex_coord.size())/dimension;
	VERTEX_INDEX num_poly = (plist.size())/numv_per_simplex;
	VERTEX_INDEX num_grid_cubes = dualiso_info.grid.num_cubes;
	VERTEX_INDEX num_non_empty_cubes = dualiso_info.scalar.num_non_empty_cubes;

	float percent = 0.0;
	if (num_grid_cubes > 0)
	{ percent = float(num_non_empty_cubes)/float(num_grid_cubes); }
	int ipercent = int(100*percent);
	cout << "  Isovalue " << output_info.isovalue[0] << ".  " 
		<< numv << " isosurface vertices.  "
		<< num_poly << " isosurface polytopes." << endl;

	if (output_info.allow_multiple_iso_vertices) {
		cout << indent4 << "# cubes with single isosurface vertex: "
			<< dualiso_info.multi_isov.num_cubes_single_isov << endl;
		cout << indent4 << "# cubes with multiple isosurface vertices: "
			<< dualiso_info.multi_isov.num_cubes_multi_isov << endl;

		if (output_info.flag_split_non_manifold) {
			cout << indent4 << "# cubes changed to 2 iso vertices to avoid non-manifold edges: "
				<< dualiso_info.multi_isov.num_non_manifold_split << endl;
		}

		if (output_info.flag_select_split) {
			cout << indent4 << "# cubes changed in selecting isosurface patch splits: "
				<< dualiso_info.multi_isov.num_1_2_change << endl;
		}

	}

	if (output_info.flag_collapse_info)
	{
		cout <<"Epsilon is set to: " << output_info.qdual_epsilon << endl;
		cout <<"Restriction Info"<<endl;
		cout <<"	Num isosurface loops: "<< dualiso_info.rs_info.restriction_BList_size << endl;
		cout <<"	Num isosurface boxes: "<< dualiso_info.rs_info.restriction_CList_size << endl;
		cout <<"\nCollapse info" << endl;
		cout <<"	Collapse across facets: "<<dualiso_info.col_info.permitted_facet_restriction<<endl;
		cout <<"	Collapse around edges: "<<dualiso_info.col_info.permitted_edge_restriction<<endl;
		cout <<"	Collapse around vertices: "<<dualiso_info.col_info.permitted_vertex_restriction<<endl;
		cout <<"	Not permitted facet restriction: "<<dualiso_info.col_info.not_permitted_facet_restriction<<endl;
		cout <<"	Not permitted edge restriction: "<<dualiso_info.col_info.not_permitted_edge_restriction<<endl;
		cout <<"	Not permitted vertex restriction: "<< dualiso_info.col_info.not_permitted_vertex_restriction<<endl;

	}
	if (output_info.flag_move_vertices)
	{
		cout <<"Move Vertices Info"<<endl;
		cout <<"	Moved from edges: "<<dualiso_info.mv_info.moveFromEdges<<endl;
		cout <<"	Moved from vertices: "<<dualiso_info.mv_info.moveFromVertices<<endl;
	}
	if(output_info.flag_cap_col)
	{
		cout <<"Cap Collapse Iinfo:"<<endl;
		cout <<"	Number of move to edge "<<dualiso_info.cp_info.moved2Edge <<endl;
		cout <<"	Number of cap quad "<<dualiso_info.cp_info.numCapQuad<<endl;
	}

}

void QDUAL::report_iso_info
	(const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
	const vector<COORD_TYPE> & vertex_coord, 
	const vector<VERTEX_INDEX> & tri_list, 
	const vector<VERTEX_INDEX> & quad_list, 
	const DUALISO_INFO & dualiso_info)
{
	const int dimension = output_info.dimension;
	const int NUM_VERT_PER_TRI(3);
	const int NUM_VERT_PER_QUAD(4);

	const char * indent4 = "    ";
	string grid_element_name = "cubes";
	if (dimension == 2) { grid_element_name = "squares"; };

	VERTEX_INDEX numv = (vertex_coord.size())/dimension;
	VERTEX_INDEX num_tri = (tri_list.size())/NUM_VERT_PER_TRI;
	VERTEX_INDEX num_quad = (quad_list.size())/NUM_VERT_PER_QUAD;
	VERTEX_INDEX num_grid_cubes = dualiso_info.grid.num_cubes;
	VERTEX_INDEX num_non_empty_cubes = dualiso_info.scalar.num_non_empty_cubes;

	float percent = 0.0;
	if (num_grid_cubes > 0)
	{ percent = float(num_non_empty_cubes)/float(num_grid_cubes); }
	int ipercent = int(100*percent);
	cout << "  Isovalue " << output_info.isovalue[0] << ".  " 
		<< numv << " isosurface vertices.  " << endl;
	cout << "    " << num_tri << " isosurface triangles.  "
		<< num_quad << " isosurface quadrilaterals."
		<< endl;

	if (output_info.allow_multiple_iso_vertices) {
		cout << indent4 << "# cubes with single isosurface vertex: "
			<< dualiso_info.multi_isov.num_cubes_single_isov << endl;
		cout << indent4 << "# cubes with multiple isosurface vertices: "
			<< dualiso_info.multi_isov.num_cubes_multi_isov << endl;

		if (output_info.flag_split_non_manifold) {
			cout << indent4 << "# cubes changed to 2 iso vertices to avoid non-manifold edges: "
				<< dualiso_info.multi_isov.num_non_manifold_split << endl;
		}

		if (output_info.flag_select_split) {
			cout << indent4 << "# cubes changed in selecting isosurface patch splits: "
				<< dualiso_info.multi_isov.num_1_2_change << endl;
		}

	}
}

// **************************************************
// REPORT TIMING INFORMATION
// **************************************************

void QDUAL::report_dualiso_time
	(const IO_INFO & io_info, const DUALISO_TIME & dualiso_time, 
	const char * mesh_type_string)
{
	cout << "CPU time to run ijkdual: " 
		<< dualiso_time.total << " seconds." << endl;
	cout << "    Time to extract " << mesh_type_string << " triangles: "
		<< dualiso_time.extract << " seconds." << endl;
	cout << "    Time to merge identical "
		<< mesh_type_string << " vertices: " 
		<< dualiso_time.merge << " seconds." << endl;
	cout << "    Time to position "
		<< mesh_type_string << " vertices: "
		<< dualiso_time.position << " seconds." << endl;
}


void QDUAL::report_time
	(const IO_INFO & io_info, const IO_TIME & io_time, 
	const DUALISO_TIME & dualiso_time, const double total_elapsed_time)
{
	const char * ISOSURFACE_STRING = "isosurface";
	const char * INTERVAL_VOLUME_STRING = "interval volume";
	const char * mesh_type_string = NULL;

	mesh_type_string = ISOSURFACE_STRING;

	cout << "Time to read file " << io_info.input_filename << ": "
		<< io_time.read_nrrd_time << " seconds." << endl;

	report_dualiso_time(io_info, dualiso_time, mesh_type_string);
	if (!io_info.nowrite_flag) {
		cout << "Time to write "
			<< mesh_type_string << ": " 
			<< io_time.write_time << " seconds." << endl;
	};
	cout << "Total elapsed time: " << total_elapsed_time
		<< " seconds." << endl;
}

// **************************************************
// USAGE/HELP MESSAGES
// **************************************************

// local namespace
namespace {

	void usage_msg(std::ostream & out)
	{
		out << "Usage: qdual [OPTIONS] {isovalue1 isovalue2 ...} {input filename}" << endl;
	}

	void options_msg()
	{
		cerr << "OPTIONS:" << endl;
		cerr << "[-subsample S] [-supersample S] [-position {centroid|cube_center}]" << endl;
		cerr << "[-sep_neg | -sep_pos] [-select_split]" << endl;
		cerr << "[-trimesh | -uniform_trimesh]" << endl;
		cerr <<	"[-no_collapse | -qt_mesh]" << endl;
		cerr << "[-move_vertex]"<<endl;
		cerr <<" [-no_res_AB | -no_res_B | -no_res_C]"<< endl;
		cerr <<" [-collapse_info]"<<endl;
		cerr << "[-off|-iv] [-o {output_filename}] [-stdout]" 
			<< endl; 
		cerr << "[-help] [-s] [-nowrite] [-time]" << endl;
	}

}

void QDUAL::usage_error()
{
	usage_msg(cerr);
	options_msg();
	exit(10);
}

void QDUAL::help()
{
	usage_msg(cout);
	cout << endl;
	cout << "qdual - Quality dual contouring isosurface generation algorithm." 
		<< endl;
	cout << endl;
	cout << "OPTIONS:" << endl;

	cout << "  -subsample S: Subsample grid at every S vertices." << endl;
	cout << "                S must be an integer greater than 1." << endl;
	cout << "  -position {centroid|cube_center}: Isosurface vertex position method." << endl;
	cout << "            centroid: Position isosurface vertices at centroid of"
		<< endl;
	cout << "                      intersection of grid edges and isosurface."
		<< endl;
	cout << "            cube_center: Position isosurface vertices at cube centers." << endl;
	cout << "  -sep_neg:     Isosurface patches separate negative grid vertices."
		<< endl;
	cout << "  -sep_pos:     Isosurface patches separate positive grid vertices."
		<< endl;
	cout << "  -select_split: Select which cube has configuration of split" 
		<< endl
		<< "     isosurface vertices where adjacent cubes share an ambiguous facet."
		<< endl;
	cout <<	"  -collapse_info  Output collapse info."<< endl;
	cout << "  -trimesh:     Output isosurface mesh of triangles." << endl;
	cout << "  -qt_mesh:	 Output quad + tri mesh, degen quads are triangulated." << endl;
	cout <<	"  -no_res_AB , -no_res_B ,-no_res_C " <<endl;
	cout << "                (Allowed only with 3D scalar data.)" << endl;
	cout << "  -no_collapse: No collapses. " << endl;
	cout << "  -off: Output in geomview OFF format. (Default.)" << endl;
	cout << "  -iv: Output in OpenInventor .iv format." << endl;
	cout << "  -o {output_filename}: Write isosurface to file {output_filename}." << endl;
	cout << "  -stdout: Write isosurface to standard output." << endl;
	cout << "  -nowrite: Don't write isosurface." << endl;
	cout << "  -time: Output running time." << endl;
	cout << "  -s: Silent mode." << endl;
	cout << "  -help: Print this help message." << endl;
	exit(20);
}


// **************************************************
// CLASS IO_INFO
// **************************************************

/// IO information
void QDUAL::IO_INFO::Init()
{
	isovalue.clear();
	isovalue_string.clear();
	input_filename = NULL;
	output_filename = NULL;
	isotable_directory = "";
	output_format = OFF;
	report_time_flag = false;
	use_stdout = false;
	nowrite_flag = false;
	flag_silent = false;
	flag_subsample = false;
	subsample_resolution = 2;
	flag_supersample = false;
	supersample_resolution = 2;
	flag_color_alternating = false;  // color simplices in alternating cubes
	region_length = 1;
}

// **************************************************
// class OUTPUT_INFO
// **************************************************

void QDUAL::OUTPUT_INFO::Init()
{
	output_filename = "";
	dimension = 3;
	num_vertices_per_isopoly = 4;
	isovalue[0] = 0;
	isovalue[1] = 0;
	nowrite_flag = false;
	use_stdout = false;
	flag_silent = false;
	output_format = OFF;
	grow_factor = 1;
	shrink_factor = 1;
	grid_spacing.resize(3,1);
}

namespace {

	void split_string(const string & s, const char c,
		string & prefix, string & suffix)
		// split string at last occurrence of character c into prefix and suffix
	{
		string::size_type i = s.rfind(c);
		if (i == string::npos) {
			prefix = s;
			suffix = "";
		}
		else {
			if (i > 0) { prefix = s.substr(0,i); }
			else { prefix = ""; };

			if (i+1 < s.length()) { suffix = s.substr(i+1, s.length()-i-1); }
			else { suffix = ""; };
		}
	}

	string construct_output_filename
		(const IO_INFO & io_info, const int i)
	{
		string prefix, suffix;

		// create output filename
		string fname = string(io_info.input_filename);

#ifndef _WIN32
		// remove path from file name
		split_string(fname, PATH_DELIMITER, prefix, suffix);
		if (suffix != "") { fname = suffix; }
#endif

		string ofilename;

		// construct output filename
		split_string(fname, '.', prefix, suffix);
		if (suffix == "nrrd" || suffix == "nhdr") { ofilename = prefix; }
		else { ofilename = string(io_info.input_filename); }

		ofilename += string(".") + string("isov=") + io_info.isovalue_string[i];

		switch (io_info.output_format) {
		case OFF: 
			ofilename += ".off";
			break;

		case IV:
			ofilename += ".iv";
			break;
		}

		return(ofilename);
	}

}

// **************************************************
// SET ROUTINES
// **************************************************

void QDUAL::set_dualiso_data
	(const IO_INFO & io_info, DUALISO_DATA & dualiso_data, DUALISO_TIME & dualiso_time)
{
	PROCEDURE_ERROR error("set_dualiso_data");

	if (!dualiso_data.IsScalarGridSet()) {
		error.AddMessage
			("Programming error. Scalar field must be set before set_dualiso_data is called.");
		throw error;
	}

	dualiso_data.Set(io_info);
}

void QDUAL::set_output_info
	(const IO_INFO & io_info, 
	const int i, OUTPUT_INFO & output_info)
{
	output_info.nowrite_flag = io_info.nowrite_flag;
	output_info.use_stdout = io_info.use_stdout;
	output_info.flag_silent = io_info.flag_silent;
	output_info.use_triangle_mesh = io_info.use_triangle_mesh;

	if (output_info.use_triangle_mesh) {
		output_info.num_vertices_per_isopoly = 3;
	}

	output_info.grow_factor = 1;
	if (io_info.flag_subsample) 
	{ output_info.grow_factor = io_info.subsample_resolution; }

	output_info.shrink_factor = 1;
	if (io_info.flag_supersample) 
	{ output_info.shrink_factor = io_info.supersample_resolution; }

	output_info.grid_spacing.clear();
	output_info.grid_spacing.resize(io_info.grid_spacing.size());
	for (unsigned int j = 0; j < io_info.grid_spacing.size(); j++)
	{ output_info.grid_spacing[j] = io_info.grid_spacing[j]; }

	output_info.output_format = io_info.output_format;
	output_info.isovalue[0] = io_info.isovalue[i];
	if (i+1 < int(io_info.isovalue.size())) 
	{ output_info.isovalue[1] = io_info.isovalue[i+1]; };

	if (io_info.output_filename != NULL) {
		output_info.output_filename = string(io_info.output_filename);
	}
	else {
		output_info.output_filename = 
			construct_output_filename(io_info, i);
	}

	output_info.Set(io_info);
}

void QDUAL::set_color_alternating
	(const DUALISO_GRID & grid, const vector<VERTEX_INDEX> & cube_list, 
	COLOR_TYPE * color)
{
	const int dimension = grid.Dimension();
	IJK::ARRAY<GRID_COORD_TYPE> coord(dimension);

	const COLOR_TYPE red[3] = { 1.0, 0.0, 0.0 };
	const COLOR_TYPE blue[3] = { 0.0, 0.0, 1.0 };
	const COLOR_TYPE green[3] = { 0.0, 1.0, 0.0 };

	VERTEX_INDEX icube = 0;
	int parity = 0;
	COLOR_TYPE * color_ptr = color;
	for (unsigned int i = 0; i < cube_list.size(); i++) {
		int new_cube = cube_list[i];
		if (icube != new_cube) {
			icube = new_cube;
			grid.ComputeCoord(icube, coord.Ptr());
			int sum = 0;
			for (int d = 0; d < dimension; d++) 
			{ sum += coord[d]; }
			parity = sum%2;
		}

		if (parity == 0) 
		{ std::copy(red, red+3, color_ptr); }
		else
		{ std::copy(blue, blue+3, color_ptr); }

		// set opacity
		color_ptr[3] = 1.0;
		color_ptr += 4;
	}

}

