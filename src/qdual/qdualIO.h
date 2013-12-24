/// \file qdualIO.h
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

#ifndef _QDUALIO_
#define _QDUALIO_

#include <ctime>
#include <string>

#include "ijk.txx"
#include "ijkgrid_nrrd.txx"

#include "qdual_types.h"
#include "qdual_datastruct.h"


namespace QDUAL {

// **************************************************
// TYPE DEFINITIONS
// **************************************************

  typedef float COLOR_TYPE;           ///< Color type.

  /// Nrrd header.
  typedef IJK::NRRD_DATA<int, QDUAL::AXIS_SIZE_TYPE> NRRD_HEADER; 


// **************************************************
// OUTPUT_FORMAT
// **************************************************

  typedef enum { OFF, IV } OUTPUT_FORMAT;   ///< Output format.

// **************************************************
// IO INFORMATION
// **************************************************

  /// IO information
  class IO_INFO:public DUALISO_DATA_FLAGS {

  protected:
    void Init();

  public:
    SCALAR_ARRAY isovalue;        ///< List of isovalues.
    std::vector<std::string> isovalue_string;
    COORD_ARRAY grid_spacing;
    char * input_filename;
    char * output_filename;
    std::string isotable_directory;
    OUTPUT_FORMAT output_format;
    bool report_time_flag;
    bool use_stdout;
    bool nowrite_flag;
    bool flag_silent;
    bool flag_subsample;
    int subsample_resolution;
    bool flag_supersample;
    int supersample_resolution;
    bool flag_color_alternating;  ///< Color simplices in alternating cubes
    int region_length;

    /// List of high resolution arguments,
    ///   e.g., "-highres {coord list}".
    std::vector<std::string> high_resolution_option;

  public:
    IO_INFO() { Init(); };
    ~IO_INFO() { Init(); };
  };

// **************************************************
// OUTPUT INFORMATION
// **************************************************

  /// Output information.
  class OUTPUT_INFO:public DUALISO_DATA_FLAGS {

  protected:
    void Init();

  public:
    std::string output_filename;
    int dimension;
    int num_vertices_per_isopoly;
    SCALAR_TYPE isovalue[2];
    bool nowrite_flag;
    bool use_stdout;
    bool flag_silent;
    OUTPUT_FORMAT output_format;
    COORD_ARRAY grid_spacing;
    int grow_factor;
    int shrink_factor;

    OUTPUT_INFO() { Init(); };
    ~OUTPUT_INFO() { Init(); };
  };

// **************************************************
// TIMING FUNCTIONS/CLASSES
// **************************************************

  /// Elapsed CPU time.
  class ELAPSED_CPU_TIME {

  protected:
    clock_t t;

  public:
    ELAPSED_CPU_TIME() { t = clock(); };

    clock_t getElapsed() {
      clock_t old_t = t;
      t = clock();
      return(t - old_t);
    };
  };

  /// Elapsed wall time.
  class ELAPSED_TIME {

  protected:
    time_t t;

  public:
    ELAPSED_TIME() { time(&t);  };

    double getElapsed() {
      time_t old_t = t;
      time(&t);
      return(difftime(t,old_t));
    };
  };

  /// IO time.
  struct IO_TIME {
    double read_nrrd_time;  ///< Wall time to read nrrd file.
    double write_time;      ///< Wall time to write output.
  };

// **************************************************
// PARSE COMMAND LINE
// **************************************************

  /// Parse the command line.
  void parse_command_line(int argc, char **argv, IO_INFO & io_info);

  /// Check input information in io_info
  bool check_input
    (const IO_INFO & io_info, 
     const DUALISO_SCALAR_GRID_BASE & scalar_grid,
     IJK::ERROR & error);

// **************************************************
// READ NEARLY RAW RASTER DATA (nrrd) FILE
// **************************************************

  /// Read a nearly raw raster data (nrrd) file.
  void read_nrrd_file
    (const char * input_filename, DUALISO_SCALAR_GRID & scalar_grid, 
     NRRD_HEADER & nrrd_header, IO_TIME & io_time);

// **************************************************
// OUTPUT ISOSURFACE
// **************************************************

  void output_dual_isosurface
    (const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
     const DUAL_ISOSURFACE & dual_isosurface,
     const DUALISO_INFO & dualiso_info, IO_TIME & io_time);

  void output_dual_isosurface
    (const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const DUALISO_INFO & dualiso_info, IO_TIME & io_time);

  /// Output isosurface of triangles.
  void output_dual_tri_isosurface
  (const OUTPUT_INFO & output_info, const DUALISO_DATA & isodual_data,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const DUALISO_INFO & dualiso_info, IO_TIME & io_time);

  void output_dual_isosurface_color
    (const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
     const DUAL_ISOSURFACE & dual_isosurface,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
     const DUALISO_INFO & dualiso_info, IO_TIME & io_time);

  void output_dual_isosurface_color
    (const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
     const DUALISO_INFO & dualiso_info, IO_TIME & io_time);

  void output_dual_isosurface_color_alternating
    (const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
     const DUAL_ISOSURFACE & dual_isosurface,
     const DUALISO_INFO & dualiso_info, IO_TIME & io_time);

// **************************************************
// RESCALE ROUTINES
// **************************************************

  /// Rescale subsampled/supersampled vertex coordinates.
  /// Also rescale to reflect grid spacing.
  void rescale_vertex_coord
    (const OUTPUT_INFO & output_info, std::vector<COORD_TYPE> & vertex_coord);

  /// Rescale vertex coordinates by grow and shrink factor and by grid_spacing.
  /// Precondition: grid_spacing.size() equals vertex dimension.
  void rescale_vertex_coord
    (const int grow_factor, const int shrink_factor,
     const COORD_ARRAY & grid_spacing, COORD_ARRAY & vertex_coord);

// **************************************************
// WRITE_DUAL_MESH
// **************************************************

  void write_dual_mesh
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist);

  void write_dual_mesh
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     IO_TIME & io_time);

  void write_dual_mesh_color
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  void write_dual_mesh_color
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist, 
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
     IO_TIME & io_time);

  /// Write dual isosurface triangular mesh.
  /// @param output_info Output information.
  /// @param vertex_coord List of vertex coordinates.
  /// @param tri_vert[] List of triangle vertices.
  ///        tri_vert[3*i+k] is k'th vertex of triangle i.
  void write_dual_tri_mesh
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert);

  /// Write dual isosurface triangular mesh.
  /// Record write time.
  void write_dual_tri_mesh
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert,
   IO_TIME & io_time);


// **************************************************
// SET ROUTINES
// **************************************************

  /// Set dualiso_data based on io_info.
  /// Precondition: Scalar field in dualiso_data must be set before
  ///   this routines is called.
  void set_dualiso_data
    (const IO_INFO & io_info, DUALISO_DATA & dualiso_data, DUALISO_TIME & dualiso_time);

  /// Set output_info based on isotable, io_info and isovalue index i.
  void set_output_info
    (const IO_INFO & io_info, 
     const int i, OUTPUT_INFO & output_info);

  /// Set simplices in alternating cubes to have different colors.
  void set_color_alternating
    (const DUALISO_GRID & grid, const std::vector<VERTEX_INDEX> & cube_list, 
     COLOR_TYPE * color);

// **************************************************
// REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
// **************************************************

  void report_num_cubes
    (const DUALISO_GRID & full_grid, const IO_INFO & io_info, 
     const DUALISO_DATA & dualiso_data);

  void report_iso_info
    (const OUTPUT_INFO & output_info, const DUALISO_DATA & dualiso_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist, 
     const DUALISO_INFO & dualiso_info);

// **************************************************
// REPORT TIMING INFORMATION
// **************************************************

  void report_dualiso_time
    (const IO_INFO & io_info, const DUALISO_TIME & dualiso_time, 
     const char * mesh_type_string);

  void report_time
    (const IO_INFO & io_info, const IO_TIME & io_time, 
     const DUALISO_TIME & dualiso_time, const double total_elapsed_time);

// **************************************************
// USAGE/HELP MESSAGES
// **************************************************

  void usage_error();
  void help();
}

#endif
