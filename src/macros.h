/*
 * Copyright 2008-2010, Hans Rullgard, Stockholm University and 
 * Lars-Goran Ofverstedt, Karolinska Institute
 *
 * This file is part of TEM Simulator.
 *
 * TEM Simulator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TEM Simulator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TEM Simulator.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef MACROS_HEADER
#define MACROS_HEADER

/* Version of the program */
#define VERSION_NUMBER "1.3"

/* The maximum length of lines (including \n) in input files */
#define MAX_LINE_LENGTH 201

/* The maximum length of names of parameters used in input files */
#define MAX_NAME_LENGTH 50

/* Character used for comments in input and output files */
#define COMMENT_CHAR '#'

/* Strings to use for input of boolean parameters */
#define YES_STRING "yes"
#define NO_STRING "no"

/* Function used to report warnings. This should be a printf-like function, taking a 
 * format string as first argument, followed by a variable number of arguments determined
 * by the format string. Defined as a macro so it is easy to include information about file
 * and line number in the warning by changing the definition. */
#define WARNING write_log_comment_disp


/* The following macros are used to fix certaing problems in Windows, mostly annoying compiler warnings.
 * FOPEN is a macro for fopen or fopen_s.
 * STRNCPY is a macro for strncpy or strncpy_s (only used inside the function strncpy_term, which
 *         works like strncpy, except that a terminating '\0' is always appended).
 * SNPRINTF is a macro for snprintf or sprintf_s.
 * VSNPRINTF is a macro for vsnprintf or vsnprintf_s. The vsnprintf_s also adds a terminating '\0'.
 * VSNPRINTFCOUNT is used to figure out how much memory to allocate for VSNPRINTF. In the Unix variant
 * this is done by calling vsnprintf with NULL and 0 as the first two arguments. Since this does not
 * work in Windows, the Windows version just makes an estimate which should be sufficient most of the
 * time. If the estimate would be to small, a subsequent call to VSNPRINTF will result in a truncated
 * string, not a buffer overflow.
 * STRTOK is a macro for strtok or strtok_win. strtok_win is a function which does the same thing as 
 * strtok, but without giving rise to compiler warnings in Windows.
 * STRTOK_R is a macro for strtok_r or strtok_s
 * BESSEL0 is the zero order bessel function.
 */

#ifdef _WIN32

#define USE_WINDOWS_MACROS

#define _USE_MATH_DEFINES

#define FOPEN(filepointer,filename,mode) fopen_s(&(filepointer), filename, mode)

#define STRNCPY(dest,destsize,source,num) strncpy_s(dest, destsize, source, num)

#define SNPRINTF sprintf_s

#define VSNPRINTF vsnprintf_win

#define VSNPRINTFCOUNT(format,vararg) (MAX_LINE_LENGTH + strlen(format))

#define STRTOK strtok_win

#define STRTOK_R strtok_s

#define BESSEL0 _j0

#else

#define FOPEN(filepointer,filename,mode) filepointer = fopen(filename, mode)

#define STRNCPY(dest,destsize,source,num) strncpy(dest, source, num)

#define SNPRINTF snprintf

#define VSNPRINTF(string,size,format,vararg) vsnprintf(string, size, format, vararg)

#define VSNPRINTFCOUNT(format,vararg) vsnprintf(NULL, 0, format, vararg)

#define STRTOK strtok

#define STRTOK_R strtok_r

#define BESSEL0 j0

#endif


/* The macro M_PI from math.h is needed below. The include must be after
 * _USE_MATH_DEFINES is defined if this should work in Windows. */
#include <math.h> 


/* Units used for certain input and output */
#define DEFOCUS_UNIT          ONE_MICROMETER
#define CS_UNIT               ONE_MILLIMETER
#define CC_UNIT               ONE_MILLIMETER
#define DET_PIXEL_UNIT        ONE_MICROMETER
#define FOCAL_LENGTH_UNIT     ONE_MILLIMETER
#define APERTURE_UNIT         ONE_MICROMETER
#define ANGLE_UNIT            ONE_DEGREE
#define COND_AP_ANGLE_UNIT    ONE_MILLIRADIAN
#define ACC_ENERGY_UNIT       ONE_KILOELECTRONVOLT
#define ENERGY_SPREAD_UNIT    ONE_ELECTRONVOLT
#define POTENTIAL_UNIT        ONE_VOLT
#define INT_FILE_POT_UNIT     ONE_MILLIVOLT
#define VOL_CONC_UNIT         ONE_PER_MICROMETER_CU
#define SURF_CONC_UNIT        ONE_PER_MICROMETER_SQ


/* Physical constants */
#define ELECTRON_MASS         (9.10938188e-31 * ONE_KILOGRAM)
#define ELEMENTARY_CHARGE     (1.60217646e-19 * ONE_COULOMB)
#define SPEED_OF_LIGHT        (299792458.0 * ONE_METER / ONE_SECOND)
#define PLANCKS_CONSTANT      (6.626068e-34 * ONE_METER * ONE_METER * ONE_KILOGRAM / ONE_SECOND)
#define PLANCKS_CONSTANT_SQ   (PLANCKS_CONSTANT * PLANCKS_CONSTANT)
#define ELEC_REST_ENERGY      (510.99891 * ONE_KILOELECTRONVOLT)
#define BOHR_RADIUS           (0.529177 * ONE_ANGSTROM)
#define ICE_DENS              (33.5 / ONE_NANOMETER_CU)
#define CARBON_DENS           (1e2 / ONE_NANOMETER_CU)

/* Basic units */
#define ONE_KILOGRAM          (1.0)
#define ONE_SECOND            (1.0)
#define ONE_NANOMETER         (1.0)
#define ONE_RADIAN            (1.0)
#define ONE_VOLT              (1.0)

/* Derived units */
#define ONE_COULOMB           (ONE_KILOGRAM * ONE_METER * ONE_METER / (ONE_SECOND * ONE_SECOND * ONE_VOLT))
#define ONE_DEGREE            (M_PI/180.0 * ONE_RADIAN)
#define ONE_MILLIRADIAN       (1e-3 * ONE_RADIAN)
#define ONE_KILOVOLT          (1e3 * ONE_VOLT)
#define ONE_MILLIVOLT         (1e-3 * ONE_VOLT)
#define ONE_KILOELECTRONVOLT  (ELEMENTARY_CHARGE * ONE_KILOVOLT)
#define ONE_ELECTRONVOLT      (ELEMENTARY_CHARGE * ONE_VOLT)
#define ONE_ANGSTROM_SQ       (ONE_ANGSTROM * ONE_ANGSTROM)
#define ONE_ANGSTROM_CU       (ONE_ANGSTROM * ONE_ANGSTROM_SQ)
#define ONE_NANOMETER_SQ      (ONE_NANOMETER * ONE_NANOMETER)
#define ONE_NANOMETER_CU      (ONE_NANOMETER * ONE_NANOMETER_SQ)
#define ONE_PER_MICROMETER_SQ (1.0/(ONE_MICROMETER * ONE_MICROMETER))
#define ONE_PER_MICROMETER_CU (ONE_PER_MICROMETER_SQ / ONE_MICROMETER)
#define ONE_ANGSTROM          (0.1 * ONE_NANOMETER)
#define ONE_MICROMETER        (1e3 * ONE_NANOMETER)
#define ONE_MILLIMETER        (1e6 * ONE_NANOMETER)
#define ONE_METER             (1e9 * ONE_NANOMETER)


/* Component types */
#define TYPE_DETECTOR             "detector"
#define TYPE_ELECTRONBEAM         "electronbeam"
#define TYPE_GEOMETRY             "geometry"
#define TYPE_OPTICS               "optics"
#define TYPE_PARTICLE             "particle"
#define TYPE_PARTICLESET          "particleset"
#define TYPE_SAMPLE               "sample"
#define TYPE_SIMULATION           "simulation"
#define TYPE_VOLUME               "volume"


/* Names of parameters used in input */
#define PAR_ACC_VOLTAGE             "acc_voltage"
#define PAR_ADD_HYDROGEN            "add_hydrogen"
#define PAR_ALIGN_ERR_ROT           "align_err_rot"
#define PAR_ALIGN_ERR_TR            "align_err_tr"
#define PAR_APERTURE                "aperture"
#define PAR_CC                      "cc"
#define PAR_COND_AP_ANGLE           "cond_ap_angle"
#define PAR_CONTRAST_IM             "contrast_im"
#define PAR_CONTRAST_RE             "contrast_re"
#define PAR_COORD_FILE_IN           "coord_file_in"
#define PAR_COORD_FILE_OUT          "coord_file_out"
#define PAR_CS                      "cs"
#define PAR_DEFOCUS_FILE_IN         "defocus_file_in"
#define PAR_DEFOCUS_FILE_OUT        "defocus_file_out"
#define PAR_DEFOCUS_NOMINAL         "defocus_nominal"
#define PAR_DEFOCUS_NONSYST_ERROR   "defocus_nonsyst_error"
#define PAR_DEFOCUS_SYST_ERROR      "defocus_syst_error"
#define PAR_DET_PIX_X               "det_pix_x"
#define PAR_DET_PIX_Y               "det_pix_y"
#define PAR_DIAMETER                "diameter"
#define PAR_DOSE_PER_IM             "dose_per_im"
#define PAR_DOSE_FILE_IN            "dose_file_in"
#define PAR_DOSE_FILE_OUT           "dose_file_out"
#define PAR_DOSE_SD                 "dose_sd"
#define PAR_DQE                     "dqe"
#define PAR_ENERGY_SPREAD           "energy_spread"
#define PAR_ERROR_FILE_IN           "error_file_in"
#define PAR_ERROR_FILE_OUT          "error_file_out"
#define PAR_FAMP                    "famp"
#define PAR_FOCAL_LENGTH            "focal_length"
#define PAR_GAIN                    "gain"
#define PAR_GEN_DOSE                "gen_dose"
#define PAR_GEN_RAND_DEFOCUS        "gen_defocus"
#define PAR_GEN_RAND_POSITIONS      "gen_rand_positions"
#define PAR_GEN_TILT_DATA           "gen_tilt_data"
#define PAR_GENERATE_MICROGRAPHS    "generate_micrographs"
#define PAR_GENERATE_PARTICLE_MAPS  "generate_particle_maps"
#define PAR_GENERATE_VOLUMES        "generate_volumes"
#define PAR_GEOM_FILE_IN            "geom_file_in"
#define PAR_GEOM_FILE_OUT           "geom_file_out"
#define PAR_GEOM_ERRORS             "geom_errors"
#define PAR_IMAGE_AXIS_ORDER        "image_axis_order"
#define PAR_IMAGE_FILE_BYTE_ORDER   "image_file_byte_order"
#define PAR_IMAGE_FILE_FORMAT       "image_file_format"
#define PAR_IMAGE_FILE_OUT          "image_file_out"
#define PAR_LOG_FILE                "log_file"
#define PAR_MAGNIFICATION           "magnification"
#define PAR_MAKE_POSITIVE           "make_positive"
#define PAR_MAP_AXIS_ORDER          "map_axis_order"
#define PAR_MAP_FILE_BYTE_ORDER     "map_file_byte_order"
#define PAR_MAP_FILE_FORMAT         "map_file_format"
#define PAR_MAP_FILE_IM_IN          "map_file_im_in"
#define PAR_MAP_FILE_IM_OUT         "map_file_im_out"
#define PAR_MAP_FILE_RE_IN          "map_file_re_in"
#define PAR_MAP_FILE_RE_OUT         "map_file_re_out"
#define PAR_MTF_A                   "mtf_a"
#define PAR_MTF_ALPHA               "mtf_alpha"
#define PAR_MTF_B                   "mtf_b"
#define PAR_MTF_BETA                "mtf_beta"
#define PAR_MTF_C                   "mtf_c"
#define PAR_MTF_P                   "mtf_p"
#define PAR_MTF_Q                   "mtf_q"
#define PAR_NTILTS                  "ntilts"
#define PAR_NUM_PARTICLES           "num_particles"
#define PAR_NX                      "nx"
#define PAR_NY                      "ny"
#define PAR_NZ                      "nz"
#define PAR_OCCUPANCY               "occupancy"
#define PAR_OFFSET_X                "offset_x"
#define PAR_OFFSET_Y                "offset_y"
#define PAR_OFFSET_Z                "offset_z"
#define PAR_PADDING                 "padding"
#define PAR_PARTICLE_CONC           "particle_conc"
#define PAR_PARTICLE_COORDS         "particle_coords"
#define PAR_PARTICLE_DIST           "particle_dist"
#define PAR_PARTICLE_TYPE           "particle_type"
#define PAR_PDB_FILE_IN             "pdb_file_in"
#define PAR_PDB_TRANSF_FILE_IN      "pdb_transf_file_in"
#define PAR_PIXEL_SIZE              "pixel_size"
#define PAR_RAND_SEED               "rand_seed"
#define PAR_SMOOTHNESS              "smoothness"
#define PAR_SOURCE                  "source"
#define PAR_THETA_INCR              "theta_incr"
#define PAR_THETA_START             "theta_start"
#define PAR_THICKNESS_CENTER        "thickness_center"
#define PAR_THICKNESS_EDGE          "thickness_edge"
#define PAR_TILT_AXIS               "tilt_axis"
#define PAR_TILT_ERR                "tilt_err"
#define PAR_TILT_MODE               "tilt_mode"
#define PAR_TOTAL_DOSE              "total_dose"
#define PAR_USE_DEFOCUS_CORR        "use_defocus_corr"
#define PAR_USE_IMAG_POT            "use_imag_pot"
#define PAR_USE_QUANTIZATION        "use_quantization"
#define PAR_VOXEL_SIZE              "voxel_size"
#define PAR_WHERE                   "where"

/* Possible values of some string parameters */
#define PAR_BYTE_ORDER__BE          "be"
#define PAR_BYTE_ORDER__LE          "le"
#define PAR_BYTE_ORDER__NATIVE      "native"
#define PAR_FILE_FORMAT__MRC        "mrc"
#define PAR_FILE_FORMAT__MRC_INT    "mrc-int"
#define PAR_FILE_FORMAT__RAW        "raw"
#define PAR_GEOM_ERRORS__FILE       "file"
#define PAR_GEOM_ERRORS__NONE       "none"
#define PAR_GEOM_ERRORS__RANDOM     "random"
#define PAR_PARTICLE_COORDS__FILE   "file"
#define PAR_PARTICLE_COORDS__GRID   "grid"
#define PAR_PARTICLE_COORDS__RANDOM "random"
#define PAR_SOURCE__MAP             "map"
#define PAR_SOURCE__PDB             "pdb"
#define PAR_SOURCE__RANDOM          "random"
#define PAR_TILT_MODE__SINGLE_PARTICLE "single_particle"
#define PAR_TILT_MODE__TILTSERIES   "tiltseries"
#define PAR_WHERE__BOTTOM           "bottom"
#define PAR_WHERE__SURFACE          "surface"
#define PAR_WHERE__TOP              "top"
#define PAR_WHERE__VOLUME           "volume"


#endif /* MACROS_HEADER */
