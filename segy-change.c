/*
 Source code of:  segy-change.c

 Copyright (C) 2009 Giuseppe Stanghellini (1), Gabriela Carrara (2).
 (1) Istituto di Scienze Marine, Geologia Marina, CNR, Bologna, Italy
 (2) LNEG – Marine Geology Department, Portugal

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.


 INSTRUCTIONS ON HOW TO COMPILE THE PROGRAM:

 gcc -lm -o segy-change segy-change.c

 to have instructions on its usage execute it without arguments.

 */

#define _FILE_OFFSET_BITS 64
// #define WITH_SDL moved into makefile
#define WITH_OPENMP
#define _XOPEN_SOURCE
#define _DEFAULT_SOURCE

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <arpa/inet.h>
#include <ctype.h>

#ifdef WITH_OPENMP
#include "omp.h"
#endif

/* bool definitions.
 */
#ifndef __cplusplus
typedef short bool;
#define true 1
#define false 0
#endif

#ifdef WITH_SDL
#include "SDL2/SDL.h"
#include "SDL2/SDL2_gfxPrimitives.h"
void DisplayData_Init();
void DisplayData_DrawWaves();
void DisplayData_CheckZoomAndShift(bool adapt);
#endif

#define my_min(a, b) ((a) < (b) ? (a) : (b))
#define my_max(a, b) ((a) > (b) ? (a) : (b))
#define my_sgn(a) ((a) < 0 ? (-1) : ((a) > 0 ? (1) : (0)))

void prerror_and_exit(const char *fmt, ...) {
	va_list args;

	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	exit(-1);
}

char *trim(char *s) // Trim blanks from both sides of string.
{
	int i = 0;
	int j = strlen(s) - 1;
	int k = 0;

	while (isspace(s[i]) && s[i] != '\0')
		i++;

	while (isspace(s[j]) && j >= 0)
		j--;

	while (i <= j)
		s[k++] = s[i++];

	s[k] = '\0';

	return s;
}

void remove_parms(int *argc, char **argv, int j, int num_args)
{
	int i;
	for(i = j + num_args; i < *argc; i++)
		argv[i - num_args] = argv[i];
	*argc -= num_args;
}

/* TAKE BACK PARAMETERS FROM ARGLIST
 */
int take_parm(int argc, char **argv, char *parm, int num_args) {
	int i;
	for (i = 0; i < argc; i++)
		if (strcmp(argv[i], parm) == 0)
		{
			if(argc > i + num_args)
			{
				for(int j = i + 1; j < i + 1 + num_args; j++)
					if(argv[j][0] == '-' && argv[j][1] != 0)
					{
						if(num_args > 1)
							prerror_and_exit("Error: switch %s requires %d additional parameters\n", parm, num_args);
						else
							prerror_and_exit("Error: switch %s requires %d additional parameter\n", parm, num_args);
					}
				return i;
			}
			else
			{
				if(num_args > 1)
					prerror_and_exit("Error: switch %s requires %d additional parameters\n", parm, num_args);
				else
					prerror_and_exit("Error: switch %s requires %d additional parameter\n", parm, num_args);
			}
		}
	return 0;

}

/* Extract a field from a string, given a separator,
 * the first field is the number one.
 */
char *get_field(char *str, int field_nr, char *buffer, char sep) {
	int i, start, end, count;
	int l = strlen(str);

	count = start = 0;
	end = -1;
	for (i = 0; i < l; i++) {
		if (str[i] == sep) {
			count++;
			if (count == field_nr) {
				end = i;
				break;
			} else if (count < field_nr)
				start = i + 1;
		}
	}

	if (end == -1 && count == field_nr - 1)
		end = l + 1;
	if (end > start) {
		strncpy(buffer, str + start, end - start);
		buffer[end - start] = 0;
		return buffer;
	}
	buffer[0] = 0;
	return NULL;
}

#define GET_IN_SAMPLE(i) get_val(&segy_file, (&segy_file)->trace_data, i)
#define SET_IN_SAMPLE(v, i) set_val(v, &segy_file, (&segy_file)->trace_data, i)

#define SET_OUT_SAMPLE(v, i) set_val(v, &out_segy_file, (&out_segy_file)->trace_data, i)
#define GET_OUT_SAMPLE(i) get_val(&out_segy_file, (&out_segy_file)->trace_data, i)

#define GET_SEGY_TRACE() get_segy_trace(&segy_file, verbose)

/* Arrays containing offset and types of values stored inside the segy and traces headers, used to perform the flipping of the endianess */
int segy_header_types[] = { 3200, 'I', 3204, 'I', 3208, 'I', 3212, 'S', 3214,
		'S', 3216, 'S', 3218, 'S', 3220, 'S', 3222, 'S', 3224, 'S', 3226, 'S',
		3228, 'S', 3230, 'S', 3232, 'S', 3234, 'S', 3236, 'S', 3238, 'S', 3240,
		'S', 3242, 'S', 3244, 'S', 3246, 'S', 3248, 'S', 3250, 'S', 3252, 'S',
		3254, 'S', 3256, 'S', 3258, 'S', -1, -1 };
char *segy_header_names[] = { "JOB_IDENTIFICATION_NUMBER", "LINE_NUMBER",
		"REEL_NUMBER", "NUMBER_OF_DATA_TRACES_PER_RECORD",
		"NUMBER_OF_AUXILLARY_TRACES_PER_RECORD",
		"SAMPLE_INTERVAL_FOR_THIS_REEL_MICROSECONDS",
		"SAMPLE_INTERVAL_FOR_ORIGINAL_FIELD_RECORDING_MICROSECONDS",
		"NUMBER_OF_SAMPLES_PER_DATA_TRACE_FOR_THIS_REEL",
		"NUMBER_OF_SAMPLES_PER_DATA_TRACE_ORIGINAL_FIELD_RECORDING",
		"DATA_SAMPLE_FORMAT_CODE", "NOMINAL_CDP_FOLD", "TRACE_SORTING_CODE",
		"NUMBER_OF_VERTICALLY_SUMMED_TRACES", "SWEEP_FREQUENCY_AT_START_HZ",
		"SWEEP_FREQUENCY_AT_END_HZ", "SWEEP_LENGTH_MILLISECONDS", "SWEEP_TYPE",
		"TRACE_NUMBER_OF_SWEEP_CHANNEL",
		"SWEEP_TAPER_LENGTH_AT_START_MILLISECONDS",
		"SWEEP_TAPER_LENGTH_AT_END_MILLISECONDS", "TAPER_TYPE",
		"CORRELATED_DATA_TRACES", "BINARY_GAIN_RECOVERED",
		"AMPLITUDE_RECOVERY_METHOD", "MEASUREMENT_SYSTEM", "IMPULSE_SIGNAL",
		"VIBRATORY_POLARITY_CODE", NULL };
int trace_header_types[] = { 0, 'I', 4, 'I', 8, 'I', 12, 'I', 16, 'I', 20, 'I',
		24, 'I', 28, 'S', 30, 'S', 32, 'S', 34, 'S', 36, 'I', 40, 'I', 44, 'I',
		48, 'I', 52, 'I', 56, 'I', 60, 'I', 64, 'I', 68, 'S', 70, 'S', 72, 'I',
		76, 'I', 80, 'I', 84, 'I', 88, 'S', 90, 'S', 92, 'S', 94, 'S', 96, 'S',
		98, 'S', 100, 'S', 102, 'S', 104, 'S', 106, 'S', 108, 'S', 110, 'S',
		112, 'S', 114, 'S', 116, 'S', 118, 'S', 120, 'S', 122, 'S', 124, 'S',
		126, 'S', 128, 'S', 130, 'S', 132, 'S', 134, 'S', 136, 'S', 138, 'S',
		140, 'S', 142, 'S', 144, 'S', 146, 'S', 148, 'S', 150, 'S', 152, 'S',
		154, 'S', 156, 'S', 158, 'S', 160, 'S', 162, 'S', 164, 'S', 166, 'S',
		168, 'S', 170, 'S', 172, 'S', 174, 'S', 176, 'S', 178, 'S', 180, 'U',
		186, 'U', 194, 'U', 198, 'S', 200, 'I', 204, 'S', 206, 'S', 208, 'S',
		210, 'S', 212, 'S', 214, 'S', 216, 'S', 218, 'S', 220, 'F', 224, 'S',
		226, 'S', 228, 'I', 232, 'I', 236, 'I', -1, -1 };
char *trace_header_names[] = { "TRACE_SEQUENCE_NUMBER_WITHIN_LINE",
		"TRACE_SEQUENCE_NUMBER_WITHIN_REEL", "ORIGINAL_FIELD_RECORD_NUMBER",
		"TRACE_NUMBER_WITHIN_FIELD_RECORD", "SOURCE_POINT_NUMBER", "CDP_NUMBER",
		"CDP_SEQUECE_NUMBER", "TRACE_IDENTIFICATION_CODE",
		"NUMBER_OF_VERTICALLY_SUMMED_TRACES",
		"NUMBER_OF_HORIZONTALLY_SUMMED_TRACES_FOLD", "DATA_USE",
		"SOURCE_RECEIVER_OFFSET_IN_FEET_OR_METERS",
		"RECEIVER_GROUP_ELEVATION_IN_FEET_OR_METERS",
		"SURFACE_ELEVATION_AT_SOURCE_FEET_OR_METERS",
		"SOURCE_DEPTH_BELOW_SURFACE", "DATUM_ELEVATION_AT_RECEIVER_GROUP",
		"DATUM_ELEVATION_AT_SOURCE", "WATER_DEPTH_AT_SOURCE",
		"WATER_DEPTH_AT_RECEIVER_GROUP",
		"ELEVATION_MULTIPLICATION_SCALAR_FOR_BYTES",
		"COORDINATE_MULTIPLICATION_SCALAR_FOR_BYTES_73_88",
		"SOURCE_X_FEET_OR_METERS_OR_LONGITUDE",
		"SOURCE_Y_FEET_OR_METERS_OR_LATITUDE",
		"RECEIVER_X_FEET_OR_METERS_OR_LONGITUDE",
		"RECEIVER_Y_FEET_OR_METERS_OR_LATITUDE", "COORDINATE_UNITS",
		"WEATHERING_VELOCITY", "SUB_WEATHERING_VELOCITY",
		"UPHOLE_TIME_AT_SOURCE_MILLISECONDS",
		"UPHOLE_TIME_AT_GROUP_MILLISECONDS",
		"SOURCE_STATIC_CORRECTION_MILLISECONDS",
		"RECEIVER_STATIC_CORRECTION_MILLISECONDS", "TOTAL_STATIC_APPLIED",
		"LAG_TIME_A_BETWEEN_TRACE_HEADER_TIME_AND_TIME",
		"LAG_TIME_B_BETWEEN_TIME_BREAK_AND_SOURCE_TIME",
		"DELAY_TIME_BETWEEN_SOURCE_AND_RECORDING_TIME",
		"BRUTE_START_TIME_MILLISECONDS", "MUTE_END_TIME_MILLISECONDS",
		"NUMBER_OF_SAMPLES_IN_THIS_TRACE", "SAMPLE_INTERVAL_MICROSECONDS",
		"GAIN_TYPE_1__FIXED_2__BINARY", "INSTRUMENT_GAIN_CONSTANT",
		"INSTRUMENT_EARLY_OR_INITIAL_GAIN", "CORRELATED_1__YES_2__NO",
		"SWEEP_FREQUENCY_AT_START_HZ", "SWEEP_FREQUENCY_AT_END_HZ",
		"SWEEP_LENGTH_MILLISECONDS", "SWEEP_TYPE_1__LINEAR_2__PARABOLIC",
		"SWEEP_TAPER_LENGTH_AT_START_MILLISECONDS",
		"SWEEP_TAPER_LENGTH_AT_END_MILLISECONDS", "SWEEP_TAPER_TYPE_1__LINEAR",
		"ALIAS_FILTER_FREQUENCY_HZ", "ALIAS_FILTER_SLOPE_DBOCTAVE",
		"NOTCH_FILTER_FREQUENCY_HZ", "NOTCH_FILTER_SLOPE_DBOCTAVE",
		"LOW_CUT_FREQUENCY_HZ", "HIGH_CUT_FREQUENCY_HZ",
		"LOW_CUT_FILTER_SLOPE_DBOCTAVE", "HIGH_CUT_FILTER_SLOPE_DBOCTAVE",
		"YEAR_DATA_RECORDED", "DAY_OF_YEAR", "HOUR_OF_DAY_24_HOUR_CLOCK",
		"MINUTE_OF_HOUR", "SECOND_OF_MINUTE_FOR_TRACE_START",
		"TIME_BASIS_CODE_1LOCAL_2GMT_3OTHER",
		"TRACE_WEIGHTING_FACTOR_2N_VOLTS_FOR_LEAST",
		"RECEIVER_GROUP_NUMBER_AT_ROLL_SWITCH_POSITION_1",
		"RECEIVER_GROUP_NUMBER_FOR_FIRST_TRACE_IN_FIELD_RECORD",
		"RECEIVER_GROUP_NUMBER_FOR_LAST_TRACE_IN_FIELD_RECORD",
		"GAP_SIZE_NUMBER_OF_RECEIVER_GROUPS_DROPPED",
		"OVERTRAVEL_ASSOCIATED_WITH_TAPER_AT_START_OR_END", NULL };

const int ascii2ebcdic[256] = { 0, 1, 2, 3, 55, 45, 46, 47, 22, 5, 37, 11, 12,
		13, 14, 15, 16, 17, 18, 19, 60, 61, 50, 38, 24, 25, 63, 39, 28, 29, 30,
		31, 64, 79, 127, 123, 91, 108, 80, 125, 77, 93, 92, 78, 107, 96, 75, 97,
		240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 122, 94, 76, 126, 110,
		111, 124, 193, 194, 195, 196, 197, 198, 199, 200, 201, 209, 210, 211,
		212, 213, 214, 215, 216, 217, 226, 227, 228, 229, 230, 231, 232, 233,
		74, 224, 90, 95, 109, 121, 129, 130, 131, 132, 133, 134, 135, 136, 137,
		145, 146, 147, 148, 149, 150, 151, 152, 153, 162, 163, 164, 165, 166,
		167, 168, 169, 192, 106, 208, 161, 7, 32, 33, 34, 35, 36, 21, 6, 23, 40,
		41, 42, 43, 44, 9, 10, 27, 48, 49, 26, 51, 52, 53, 54, 8, 56, 57, 58,
		59, 4, 20, 62, 225, 65, 66, 67, 68, 69, 70, 71, 72, 73, 81, 82, 83, 84,
		85, 86, 87, 88, 89, 98, 99, 100, 101, 102, 103, 104, 105, 112, 113, 114,
		115, 116, 117, 118, 119, 120, 128, 138, 139, 140, 141, 142, 143, 144,
		154, 155, 156, 157, 158, 159, 160, 170, 171, 172, 173, 174, 175, 176,
		177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190,
		191, 202, 203, 204, 205, 206, 207, 218, 219, 220, 221, 222, 223, 234,
		235, 236, 237, 238, 239, 250, 251, 252, 253, 254, 255 };

const int ebcdic2ascii[256] = { 0, 1, 2, 3, 156, 9, 134, 127, 151, 141, 142, 11,
		12, 13, 14, 15, 16, 17, 18, 19, 157, 133, 8, 135, 24, 25, 146, 143, 28,
		29, 30, 31, 128, 129, 130, 131, 132, 10, 23, 27, 136, 137, 138, 139,
		140, 5, 6, 7, 144, 145, 22, 147, 148, 149, 150, 4, 152, 153, 154, 155,
		20, 21, 158, 26, 32, 160, 161, 162, 163, 164, 165, 166, 167, 168, 91,
		46, 60, 40, 43, 33, 38, 169, 170, 171, 172, 173, 174, 175, 176, 177, 93,
		36, 42, 41, 59, 94, 45, 47, 178, 179, 180, 181, 182, 183, 184, 185, 124,
		44, 37, 95, 62, 63, 186, 187, 188, 189, 190, 191, 192, 193, 194, 96, 58,
		35, 64, 39, 61, 34, 195, 97, 98, 99, 100, 101, 102, 103, 104, 105, 196,
		197, 198, 199, 200, 201, 202, 106, 107, 108, 109, 110, 111, 112, 113,
		114, 203, 204, 205, 206, 207, 208, 209, 126, 115, 116, 117, 118, 119,
		120, 121, 122, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220,
		221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 123, 65, 66, 67,
		68, 69, 70, 71, 72, 73, 232, 233, 234, 235, 236, 237, 125, 74, 75, 76,
		77, 78, 79, 80, 81, 82, 238, 239, 240, 241, 242, 243, 92, 159, 83, 84,
		85, 86, 87, 88, 89, 90, 244, 245, 246, 247, 248, 249, 48, 49, 50, 51,
		52, 53, 54, 55, 56, 57, 250, 251, 252, 253, 254, 255 };

/* SEGY HEADER.
 */
typedef struct {
	unsigned char EBCDIC[3200];
	unsigned char BINARY[400];

#define GET_SEGYH_Line_number(header) get_int((void*)header + 3204)
#define GET_SEGYH_Reel_number(header) get_int((void*)header + 3208)
#define GET_SEGYH_Number_of_samples_per_datatrace_for_this_reel(header) get_short((void*)header + 3220)
#define GET_SEGYH_Sample_interval_for_this_reel(header) get_short((void*)header + 3216)
#define GET_SEGYH_Measurements_system(header) get_short((void*)header + 3254)
#define GET_SEGYH_Number_of_data_traces_per_record(header) get_short((void*)header + 3212)
#define GET_SEGYH_Data_sample_format_code(header) get_short((void*)header + 3224)
} SEGY_header;

/* TRACE HEADER.
 */
typedef struct {
	unsigned char HEADER[240];

#define GET_SEGYTRACEH_Field(header, offset) get_int((void*)header + offset)
#define GET_SEGYTRACEH_Trace_sequence_number_within_line(header) get_int((void*)header)
#define GET_SEGYTRACEH_Trace_sequence_number_within_reel(header) get_int((void*)header + 4)
#define GET_SEGYTRACEH_Original_field_record_number(header) get_int((void*)header + 8)
#define GET_SEGYTRACEH_Trace_number_within_field_record(header) get_int((void*)header + 12)
#define GET_SEGYTRACEH_Number_of_samples_in_this_trace(header) get_short((void*)header + 114)
#define GET_SEGYTRACEH_USHORT_Number_of_samples_in_this_trace(header) get_ushort((void*)header + 114)
#define GET_SEGYTRACEH_Sample_interval(header) get_short((void*)header + 116)
#define GET_SEGYTRACEH_Delay_time_between_source_and_recording_time(header) get_short((void*)header + 108)

/* Coordinate entries.
 */
#define GET_SEGYTRACEH_Source_x_or_lon(header) get_int((void*)header + 72)
#define GET_SEGYTRACEH_Source_y_or_lat(header) get_int((void*)header + 76)
#define GET_SEGYTRACEH_Receiver_x_or_lon(header) get_int((void*)header + 80)
#define GET_SEGYTRACEH_Receiver_y_or_lat(header) get_int((void*)header + 84)
#define GET_SEGYTRACEH_Coord_mult_scalar(header) get_short((void*)header + 70)
#define GET_SEGYTRACEH_Coord_units_1_fm_2_arcsec(header) get_short((void*)header + 88)

} SEGY_trace_header;

typedef struct {
	SEGY_header header;
	SEGY_trace_header trace_header;
	unsigned char *trace_data;
	double *trace_data_double;
	FILE *fp;
	char *fname;
} SEGY_file;

/* PAPER SIZES for postscript plot.
 */
int paper_size_x[] = { 841, 594, 420, 297, 210, 500 };
int paper_size_y[] = { 1189, 841, 594, 420, 297, 1000 };

/* Static variables.
 */
FILE *file_of_fields;
SEGY_file out_segy_file;
SEGY_file segy_file;

bool only_traces_with, replace_ebcdic, plot_data, dump_fields, change_fields,
		dump_header_fields, change_header_fields, shot_renumber, trace_renumber,
		no_header, flip_endianess, add_xy, source_1_or_receiver_2,
		print_rec_seq_num, dump_xy, use_names, scan, apply_correction,
		no_EBCDIC_stamp, enable_X11;
char *fields_to_dump, *fields_to_change_fname, *add_coordinates_fname,
		*header_fields_to_change, *header_fields_to_dump,
		*traces_fields_valid_values;
unsigned char convert_to = ' ', correction_op = '+';
double correction_val;
int _i, _j, _n, i, j, k, all_file, dump, count, vertical_stack;
int initial_record, initial_trace_seq, trace_offset, trace_start, trace_end,
		rec_start, rec_end, output_segy, verbose, n_traces, n_samples, reccnt,
		current_trace, processed_traces, skip_ntraces, only_ntraces,
		skip_nsamples, only_nsamples, current_record, num_traces_per_cm, page_format;
off_t initial_seek;
double trace_scale;
double actual_row;
long actual_record;
char my_ebcdic[3200];
double trace_min_val, trace_max_val;
long total_records, total_traces, min_traces_per_records, max_traces_per_records, min_num_samples, max_num_samples;
long min_num_samples_rec_num, max_num_samples_rec_num;
long min_num_samples_trace_num, max_num_samples_trace_num;
long prev_record_nr;
short delay_time;

typedef struct {
	int original_field_record, trace_seq_within_reel,
			trace_seq_within_field_record;
	double lat_or_y, lon_or_x, elevation, depth;
	int unit_of_measure_1_feet_or_meters_2_arcsec;
} COORDS;

COORDS *coords;
int num_coords;
int scale_x, scale_y;
double max_x, min_x, max_y, min_y, max_z, min_z;

/* FUNCTIONS TO TEST/ADJUST DATA TO CORRECT ENDIANESS.
 */
void swap2(unsigned short* x) {
	*x = ((*x) >> 8) | ((*x) << 8);
}

void swap4(unsigned int* x) {
	*x = ((*x) >> 24) | (((*x) << 8) & 0x00FF0000) | (((*x) >> 8) & 0x0000FF00)
			| ((*x) << 24);
}

/* return the parameter index given its name or -1 if wrong name.
 */
int get_parameter_index_by_name(char **names, char *name) {
	int r = 0;
	while (names[r] != NULL) {
		if (strcmp(names[r], name) == 0) {
			break;
		}
		r++;
	}
	if (names[r] != NULL)
		return r;
	else {
		prerror_and_exit("FATAL ERROR: wrong parameter name.\n");
		return -1;
	}
}

int get_parameter_index_by_offset(int *header_types, int offset) {
	int r = 0;
	while (header_types[r * 2] != offset) {
		r++;
	}
	if (header_types[r * 2] != -1)
		return r;
	else {
		prerror_and_exit("FATAL ERROR: wrong parameter offset.\n");
		return -1;
	}
}

/* Fill the given offsets string with offsets and types of the parameters given into my_names.
 * require also the names variable that must be either header_names or trace_names.
 */
char my_header_fields_offsets[2048];
char my_traces_fields_offsets[2048];
char my_valid_traces_fields_offsets[2048];
void from_names_to_offsets(char **names, char *my_names, char *my_offsets) {
	char n[2048];
	int j = 0, c = 0;
	my_offsets[0] = 0;
	bool next_field_is_value = false;
	bool again = true;
	while (again) {
		if (my_names[j] == ',' || my_names[j] == ':' || my_names[j] == 0) {
			strncpy(n, my_names + c, j - c);
			n[j - c] = 0;
			if (next_field_is_value) {
				strcat(my_offsets, ":");
				strcat(my_offsets, n);
			} else {
				int p = get_parameter_index_by_name(names, n);
				if (p == -1)
					prerror_and_exit("FATAL ERROR: Parameter name not found");
				if (segy_header_names == names) /* the given names are the header_names */
				{
					sprintf(n, "%d:%c", segy_header_types[p * 2],
							segy_header_types[p * 2 + 1]);
					strcat(my_offsets, n);
				} else /* the given names are trace_names */
				{
					sprintf(n, "%d:%c", trace_header_types[p * 2],
							trace_header_types[p * 2 + 1]);
					strcat(my_offsets, n);
				}
			}

			c = j + 1;
			if (my_names[j] == ',')
				strcat(my_offsets, ",");
			if (my_names[j] == ':')
				next_field_is_value = true;
			else
				next_field_is_value = false;
		}
		if (my_names[j] != 0)
			j++;
		else
			again = false;
	}
}

int is_little_endian() {
	int a;
	char *b;
	a = 1;
	b = (char*) &a;
	if (b[0] == 0 && b[sizeof(int) - 1] == 1)
		return 0;
	else
		return 1;
}

short get_short(void *i) {
	short j;
	if (is_little_endian())
		swab(i, (unsigned char*) (&j), sizeof(short));

	else
		j = *((short*) i);
	return j;
}

unsigned short get_ushort(void *i) {
	unsigned short j;
	if (is_little_endian())
		swab(i, (unsigned char*) (&j), sizeof(short));

	else
		j = *((unsigned short*) i);
	return j;
}

int get_int(void *i) {
	unsigned char *_i, *_j;
	int j;
	if (is_little_endian()) {
		_j = (unsigned char*) &j;
		_i = (unsigned char*) i;
		_j[0] = _i[3];
		_j[1] = _i[2];
		_j[2] = _i[1];
		_j[3] = _i[0];
	} else
		j = *((int*) i);
	return j;
}

float get_ieee(void *i) {
	unsigned char *_i, *_j;
	float j;
	if (is_little_endian()) {
		_j = (unsigned char*) &j;
		_i = (unsigned char*) i;
		_j[0] = _i[3];
		_j[1] = _i[2];
		_j[2] = _i[1];
		_j[3] = _i[0];
	} else
		j = *((float*) i);
	return j;
}

void set_short(short value, void *i) {
	if (is_little_endian())
		swab((const char*) &value, (unsigned char*) i, sizeof(short));
	else
		*((short*) i) = value;
}

void set_ushort(unsigned short value, void *i) {
	if (is_little_endian())
		swab((const char*) &value, (unsigned char*) i, sizeof(unsigned short));
	else
		*((unsigned short*) i) = value;
}

void set_int(int value, void *i) {
	unsigned char *_i, *_j;
	if (is_little_endian()) {
		_j = (unsigned char*) i;
		_i = (unsigned char*) &value;
		_j[0] = _i[3];
		_j[1] = _i[2];
		_j[2] = _i[1];
		_j[3] = _i[0];
	} else
		*((int*) i) = value;
}

void set_ieee(float value, void *idx) {
	*((float*) idx) = value;
	swap4(idx);
	//printf("%lf==%lf, ", value, get_ieee(idx));
}

void set_ibm(double value, void *idx) {
	int v, mantissa, tmp;
	float *fint;
	fint = (float*) &v;
	*fint = (float) value;

	if (v != 0) {
		mantissa = (0x007fffff & v) | 0x00800000;
		tmp = ((0x7f800000 & v) >> 23) - 126;
		while (tmp & 0x3) {
			mantissa = mantissa >> 1;
			tmp++;
		}
		v = (0x80000000 & v) | (((tmp >> 2) + 64) << 24) | mantissa;
	}

	unsigned char *bytes = (unsigned char *) idx;

	bytes[0] = v >> 24;
	bytes[1] = ((v & 0xff0000) >> 16);
	bytes[2] = ((v & 0xff00) >> 8);
	bytes[3] = (v & 0xff);

	return;
}

/* Given the trace header and a type, return the string
 * representation of the value stored at the given offset.
 */
char *get_str_val(unsigned char *buf_data, char *value, int offset,
		unsigned char type) {
	switch (type) {
	case 'S': /* Short (2 bytes) */
	{
		short _p;
		_p = *((short *) (buf_data + offset));
		sprintf(value, "%d", get_short(&_p));
		return value;
		break;
	}
	case 'I': /* Integer (4 bytes) */
	{
		int _p;
		_p = *((int*) (buf_data + offset));
		sprintf(value, "%d", get_int(&_p));
		return value;
		break;
	}
	case 'F': /* IEEE (4 bytes) */
	{
		float _p;
		_p = *((float*) (buf_data + offset));
		sprintf(value, "%f", get_ieee(&_p));
		return value;
		break;
	}
	default:
		strcpy(value, "Unknown type");
		return value;
	}
}

/* Given the trace header, a value type and a string
 * representation of that value,
 * stores at the given offset of the trace header that value.
 */
void set_str_val(unsigned char *buf_data, char *value, int offset,
		unsigned char type) {
	int s_zero, s_uno;
	int zero, uno, due, tre;

	if (is_little_endian()) {
		s_zero = 0;
		s_uno = 1;
		zero = 0;
		uno = 1;
		due = 2;
		tre = 3;
	} else {
		s_zero = 1;
		s_uno = 0;
		zero = 3, uno = 2, due = 1, tre = 0;
	}

	switch (type) {
	case 'S': /* Short (2 bytes) */
	{
		short lvalue;
		unsigned char *_p;
		lvalue = atoi(value);
		_p = (unsigned char*) &lvalue;

		buf_data[offset] = _p[s_uno];
		buf_data[offset + 1] = _p[s_zero];
		break;
	}
	case 'I': /* Integer (4 bytes) */
	{
		int lvalue;
		unsigned char *_p;
		lvalue = atol(value);
		_p = (unsigned char*) &lvalue;

		buf_data[offset] = _p[tre];
		buf_data[offset + 1] = _p[due];
		buf_data[offset + 2] = _p[uno];
		buf_data[offset + 3] = _p[zero];

		break;
	}
	case 'F': /* IEEE (4 bytes) */
	{
		float lvalue;
		unsigned char *_p;
		lvalue = atof(value);
		_p = (unsigned char*) &lvalue;

		buf_data[offset] = _p[tre];
		buf_data[offset + 1] = _p[due];
		buf_data[offset + 2] = _p[uno];
		buf_data[offset + 3] = _p[zero];

		break;
	}
	default:
		prerror_and_exit("Unknown type.");
		break;
	}
}

/* SHOW USAGE
 */
void print_usage(int argc, char **argv) {
	fprintf(stdout, "\nsegy-change: a progam to read and optionally change\n"
			"the contents of segy seismic data files.\n\n");
	fprintf(stdout,
			"If you use this software for scientific publication, please cite:\n"
			"\n"
			"G. Stanghellini, G. Carrara\n"
			"Segy-change: The swiss army knife for the SEG-Y files\n"
			"SoftwareX, Volume 6, 2017, Pages 42–47, DOI:10.1016/j.softx.2017.01.003\n"

			"Copyright (C) 2016, Giuseppe Stanghellini (1), Gabriela Carrara (2).\n"
					"  (1) Istituto di Scienze Marine, Geologia Marina, CNR, Bologna, Italy\n"
					"  (2) Independent Researcher, Monte San Pietro, Italy\n\n"
					"This program comes with ABSOLUTELY NO WARRANTY! see GPLV3.\n"
					"This is free software, and you are welcome to redistribute it\n"
					"under certain conditions; for details see the GPLV3 license\n"
					"as published in 'http://www.gnu.org/licenses/'\n\n");

	fprintf(stdout, "Usage :");
	fprintf(stdout, "\tsegy-change [options] -f input_file [-o output_file]\n");
	fprintf(stdout,
			"\n"
					" Input & output:\n"
					"   -f input_file         : The input file (use \"-\" for stdin)\n"
					"   -o output_file        : The output file (use \"-\" for stdout)\n"
					"\n"
					" Switches:\n"
					"   -flip_endianess       : Flip endianess, useful to make\n"
					"                           the file conformant to SEGY standard.\n"
					"   -all                  : Process the whole file. If you do not use this one\n"
					"                           you should use the '-record' and/or '-trace' switches\n"
					"                           to set the part of the segy file that will be\n"
					"                           processed.\n"
					"   -change_header_fields : Change SEGY header fields given an offset,\n"
					"                           Use the following syntax:\n"
					"                           offset0:type0:value0,...\n"
					"                           where type can be S for short (2 bytes) or\n"
					"                           I for int (4 bytes)\n"
					"                           or\n"
					"                           parameter0_name:value,parameter1_name:value,...\n"
					"                           To know at which offsets a field is stored\n"
					"                           or which parameters are available, use the\n"
					"                           '-segy_info and -use_names switches.'\n"
					"   -change_trace_fields  : Change trace header fields given a file with values.\n"
					"                           The file must be in the same format as the output\n"
					"                           obtained by the -dump_trace_fields switch,\n"
					"                           with or without the -use_names switch.\n"
					"   -dump_xy              : Print x and y location to stdout.\n"
					"    SOURCE|RECEIVER        the word SOURCE or RECEIVER must be given, according\n"
					"                           to locations to print.\n"
					"   -vertical_stack num   : Do a vertical stack, summing num samples together\n"
					"   -add_xy               : Add coordinates (in m, feet or arcsec) to trace\n"
					"    fname,SOURCE|RECEIVER  header, given a filename and one of the two words\n"
					"                           SOURCE or RECEIVER, depending on whose entries of the\n"
					"                           trace headers we want to update.\n"
					"                           The specified file must contain for each line:\n"
					"                           the original_field_record_number,\n"
					"                           the trace_sequence_number_within_reel,\n"
					"                           the trace_number_within_field_record,\n"
					"                           the x coordinate,\n"
					"                           the y coordinate,\n"
					"                           the z elevation,\n"
					"                           the unit of measure,\n"
					"                           ie:\n"
					"                           1 1 1 3600 4500 100 meters\n"
					"                           2 1 1 3602 4400 101 arcsec\n"
					"                           2 1 2 3602 4400 103 meters\n"
					"                           5 1 1 3608 4200 105 meters\n"
					"                           ...\n"
					"                           NOTE: to have the lines with the first three numbers\n"
					"                           of a given file, execute:\n"
					"                           segy-change -print_rec_seq_num -all -f segyfile\n"
					"                           you can then just add coordinates and the unit of\n"
					"                           measure in the proper format for each line.\n"
					"   -do_op   +|-|*|/:VAL  : Do the given operation, for example to multiply all\n"
					"                           trace values by 3.0 use: -do_op *:3.0\n"
					"   -convert S|I|F|E      : Convert the trace encoding to the given format, where\n"
					"                           S: short\n"
					"                           I: integer\n"
					"                           F: IBM floating point\n"
					"                           E: IEEE 754.\n"
					"   -do_ps siz,n,sc       : Plot the segy to a postscript file where:\n"
					"                           siz = paper size (for example A4)\n"
//			        "                           custom dimensions can be specified with AC=600x2000\n"
					"                           n = number of traces per cm to plot\n"
					"                           sc = factor to multiply trace values, darkening\n"
					"                                (if>1) or lightening (if<1) the plot.\n"
					"                           Valid paper sizes are: A0,A1,A2,A3,A4.\n"
					"   -print_rec_seq_num    : Print to stdout the original_field_record_number,\n"
					"                           trace_sequence_within_line and\n"
					"                           trace_sequence_within_field_record terna.\n"
					"   -scan                 : Scan SEGY and print some info.\n"
					"   -dump                 : Dump traces values to stdout.\n"
					"   -dump_header_fields   : Dump header fields stored into the traces header.\n"
					"                           Use it with the following syntax:\n"
					"                           field_offset0:field_type0,offset1:type1,...\n"
					"                           where type can be S for short (2 bytes) or \n"
					"                           I for int (4 bytes)\n"
					"                           or\n"
					"                           parameter0_name,parameter1_name,...\n"
					"                           To know at which offsets a field is stored\n"
					"                           or which parameters are available, use the\n"
					"                           '-segy_info and -use_names switches.'\n"
					"   -dump_trace_fields    : Dump fields from trace header given an offset,\n"
					"                           syntax is field_offset0:field_type0,...\n"
					"                           where type can be S for short (2 bytes),\n"
					"                           I for int (4 bytes)\n"
					"                           or\n"
					"                           parameter0_name,parameter1_name,...\n"
					"                           To know at which offsets a field is stored\n"
					"                           or which parameters are available, use the\n"
					"                           '-segy_info and -use_names switches.'\n"
					"   -only_traces_with     : Keep only traces with fields from trace header\n"
					"                           with a given value\n"
					"                           syntax is field_offset0:field_type0:value,...\n"
					"                           where type can be S for short (2 bytes),\n"
					"                           I for int (4 bytes)\n"
					"                           or\n"
					"                           parameter0_name:value,parameter1_name:value,...\n"
					"                           To know at which offsets a field is stored\n"
					"                           or which parameters are available, use the\n"
					"                           '-segy_info and -use_names switches.'\n"
					"   -use_names            : use names instead of offsets.\n"
					"                           to have a list of available parameters names use:\n"
					"                           segy-change -segy_info -use_names\n"
					"   -info                 : Write the SEGY header and exits.\n"
					"   -EBCDIC file          : Replace the EBCDIC header section with the content\n"
					"                           of the specified ASCII file, converting it to EBCDIC\n"
					"   -irc   num            : Do a record renumbering starting from num.\n"
					"   -itc   num            : Do a trace renumbering starting from num.\n"
					"   -no_header            : Do not write the SEGY header into the output file.\n"
					"                           Useful to concatenate SEGY files or to make a file\n"
					"                           in SU format to be processed with Seismic Unix.\n"
					"   -segy_info            : Write the SEGY header formats and exits.\n"
					"   -record num num       : The record interval to process. Tested against\n"
					"                           the 'Original_field_record_number' field\n\n"
					"   -trace num num        : Set the trace interval to process. Tested against\n"
					"                           the 'Trace_number_within_field_record' field.\n"
					"   -only_n_traces num     : Consider the file long at most 'num' traces.\n"
					"   -skip_n_traces num     : Skip the first 'num' traces.\n"
					"   -only_n_samples num    : Consider the trace only long at most 'num' samples.\n"
					"   -skip_n_samples num    : Skip the first 'num' samples.\n"
					"   -num_trace_offset num : Set the offset at which, into the trace header,\n"
					"                           there is the field to use for the trace number.\n"
					"                           Useful when we don't want to use the \n"
					"                           'Trace_number_within_field_record' field.\n"
					"   -traces_per_record num: If the traces per record value stored into the segy\n"
					"                           header is wrong you can use this to override it.\n"
					"                           Please note that in this way you do not change the\n"
					"                           actual value inside the header, to do this you must\n"
					"                           use -change_header_fields switch.\n"
					"   -samples_per_trace    : If the samples per trace value stored into the traces\n"
					"                           header is wrong you can use this to override it.\n"
					"                           Please note that in this way you do not change the\n"
					"                           actual value inside the header, to do this you must\n"
					"                           use -change_trace_fields switch.\n"
					"   -v     num            : Verbosity level: 1 | 2 | 3\n"
					"   -x     num            : Skip num bytes at the beginning of input file.\n"
					"   -view                 : Display segy data.\n"
					"\n\n"
					" Examples:\n\n"
					" - To convert a whole file to a 4-byte ibm floating point format:\n\n"
					"   segy-change -all -f in.segy -o out.segy -convert F\n"
					"\n"
					" - To change the value of 'Number of samples per data trace for this reel'\n"
					"   into the segy header:\n\n"
					"   segy-change -all -f in.segy -o out.segy -change_header_fields 3220:S:3000\n"
					"\n"
					" - To dump some trace fields from all the trace headers:\n\n"
					"   segy-change -all -f in.sgy -dump_trace_fields 0:I,8:I\n"
					"\n"
					" - To change the 'Traces per record field' into the segy header:\n\n"
					"   segy-change -all -f in.sgy -o out.sgy -change_header_fields 3212:S:3336\n"
					"               -traces_per_record 3336\n\n"
					" - To extract traces containing exactly 10000 samples:\n"
					"\n"
					"   segy-change -f in.segy -all -o out.segy -use_names\n"
					"               -only_traces_with NUMBER_OF_SAMPLES_IN_THIS_TRACE:10000\n"
					" - To make a postscript plot of the first two shot-gather, with 25 traces\n"
			        "   per cm:\n"
					"\n"
					"   segy-change -f in.segy -record 100 101 -trace 1 120\n"
					"               -do_ps A4,25,0.01\n"
					"\n"
					" - To view interactively the segy:\n"
					"\n"
					"   segy-change -f in.segy -view\n"
					"\n"
			);

	printf("For your information this machine is ");
	if (is_little_endian())
		printf("LITTLE_ENDIAN ");
	else
		printf("BIG_ENDIAN ");

	printf("and the size of int appears\nto be %d bytes", (int) sizeof(int));
	if (sizeof(int) == 4)
		printf(", as it should be.\n\n");
	else
		printf(" and that's wrong! it should be 4 bytes.\n\n");

	printf( "***********************************************************************\n"
			"If you use this software for scientific publication, please cite:\n"
			"G. Stanghellini, G. Carrara\n"
			"Segy-change: The swiss army knife for the SEG-Y files\n"
			"SoftwareX, Volume 6, 2017, Pages 42–47, DOI:10.1016/j.softx.2017.01.003\n"
			"***********************************************************************\n"
			);


	exit(0);
}

void get_source_xy(SEGY_file *file, double *x, double *y,
		int *unit_of_measure_1_feet_or_meters_2_arcsec) {
	double mult_scalar;
	SEGY_trace_header *header = &file->trace_header;
	mult_scalar = GET_SEGYTRACEH_Coord_mult_scalar(header);
	if (mult_scalar == 0) mult_scalar = 1;
	*unit_of_measure_1_feet_or_meters_2_arcsec =
			GET_SEGYTRACEH_Coord_units_1_fm_2_arcsec(header);
	if (mult_scalar < 0)
		mult_scalar = -1 / mult_scalar;
	*x = GET_SEGYTRACEH_Source_x_or_lon(header) * mult_scalar;
	*y = GET_SEGYTRACEH_Source_y_or_lat(header) * mult_scalar;
}

void get_receiver_xy(SEGY_file *file, double *x, double *y,
		int *unit_of_measure_1_feet_or_meters_2_arcsec) {
	double mult_scalar;
	SEGY_trace_header *header = &file->trace_header;
	mult_scalar = GET_SEGYTRACEH_Coord_mult_scalar(header);
	if (mult_scalar == 0) mult_scalar = 1;
	*unit_of_measure_1_feet_or_meters_2_arcsec =
			GET_SEGYTRACEH_Coord_units_1_fm_2_arcsec(header);
	if (mult_scalar < 0)
		mult_scalar = -1.0 / mult_scalar;
	*x = GET_SEGYTRACEH_Receiver_x_or_lon(header) * mult_scalar;
	*y = GET_SEGYTRACEH_Receiver_y_or_lat(header) * mult_scalar;
}

int coordinates_scaling_factor;

void read_xy(char *fname) {
	num_coords = 0;
	max_x = max_y = -1e100;
	min_x = min_y = 1e100;

	FILE *fp = fopen(fname, "r");
	if (fp == NULL)
		prerror_and_exit("Cannot open '%s' file, aborting.\n", fname);

	int _idum;
	double _fdum;
	char tmp_c[1024000];

	/* Find out how many records we have.
	 */
	while (!feof(fp)) {
		fscanf(fp, "%d %d %d %lf %lf %s", &_idum, &_idum, &_idum, &_fdum,
				&_fdum, tmp_c);
		num_coords++;
	}

	// Rewind the file and allocate the needed memory
	rewind(fp);
	coords = (COORDS*) calloc(sizeof(COORDS), num_coords + 1);

	/* Read the actual data, and find the max and min in order to properly scale
	 * the numbers.
	 */
	num_coords = 0;
	while (!feof(fp)) {
		fscanf(fp, "%d %d %d %lf %lf %s",
				&(coords[num_coords].original_field_record),
				&(coords[num_coords].trace_seq_within_reel),
				&(coords[num_coords].trace_seq_within_field_record),
				&(coords[num_coords].lon_or_x), &(coords[num_coords].lat_or_y),
				tmp_c);

		if (strcmp(tmp_c, "feet") == 0 || strcmp(tmp_c, "meters") == 0)
			coords[num_coords].unit_of_measure_1_feet_or_meters_2_arcsec = 1;
		if (strcmp(tmp_c, "arcsec") == 0)
			coords[num_coords].unit_of_measure_1_feet_or_meters_2_arcsec = 2;

		max_x = my_max(max_x, coords[num_coords].lon_or_x);
		max_y = my_max(max_y, coords[num_coords].lat_or_y);
		min_x = my_min(min_x, coords[num_coords].lon_or_x);
		min_y = my_min(min_y, coords[num_coords].lat_or_y);
		num_coords++;
	}

	double _max = my_max(fabs(my_max(max_x, max_y)),
			fabs(my_min(min_x, min_y)));
	double scaling = 2147483647.0 / _max;
	if (scaling > 0.0001)
		coordinates_scaling_factor = 10000;
	if (scaling > 0.001)
		coordinates_scaling_factor = 1000;
	if (scaling > 0.01)
		coordinates_scaling_factor = 100;
	if (scaling > 0.1)
		coordinates_scaling_factor = 10;
	if (scaling > 1)
		coordinates_scaling_factor = 1;
	if (scaling > 10)
		coordinates_scaling_factor = -10;
	if (scaling > 100)
		coordinates_scaling_factor = -100;
	if (scaling > 1000)
		coordinates_scaling_factor = -1000;
	if (scaling > 10000)
		coordinates_scaling_factor = -10000;

	/* Compute optimal scaling factor for coordinates.
	 */

	if (verbose == 1) {
		printf("%d entries read from '%s' file\n", num_coords, fname);
		printf("%lf max magnitude found.\n", _max);
		printf("scaling factor = %d.\n", coordinates_scaling_factor);
	}

	/* Scale the values.
	 */
	if (coordinates_scaling_factor < 0)
		scaling = -1.0 / coordinates_scaling_factor;
	else
		scaling = coordinates_scaling_factor;
	int i;
	for (i = 0; i < num_coords; i++) {
		coords[i].lon_or_x = coords[i].lon_or_x / scaling; // divide because when read the scaling is a multiply factor
		coords[i].lat_or_y = coords[i].lat_or_y / scaling; // divide because when read the scaling is a multiply factor
	}

	/* all done, we are ready to use the data.
	 */
}

/* SHOW SEGY HEADER INFORMATION.
 */
void print_segy_info() {
	if (use_names) {
		int j = 0;
		fprintf(stdout,
				"--------------------------------------------------------------\n"
						"                        SEG Y FORMAT\n"
						"--------------------------------------------------------------\n"
						"SEG Y format was developed by the Society of\n"
						"Exploration Geophysicists (SEG) as an exchange format for\n"
						"demultiplexed seismic data on 9-track tape.\n"
						"A tape reel is divided into a tape identification header\n"
						"followed by multiple trace data blocks. The tape identification\n"
						"header can also be omitted for data residing on disks.\n"
						"\n\n"
						"--------------------------------------------------------------\n"
						"             SEG Y tape Identification Header\n"
						"--------------------------------------------------------------\n"
						"The SEGY tape identification header is made by an EBCDIC\n"
						"block and by a binary block. The EBCDIC block is\n"
						"composed of 40 cards (80 bytes per card) for a total of\n"
						"3200 bytes. Each card begins with 'C' and bytes at offsets\n"
						"ranging from 29 to 39 are left unassigned for optional use.\n"
						"Unused cards are filled with EBCDIC blank code.\n"
						"After the EBCDIC block begins a 400 byte long binary block.\n\n"
						"Offset Type           Parameter\n");
		while (segy_header_names[j]) {
			fprintf(stdout, "%5d ", segy_header_types[j * 2]);
			switch (segy_header_types[j * 2 + 1]) {
			case 'I':
				fprintf(stdout, " Integer ");
				break;
			case 'S':
				fprintf(stdout, " Short   ");
				break;
			}
			fprintf(stdout, " %s\n", segy_header_names[j]);
			j++;
		}
		fprintf(stdout,
				"\n--------------------------------------------------------------\n"
						"                      SEG Y Trace Header\n"
						"--------------------------------------------------------------\n\n"
						"Offset Type           Parameter\n");
		j = 0;
		while (trace_header_names[j]) {
			fprintf(stdout, "%5d ", trace_header_types[j * 2]);
			switch (trace_header_types[j * 2 + 1]) {
			case 'I':
				fprintf(stdout, " Integer ");
				break;
			case 'S':
				fprintf(stdout, " Short   ");
				break;
			}
			fprintf(stdout, " %s\n", trace_header_names[j]);
			j++;
		}

		return;
	}
	fprintf(stdout,
			"--------------------------------------------------------------\n"
					"                        SEG Y FORMAT\n"
					"--------------------------------------------------------------\n"
					"SEG Y format was developed by the Society of\n"
					"Exploration Geophysicists (SEG) as an exchange format for\n"
					"demultiplexed seismic data on 9-track tape.\n"
					"A tape reel is divided into a tape identification header\n"
					"followed by multiple trace data blocks. The tape identification\n"
					"header can also be omitted for data residing on disks.\n"
					"\n\n"
					"--------------------------------------------------------------\n"
					"             SEG Y tape Identification Header\n"
					"--------------------------------------------------------------\n"
					"The SEGY tape identification header is made by an EBCDIC\n"
					"block and by a binary block. The EBCDIC block is\n"
					"composed of 40 cards (80 bytes per card) for a total of\n"
					"3200 bytes. Each card begins with 'C' and bytes at offsets\n"
					"ranging from 29 to 39 are left unassigned for optional use.\n"
					"Unused cards are filled with EBCDIC blank code.\n"
					"After the EBCDIC block begins a 400 byte long binary block.\n"
					"\n"
					"Offset Type     Description\n"
					"0      EBCDIC   EBCDIC block.\n"
					"3200   Integer  Job Identification Number\n"
					"3204   Integer  Line number (1 line per reel)\n"
					"3208   Integer  Reel number\n"
					"3212   Short    Number of data traces per record\n"
					"3214   Short    Number of auxiliary traces per record\n"
					"3216   Short    Sample interval for this reel (microseconds)\n"
					"3218   Short    Sample interval for original field recording\n"
					"                (microseconds)\n"
					"3220   Short    Number of samples per data trace for this\n"
					"                reel\n"
					"3222   Short    Number of samples per data trace (original\n"
					"                field recording)\n"
					"3224   Short    Data sample format code: 1 = 32-bit float;\n"
					"                2 = 32-bit integer\n"
					"                3 = 16-bit integer; 5 = ieee 754 single precision\n"
					"3226   Short    Nominal CDP fold\n"
					"3228   Short    Trace sorting code: \n"
					"                1 = as recorded (no sorting)\n"
					"                2 = CDP ensemble; \n"
					"                3 = single fold continuous profile\n"
					"                4 = horizontally stacked\n"
					"3230   Short    Number of vertically summed traces\n"
					"3232   Short    Sweep frequency at start (Hz)\n"
					"3234   Short    Sweep frequency at end (Hz)\n"
					"3236   Short    Sweep length (ms)\n"
					"3238   Short    Sweep type: \n"
					"                1 = linear; 2 = parabolic\n"
					"                3 = exponential; 4 = other\n"
					"3240   Short    Trace number of sweep channel\n"
					"3242   Short    Sweep taper length at start (ms)\n"
					"3244   Short    Sweep taper length at end (ms)\n"
					"3246   Short    Taper type: \n"
					"                1 = linear 2 = cosine squared; 3 = other\n"
					"3248   Short    Correlated data traces: 1 = yes, 2 = no\n"
					"3250   Short    Binary gain recovered: 1 = yes, 2 = no\n"
					"3252   Short    Amplitude recovery method: 1 = none\n"
					"                2 = spherical divergence; 3 = AGC; 4 = none\n"
					"3254   Short    Measurement system: 1 = feet, 2 = meters\n"
					"3256   Short    Impulse signal: 1 = negative amplitude means\n"
					"                increased pressure or upward\n"
					"                movement; 2 = positive amplitude means\n"
					"                increased pressure or upward movement\n"
					"3258   Short    Vibratory polarity code (seismic signal lags\n"
					"                pilot trace by):\n"
					"                1 = 337.5 to 22.5 degrees, \n"
					"                2 = 22.5 to 67.5 degrees\n"
					"                3 = 67.5 to 112.5 degrees, \n"
					"                4 = 112.5 to 157.5 degrees\n"
					"                5 = 157.5 to 202.5 degrees, \n"
					"                6 = 202.5 to 247.5 degrees\n"
					"                7 = 247.5 to 292.5 degrees, \n"
					"                8 = 292.5 to 337.5 degrees\n"
					"3260...3599      Unassigned\n\n"
					"--------------------------------------------------------------\n"
					"                      SEG Y Trace Header\n"
					"--------------------------------------------------------------\n"
					"Offset Type     Standard SEG Y\n"
					"0      Integer  Trace sequence number within line\n"
					"4      Integer  Trace sequence number within reeln\n"
					"8      Integer  Original field record number\n"
					"12     Integer  Trace number within field record\n"
					"16     Integer  Source point number\n"
					"20     Integer  CDP number\n"
					"24     Integer  CDP sequece number\n"
					"28     Short    Trace identification code:\n"
					"                1 = seismic   4 = time break  7 = timing\n"
					"                2 = dead      5 = uphole      8 = water break\n"
					"                3 = dummy     6 = sweep       9 = option use\n"
					"                                                             \n"
					"30     Short    Number of vertically summed traces\n"
					"32     Short    Number of horizontally summed traces (fold)\n"
					"34     Short    Data use:\n"
					"                1 = production; 2 = test\n"
					"36     Integer  Source-receiver offset in feet or meters\n"
					"40     Integer  Receiver group elevation in feet or meters;\n"
					"                positive above sea level\n"
					"44     Integer  Surface elevation at source (feet or meters)\n"
					"48     Integer  Source depth below surface (>0)\n"
					"52     Integer  Datum elevation at receiver group\n"
					"56     Integer  Datum elevation at source\n"
					"60     Integer  Water depth at source\n"
					"64     Integer  Water depth at receiver group\n"
					"68     Short    Elevation multiplication scalar for bytes at\n"
					"                offsets 40-67 = 1, 10, 100, 1000 or 10,000\n"
					"                positive indicates multiplication; negative, division\n"
					"70     Short    Coordinate multiplication scalar for bytes at offsets\n"
					"                72-87 (see bytes 68-69)\n"
					"72     Integer  Source X (feet or meters) or\n"
					"                longitude (seconds of arc; + is E, - is W)\n"
					"76     Integer  Source Y (feet or meters) or\n"
					"                latitude (seconds of arc; + is N, - is S)\n"
					"80     Integer  Receiver X (feet or meters) or\n"
					"                longitude (seconds of arc; + is E, - is W)\n"
					"84     Integer  Receiver Y (feet or meters) or\n"
					"                latitude (seconds of arc; + is N, - is S)\n"
					"88     Short    Coordinate units\n"
					"                1 = feet or meters; 2 = seconds of arc\n"
					"90     Short    Weathering velocity\n"
					"92     Short    Sub-weathering velocity\n"
					"94     Short    Uphole time at source (ms)\n"
					"96     Short    Uphole time at group (ms)\n"
					"98     Short    Source static correction (ms)\n"
					"100    Short    Receiver static correction (ms)\n"
					"102    Short    Total static applied\n"
					"104    Short    Lag time A between trace header time and time\n"
					"                break (ms)\n"
					"106    Short    Lag time B between time break and source time\n"
					"                (ms)\n"
					"108    Short    Delay time between source and recording time\n"
					"                (ms)\n"
					"110    Short    Brute start time (ms)\n"
					"112    Short    Mute end time (ms)\n"
					"114    Short    Number of samples in this trace\n"
					"116    Short    Sample interval (microseconds)\n"
					"118    Short    Gain type: 1 = fixed; 2 = binary\n"
					"                3 = floating point; 4 = optional use\n"
					"120    Short    Instrument Gain Constant \n"
					"122    Short    Instrument early or initial gain\n"
					"124    Short    Correlated; 1 = yes, 2 = no\n"
					"126    Short    Sweep frequency at start (Hz)\n"
					"128    Short    Sweep frequency at end (Hz)\n"
					"130    Short    Sweep length (ms)\n"
					"132    Short    Sweep type: 1 = linear; 2 = parabolic\n"
					"                3 = exponential; 4 = other\n"
					"134    Short    Sweep taper length at start (ms)\n"
					"136    Short    Sweep taper length at end (ms)\n"
					"138    Short    Sweep taper type: 1 = linear          \n"
					"                2 = cosine squared; 3 = other\n"
					"140    Short    Alias filter frequency (Hz), if used  \n"
					"142    Short    Alias filter slope (dB/octave)        \n"
					"144    Short    Notch filter frequency (Hz), if used  \n"
					"146    Short    Notch filter slope (dB/octave)    \n"
					"148    Short    Low-cut frequency (Hz), if used   \n"
					"150    Short    High-cut frequency (Hz), if used  \n"
					"152    Short    Low-cut filter slope (dB/octave)  \n"
					"154    Short    High-cut filter slope (dB/octave) \n"
					"156    Short    Year data recorded                \n"
					"158    Short    Day of year                       \n"
					"160    Short    Hour of day (24 hour clock)       \n"
					"162    Short    Minute of hour                    \n"
					"164    Short    Second of minute (for trace start)\n"
					"166    Short    Time basis code: 1=local 2=GMT 3=other\n"
					"168    Short    Trace weighting factor: 2*N volts for least\n"
					"                significant bit where N = 0,1,...,32.767\n"
					"170    Short    Receiver group number at roll switch position 1      \n"
					"172    Short    Receiver group number for first trace in field record\n"
					"174    Short    Receiver group number for last trace in field record \n"
					"176    Short    Gap size (number of receiver groups dropped)         \n"
					"178    Short    Overtravel associated with taper at start or end     \n"
					"                of line: 1 = down (or behind), 2 = up (or ahead)\n"
					"180-236         Unassigned\n"
					"\n"
					"                              \n");
}

int trace_sample_length(SEGY_file *segy_file) {
	int nb;
	if (GET_SEGYH_Data_sample_format_code(&segy_file->header) == 5)
		nb = 4;
	if (GET_SEGYH_Data_sample_format_code(&segy_file->header) == 2)
		nb = 4;
	if (GET_SEGYH_Data_sample_format_code(&segy_file->header) == 3)
		nb = 2;
	if (GET_SEGYH_Data_sample_format_code(&segy_file->header) == 1)
		nb = 4;
	return nb;
}

int trace_data_length(SEGY_file *segy_file) {
	return trace_sample_length(segy_file)
			* GET_SEGYTRACEH_USHORT_Number_of_samples_in_this_trace(
					&segy_file->trace_header);
}

void ibm2ieee(void *to, const void *from, int len) {
	register unsigned fr;
	register int exp;
	register int sgn;

		fr = ntohl(*(long *) from);
		sgn = fr >> 31;
		fr <<= 1;
		exp = fr >> 25;
		fr <<= 7;

		if (fr == 0) {
			exp = 0;
			goto done;
		}

		exp = (exp << 2) - 130;

		while (fr < 0x80000000) {
			--exp;
			fr <<= 1;
		}

		if (exp <= 0) {
			if (exp < -24)
				fr = 0;
			else
				fr >>= -exp;
			exp = 0;
		} else if (exp >= 255) {
			fr = 0;
			exp = 255;
		} else {
			fr <<= 1;
		}

		done:
		*(unsigned *) to = (fr >> 9) | (exp << 23) | (sgn << 31);
}

/* Convert from segy_file->format to double
 */
double get_val(SEGY_file *segy_file, void *addr, int idx) {
	float value;

	if (GET_SEGYH_Data_sample_format_code(&segy_file->header) == 5)
		return (double) get_ieee((float*) addr + idx); /* IEEE 754 */
	if (GET_SEGYH_Data_sample_format_code(&segy_file->header) == 3)
		return (double) get_short((short*) addr + idx); /* SHORT */
	if (GET_SEGYH_Data_sample_format_code(&segy_file->header) == 2)
		return (double) get_int((int*) addr + idx); /* LONG */

	if (GET_SEGYH_Data_sample_format_code(&segy_file->header) == 1) {
		ibm2ieee(&value, addr + idx * 4, 1);
		return value;
//		chr = chr + idx * 4;
//		power = ((*chr) & 0x7f) - 64;
//		factor = pow(16.0, (double) power);
//		mantissa = 0.0;
//		a16 = 16.0;
//		for (i = 1; i < 4; i++) {
//			unsigned char tmp_uc;
//			tmp_uc = chr[i];
//			mantissa += (double) ((tmp_uc & 0xf0) >> 4) / a16;
//			a16 *= 16.0;
//			mantissa += (double) (tmp_uc & 0x0f) / a16;
//			a16 *= 16.0;
//		}
//		if (chr[0] & 0x80)
//			return -mantissa * factor;
//		else
//			return mantissa * factor;
	} else {
		prerror_and_exit("Error: Unknown Format.\n");
		return 0;
	}
}

/* Convert from double to segy_file->format and store the conversion into *addr + idx.
 */
void set_val(double value_to_convert, SEGY_file *segy_file, void *addr, int idx) {

	switch (GET_SEGYH_Data_sample_format_code(&segy_file->header)) {
	case 1:
		set_ibm((double) value_to_convert, (float*) addr + idx);
		return;
	case 2:
		set_int((int) value_to_convert, (int*) addr + idx); /* INT */
		return;
	case 3:
		set_short((short) value_to_convert, (short*) addr + idx); /* SHORT */
		return;
	case 5:
		set_ieee((float) value_to_convert, (float*) addr + idx); /* IEEE 754 */
		return;
	}

	prerror_and_exit("Error: Unknown Format.\n");
}

void flip_trace_data_endianess(SEGY_file *segy_file) {
	if (flip_endianess) {
		switch (GET_SEGYH_Data_sample_format_code(&segy_file->header)) {
		case 1:
		case 2:
		case 5:
			for (i = 0; i < n_samples; i++)
				swap4((unsigned int*) (segy_file->trace_data + i * 4));
			break;
		case 3:
			for (i = 0; i < n_samples; i++)
				swap2((unsigned short*) (segy_file->trace_data + i * 2));
			break;
		}
	}
}

void flip_trace_header_endianess(SEGY_file *segy_file) {
	/* Flip the endianess if it's the case.
	 */
	if (flip_endianess) {
		int kk = 0;
		while (trace_header_types[kk * 2] != -1) {
			switch (trace_header_types[kk * 2 + 1]) {
			case 'I':
			case 'F':
				swap4(
						(unsigned int *) ((unsigned char*) (&(segy_file->trace_header.HEADER))
								+ trace_header_types[kk * 2]));
				break;
			case 'S':
				swap2(
						(unsigned short *) ((unsigned char *) (&(segy_file->trace_header.HEADER))
								+ trace_header_types[kk * 2]));
				break;
			}
			kk++;
		}
	}
}

void flip_header_endianess(SEGY_file *segy_file) {
	/* Flip the endianess if it's the case.
	 */
	if (flip_endianess) {
		int kk = 0;
		while (segy_header_types[kk * 2] != -1) {
			switch (segy_header_types[kk * 2 + 1]) {
			case 'I':
			case 'F':
				swap4(
						(unsigned int *) ((unsigned char *) (&segy_file->header)
								+ segy_header_types[kk * 2]));
				break;
			case 'S':
				swap2(
						(unsigned short *) ((unsigned char *) (&segy_file->header)
								+ segy_header_types[kk * 2]));
				break;
			}
			kk++;
		}
	}
}

void copy_segy_header(SEGY_file *source, SEGY_file *dest) {
	memcpy(&dest->header, &source->header, sizeof(source->header));
//	if (only_nsamples != -1) {
//		set_short(only_nsamples, dest->header.EBCDIC + 3220);
//	}
}

void copy_segy_trace_header(SEGY_file *source, SEGY_file *dest) {
	memcpy(&dest->trace_header, &source->trace_header,
			sizeof(source->trace_header));
//	if (only_nsamples != -1)
//		set_short(only_nsamples, dest->trace_header.HEADER + 114);
}

void copy_segy_trace_data(SEGY_file *source, SEGY_file *dest) {
	int i;
	if ((dest->trace_data = (unsigned char *) realloc(dest->trace_data,
			trace_data_length(dest)))) {
		int start = 0, end = n_samples;
//		if (skip_nsamples != -1)
//			start = skip_nsamples;
//		if (only_nsamples != -1)
//			end = start + only_nsamples;
		int source_nsamples =
				GET_SEGYTRACEH_USHORT_Number_of_samples_in_this_trace(
						&source->trace_header);
		if (vertical_stack > 1) { // DO THE VERTICAL STACKING
			int j, k = 0;
			for (i = start; i < end; i += vertical_stack) {
				double val = 0;
				for (j = 0; j < vertical_stack; j++)
					val += GET_IN_SAMPLE(i + j);

				SET_OUT_SAMPLE(val, k);
				k++;
			}

		} else {
			memset(dest->trace_data, 0, trace_data_length(dest));
			if (GET_SEGYH_Data_sample_format_code(
					&source->header) != GET_SEGYH_Data_sample_format_code(&dest->header)) {
				dest->trace_data = (unsigned char *) realloc(dest->trace_data,
						trace_data_length(dest));
				for (i = start; i < my_min(end, source_nsamples); i++) {
					set_val(get_val(source, source->trace_data, i), dest,
							dest->trace_data, i - start);
				}
			} else {
//				if (skip_nsamples != -1)
//					memcpy(dest->trace_data,
//							source->trace_data
//									+ skip_nsamples
//											* trace_sample_length(source),
//							trace_data_length(dest));
//				else
					memcpy(dest->trace_data, source->trace_data,
							trace_data_length(dest));
			}
		}
	} else {
		prerror_and_exit("Error: Cannot allocate memory for TRACE DATA.\n");
	}
}

void write_segy_header(SEGY_file *dest) {
	if (replace_ebcdic) {
		memcpy(&dest->header, my_ebcdic, 3200);
	}

	fwrite(&dest->header, sizeof(dest->header), 1, dest->fp);
}

void write_segy_trace_header(SEGY_file *dest) {
	fwrite(&dest->trace_header, sizeof(dest->trace_header), 1, dest->fp);
}

void write_segy_trace_data(SEGY_file *dest) {
	fwrite(dest->trace_data, 1, trace_data_length(dest), dest->fp);
}

/* OPEN THE SEGY FILE.
 * RETURN: 0 ON SUCCESS, -1 OTHERWISE.
 */
int open_segy(SEGY_file *segy_file, char *fname, char *mode, off_t initial_seek) {
	if (fname[0] == '-' && fname[1] == 0 && mode[0] == 'r') {
		segy_file->fp = stdin;
		segy_file->fname = (char *) malloc(6);
		strncpy(segy_file->fname, "stdin", 5);
		segy_file->fname[5] = 0;
	} else {
		if (fname[0] == '-' && fname[1] == 0  && mode[0] == 'w') {
			segy_file->fp = stdout;
			segy_file->fname = (char *) malloc(7);
			strncpy(segy_file->fname, "stdout", 6);
			segy_file->fname[6] = 0;
		} else {
			segy_file->fp = (FILE*)fopen(fname, mode);
			segy_file->fname = (char *) malloc(strlen(fname) + 1);
			strncpy(segy_file->fname, fname, strlen(fname));
			segy_file->fname[strlen(fname)] = 0;
		}
	}

	if (segy_file->fp != NULL) {
		fseek(segy_file->fp, initial_seek, SEEK_SET);
		return 0;
	} else
		return -1;
}

void do_change_header() {
	if (change_header_fields) {
		int field_segy_offset;
		char field_segy_type;
		int m_field_nr = 1;
		char m_field[1000], m_fields[1000], value[1000];
		while (get_field(header_fields_to_change, m_field_nr, m_fields, ',')) {
			get_field(m_fields, 1, m_field, ':');
			field_segy_offset = atoi(m_field);
			get_field(m_fields, 2, m_field, ':');
			field_segy_type = m_field[0];
			get_field(m_fields, 3, value, ':');
			if (output_segy)
				set_str_val((unsigned char*) (&out_segy_file.header), value,
						field_segy_offset, field_segy_type);
			set_str_val((unsigned char*) (&segy_file.header), value,
					field_segy_offset, field_segy_type);

			m_field_nr++;
		}
	}

	/* At 3224-3225 is stored the data_sample_format_code
	 * see below in print_segy_info */
	switch (convert_to) {
	case 'S': /* Short int */
		set_str_val((unsigned char*) &out_segy_file.header, "3", 3224, 'S');
		break;
	case 'I': /* int */
		set_str_val((unsigned char*) &out_segy_file.header, "2", 3224, 'S');
		break;
	case 'F': /* ibm floating point */
		set_str_val((unsigned char*) &out_segy_file.header, "1", 3224, 'S');
		break;
	case 'E': /* IEEE 754 floating point */
		set_str_val((unsigned char*) &out_segy_file.header, "5", 3224, 'S');
		break;
	default:
		/* do nothing, no conversion. */
		break;
	}

	if (vertical_stack > 1 && output_segy) {
		set_short(
				get_short((void*) (&out_segy_file.header) + 3216)
						* vertical_stack, (void*) &out_segy_file.header + 3216);
		set_short(
				get_short((void*) (&out_segy_file.header) + 3220)
						/ vertical_stack, (void*) &out_segy_file.header + 3220);
	}
}

/* READ FROM file THE SEGY HEADER
 * Dump header fields if dump_header_fields static variable is true.
 * RETURN: 0 ON SUCCESS, -1 OTHERWISE.
 */
int get_segy_header(SEGY_file *segy_file, int verbose) {
	int bytes_read, i;

	/* READ 3200 + 400 BYTES
	 */
	bytes_read = fread(&(segy_file->header), sizeof(segy_file->header), 1,
			segy_file->fp);
	if (bytes_read != 1)
		prerror_and_exit("Cannot read the SEGY HEADER.\n");

	flip_header_endianess(segy_file);

	if (n_traces == 0)
		n_traces = GET_SEGYH_Number_of_data_traces_per_record(
				&segy_file->header);

	/* Check and dump segy_header fields.
	 */
	if (dump_header_fields) {
		int field_segy_offset, jjj;
		char field_segy_type;
		int m_field_nr = 1;
		char m_field[1000], m_fields[1000], value[1000];

		printf("Segy header fields:\n");
		while (get_field(header_fields_to_dump, m_field_nr, m_fields, ',')) {
			if (m_field_nr > 1)
				printf("; ");
			get_field(m_fields, 1, m_field, ':');
			if (use_names)
				jjj = get_parameter_index_by_offset(segy_header_types,
						atoi(m_field));
			field_segy_offset = atoi(m_field);
			get_field(m_fields, 2, m_field, ':');
			field_segy_type = m_field[0];
			if (use_names)
				printf("%s,%s", segy_header_names[jjj],
						get_str_val((unsigned char*) (&segy_file->header),
								value, field_segy_offset, field_segy_type));
			else
				printf("%d,%c,%s", field_segy_offset, field_segy_type,
						get_str_val((unsigned char*) (&segy_file->header),
								value, field_segy_offset, field_segy_type));
			m_field_nr++;
		}

		printf("\n");
	}

	if (verbose >= 1 || scan) {
		fprintf(stderr,
				"================================== EBCDIC dump ================================\n");
		for (i = 0; i < 3200; i++) {
			if (i % 80 == 0)
				fprintf(stderr, "\n");
			fprintf(stderr, "%c", ebcdic2ascii[segy_file->header.EBCDIC[i]]);
		}
		fprintf(stderr,
				"\n================================= SEGY HEADER =================================\n");
		fprintf(stderr, "                                        Line nr. : %d\n",
				GET_SEGYH_Line_number(&segy_file->header));
		fprintf(stderr, "                                        Reel nr. : %d\n",
				GET_SEGYH_Reel_number(&segy_file->header));
		fprintf(stderr, "                Number of data traces per record : %d\n",
				GET_SEGYH_Number_of_data_traces_per_record(&segy_file->header));
		fprintf(stderr, "    Sample interval for this reel (microseconds) : %d\n",
				GET_SEGYH_Sample_interval_for_this_reel(&segy_file->header));
		fprintf(stderr, "  Number of samples per data trace for this reel : %d\n",
				GET_SEGYH_Number_of_samples_per_datatrace_for_this_reel(
						&segy_file->header));
		int dsfc = GET_SEGYH_Data_sample_format_code(&segy_file->header);
		fprintf(stderr,"                         Data sample format code : %d ", dsfc);
		switch(dsfc)
		{
		case 1:
			fprintf(stderr, "(IBM floating point)\n");
			break;
		case 2:
			fprintf(stderr, "(32-bit integer)\n");
			break;
		case 3:
			fprintf(stderr, "(16-bit integer)\n");
			break;
		case 5:
			fprintf(stderr, "(IEEE 754 single precision)\n");
			break;
		default:
			fprintf(stderr, "(unknown)\n");
			break;
		}
		int uom = GET_SEGYH_Measurements_system(&segy_file->header);
		if(uom == 1) printf("        Measurement system: 1 = feet, 2 = meters : %d (feet)\n", uom);
		else if(uom == 2) printf("        Measurement system: 1 = feet, 2 = meters : %d (meters)\n", uom);
		else printf("        Measurement system: 1 = feet, 2 = meters : %d (unknown)\n", uom);
		fprintf(stderr, "===============================================================================\n");
	}

	if (!no_EBCDIC_stamp) {
		char *s =
				"This segy was processed with segy-change by Giuseppe Stanghellini @ Ismar-CNR";
		char stamp[100];
		strcpy(stamp, s);

		int ii = 0;
		while (stamp[ii]) {
			stamp[ii] = ascii2ebcdic[(int)stamp[ii]];
			ii++;
		}

		memcpy(((unsigned char*) &(segy_file->header)) + 3120, stamp,
				strlen(stamp));
	}

	switch (GET_SEGYH_Data_sample_format_code(&segy_file->header)) {
	case 1:
	case 2:
	case 3:
	case 5:
		break;
	default:
		fprintf(stderr,
				"\nWARNING ! serious error detected on the segy header,\nthe Data_sample_format_code is unknown, please correct it with:\n\n    -change_header_fields 3224:S:value\n\nbefore continuing.\n\n");
		break;
	}

	if (n_traces == 0)
	{
		fprintf(stderr,
				"\n\nWARNING ! the number of data traces per record in segy header is set to 0\nplease fix it with:\n\n  -traces_per_record ntraces  -change_header_fields 3212:S:ntraces\n\nbefore continuing.\nFrom now on, 1 trace per record will be used.\n");
		n_traces = 1;
	}
	return 0;
}

/* check if the actual trace is inside the given ranges and fields_values.
 */
bool keep_trace(SEGY_file *segy_file, int rec_start, int rec_end,
		int trace_start, int trace_end, int trace_offset) {
	int field_segy_offset;
	char field_segy_type;
	int m_field_nr = 1;
	char m_field[1000], m_fields[1000], value[1000], expected_value[1000];

	int tr_num = GET_SEGYTRACEH_Field(&segy_file->trace_header, trace_offset);
	//tr_num =  GET_SEGYTRACEH_Trace_number_within_field_record(&segy_file->trace_header);
	int sh_num = GET_SEGYTRACEH_Original_field_record_number(
			&segy_file->trace_header);
	if (tr_num >= trace_start && tr_num <= trace_end && sh_num >= rec_start
			&& sh_num <= rec_end) {
		if(only_traces_with)
		{
		while (get_field(traces_fields_valid_values, m_field_nr, m_fields, ',')) {
			get_field(m_fields, 1, m_field, ':');
			field_segy_offset = atoi(m_field);
			get_field(m_fields, 2, m_field, ':');
			field_segy_type = m_field[0];
			get_field(m_fields, 3, expected_value, ':');
			if (atof(get_str_val(segy_file->trace_header.HEADER, value,
							field_segy_offset, field_segy_type))
					!= atof(expected_value))
			{
				if(verbose >= 2) fprintf(stderr, "Discarding trace #%d of record #%d because it contains %f and not %f.\n",
   					tr_num, sh_num,
						atof(get_str_val(segy_file->trace_header.HEADER, value,	field_segy_offset, field_segy_type))
					    , atof(expected_value));
				return false;
			}
			m_field_nr++;
		}
		}
		return true;
	}
	else
		return false;
}

void do_skip_ntraces() {
	if (skip_ntraces == 0)
		return;
	long nn = GET_SEGYH_Number_of_samples_per_datatrace_for_this_reel(
			&segy_file.header);

	if (GET_SEGYH_Data_sample_format_code(&segy_file.header) == 2)
		nn *= 4;
	else if (GET_SEGYH_Data_sample_format_code(&segy_file.header) == 3)
		nn *= 2;
	else if (GET_SEGYH_Data_sample_format_code(&segy_file.header) == 1)
		nn *= 4;
	else
		prerror_and_exit("Error: Unknown data sample format code\n");
	fseek(segy_file.fp, (240 + nn) * skip_ntraces, SEEK_CUR);
}

int get_segy_trace(SEGY_file *segy_file, int verbose) {
	int i;

	int bytes_read;
	bytes_read = fread(&segy_file->trace_header, 1, 240, segy_file->fp);
	if (bytes_read == 0)
		return 0;

	if (bytes_read != 240) {
		printf(
				"Fatal error: End of file found too early, cannot read the trace header:\n"
						"expected 240 bytes but found only %d bytes.\n",
				bytes_read);
		return 1;
	}

	flip_trace_header_endianess(segy_file);

	n_samples = GET_SEGYTRACEH_Number_of_samples_in_this_trace(
				&segy_file->trace_header);
	if (n_samples == 0)
		prerror_and_exit(
				"Warning! the number of samples within the trace header is set to 0\n"
						"this can cause troubles, you can override that with the -samples_per_trace switch.\n");
	if (n_samples < 0) {
		fprintf(stderr,
				"WARNING: the field NUMBER_OF_SAMPLES_IN_THIS_TRACE is below zero,\n"
				"         trying to continue casting the value to unsigned short ...\n");
		n_samples = GET_SEGYTRACEH_USHORT_Number_of_samples_in_this_trace(
				&segy_file->trace_header);
		fprintf(stderr,
				"         NUMBER_OF_SAMPLES_IN_THIS_TRACE now equals to: %d\n",
				n_samples);
	}

	long record_nr = GET_SEGYTRACEH_Original_field_record_number(&segy_file->trace_header);
	long tr_nr_in_reel = GET_SEGYTRACEH_Trace_sequence_number_within_reel(&segy_file->trace_header);
	long tr_nr_in_record = GET_SEGYTRACEH_Trace_number_within_field_record(&segy_file->trace_header);
	long tr_nr_in_line = GET_SEGYTRACEH_Trace_sequence_number_within_line(&segy_file->trace_header);
	long sample_interval = GET_SEGYTRACEH_Sample_interval(&segy_file->trace_header);

	if ((verbose == 2 || verbose == 3) && !feof(segy_file->fp)
			&& keep_trace(segy_file, rec_start, rec_end, trace_start,
					trace_end, trace_offset)) {
		printf("--------------------------- TRACE HEADER ------------------\n");
		printf("Original field record number : %ld\n", record_nr);
		printf("Trace number within field record: %ld\n", tr_nr_in_record);
		printf("Trace sequence number within line : %ld\n", tr_nr_in_line);
		printf("Trace sequence number within reel : %ld\n", tr_nr_in_reel);
		printf("Number of samples in this trace: %d\n", n_samples);
		printf("Sample interval (microseconds) : %ld\n", sample_interval);
	}

	segy_file->trace_data = (unsigned char *) realloc(segy_file->trace_data,
			trace_data_length(segy_file));
	segy_file->trace_data_double = (double *) realloc(
			segy_file->trace_data_double, n_samples * sizeof(double));

	if ((bytes_read = fread(segy_file->trace_data, trace_data_length(segy_file),
			1, segy_file->fp)) != 1) {
		printf(
				"Fatal error: End of file found too early, cannot read the trace data.\n");
		return 1;
	}

	flip_trace_data_endianess(segy_file);

	/* skip_n_samples....
	 */

	if(skip_nsamples != -1 || only_nsamples != -1)
	{
		int start, end;
		start = 0;
		if (skip_nsamples != -1)
			start = skip_nsamples;
		end = n_samples;
		if (only_nsamples != -1)
			end = start + only_nsamples;
		if(end > n_samples) end = n_samples;

		int sl = trace_sample_length(segy_file);
		int ns = end - start;

		int ii = 0;
		for(int i = start; i < end; i++)
			segy_file->trace_data[ii++ * sl] = segy_file->trace_data[i * sl];

		/* Change num_samples and delay time accordingly
		 */
		set_short(ns, (void*)(&segy_file->header) + 3220);
		set_short(ns, (void*)(&segy_file->trace_header) + 114);
		set_short((skip_nsamples * sample_interval) / 1000 + get_short((void*)(&segy_file->trace_header) + 108), (void*)(&segy_file->trace_header) + 108);
		n_samples = ns;
	}
	delay_time = GET_SEGYTRACEH_Delay_time_between_source_and_recording_time(&segy_file->trace_header);

	if (verbose == 2 && keep_trace(segy_file, rec_start, rec_end, trace_start, trace_end, trace_offset))
		printf("Rec/Seq/Num = %ld/%ld/%ld\n", record_nr, tr_nr_in_reel, tr_nr_in_record);

	if (print_rec_seq_num &&
			keep_trace(segy_file, rec_start, rec_end, trace_start,	trace_end, trace_offset))
		printf("%ld %ld %ld\n", record_nr, tr_nr_in_reel, tr_nr_in_record);

	total_traces++;

	if(min_num_samples > n_samples) {
		min_num_samples = n_samples;
		min_num_samples_rec_num = record_nr;
		min_num_samples_trace_num = tr_nr_in_record;
	}

	if(max_num_samples < n_samples) {
		max_num_samples = n_samples;
		max_num_samples_rec_num = record_nr;
		max_num_samples_trace_num = tr_nr_in_record;
	}

	if(prev_record_nr != record_nr)
	{
		total_records++;
		prev_record_nr = record_nr;
	}

	/* Decode the trace data in double.
	 */
	if (dump || plot_data || scan || apply_correction || enable_X11) {
		for (i = 0; i < n_samples; i++) {
			segy_file->trace_data_double[i] = get_val(segy_file,
					segy_file->trace_data, i);
			if (trace_min_val > segy_file->trace_data_double[i])
				trace_min_val = segy_file->trace_data_double[i];
			if (trace_max_val < segy_file->trace_data_double[i])
				trace_max_val = segy_file->trace_data_double[i];
		}
	}

	double x, y;
	int uom;
	if (dump_xy
			&& keep_trace(segy_file, rec_start, rec_end, trace_start,
					trace_end, trace_offset)) {
		printf("%d %d %d ",
				GET_SEGYTRACEH_Original_field_record_number(
						&segy_file->trace_header),
				GET_SEGYTRACEH_Trace_sequence_number_within_reel(
						&segy_file->trace_header),
				GET_SEGYTRACEH_Trace_number_within_field_record(
						&segy_file->trace_header));

		if (source_1_or_receiver_2 == 1) {
			get_source_xy(segy_file, &x, &y, &uom);
			printf("%lf %lf", x, y);
		}

		if (source_1_or_receiver_2 == 2) {
			get_receiver_xy(segy_file, &x, &y, &uom);
			printf("%lf %lf", x, y);
		}
		switch (uom) {
		case 2:
			printf(" arcsec\n");
			break;
		case 1:
			switch (GET_SEGYH_Measurements_system(&segy_file->header)) {
			case 1:
				printf(" feet\n");
				break;
			case 2:
				printf(" meters\n");
				break;
			default:
				printf(" unit_of_measure_not_set\n");
				break;
			}
			break;
		default:
			printf(" unit_of_measure_not_set\n");
			break;
		}
	}

	if (dump_fields) {
		int field_segy_offset, jjj;
		char field_segy_type;
		int m_field_nr = 1;
		char m_field[1000], m_fields[1000], value[1000];

		printf("Rec/Seq/Num = %d/%d/%d : fields = ",
				GET_SEGYTRACEH_Original_field_record_number(
						&segy_file->trace_header),
				GET_SEGYTRACEH_Trace_sequence_number_within_reel(
						&segy_file->trace_header),
				GET_SEGYTRACEH_Trace_number_within_field_record(
						&segy_file->trace_header));
		while (get_field(fields_to_dump, m_field_nr, m_fields, ',')) {
			if (m_field_nr > 1)
				printf("; ");
			get_field(m_fields, 1, m_field, ':');
			jjj = get_parameter_index_by_offset(trace_header_types,
					atoi(m_field));
			field_segy_offset = atoi(m_field);
			get_field(m_fields, 2, m_field, ':');
			field_segy_type = m_field[0];
			if (use_names)
				printf("%s,%s", trace_header_names[jjj],
						get_str_val(segy_file->trace_header.HEADER, value,
								field_segy_offset, field_segy_type));
			else
				printf("%d,%c,%s", field_segy_offset, field_segy_type,
						get_str_val(segy_file->trace_header.HEADER, value,
								field_segy_offset, field_segy_type));
			m_field_nr++;
		}

		printf("\n");
	}

	processed_traces++;
	return 0;
}

int do_dump_trace() {
	int i;
	if (!dump)
		return 1;

	printf("Rec/Seq/Num = %d/%d/%d : ",
			GET_SEGYTRACEH_Original_field_record_number(
					&segy_file.trace_header),
			GET_SEGYTRACEH_Trace_sequence_number_within_reel(
					&segy_file.trace_header),
			GET_SEGYTRACEH_Trace_number_within_field_record(
					&segy_file.trace_header));

	printf("%d\n", n_samples);
	for (i = 0; i < n_samples; i++) {
		printf("%lf,", segy_file.trace_data_double[i]);
	}
	printf("\n");

	return 0;
}

/* initialize all variables.
 */
void setup() {
	total_records = total_traces = max_num_samples = 0;
	min_num_samples = 9999999;
	min_num_samples_rec_num = max_num_samples_rec_num = -1;
	min_num_samples_trace_num = max_num_samples_trace_num = -1;
	prev_record_nr = -1;
	vertical_stack = all_file = i = j = k = source_1_or_receiver_2 = 0;
	initial_seek = rec_start = trace_start = 0;
	trace_end = rec_end = output_segy = verbose = n_traces = n_samples =
			reccnt = 0;
	dump = false;
	print_rec_seq_num = dump_header_fields = change_header_fields = add_xy =
	false;
	dump_xy = shot_renumber = trace_renumber = false;
	processed_traces = use_names = 0;
	enable_X11 = no_header = plot_data = false;
	trace_scale = 1.0;
	page_format = 4;
	trace_min_val = 1e100;
	trace_max_val = -1e100;
	scan = 0;
	correction_op = '+';
	correction_val = 0.0;

	convert_to = ' ';
	count = 1;

	fields_to_dump = fields_to_change_fname = add_coordinates_fname = NULL;
	only_traces_with = dump_fields = change_fields = false;

	out_segy_file.trace_data = NULL;
	out_segy_file.trace_data_double = NULL;
	segy_file.trace_data = NULL;
	segy_file.trace_data_double = NULL;
	initial_record = 1;
	initial_trace_seq = 1;
	actual_row = 1;
	actual_record = -1;
	memset(my_ebcdic, 0, 3200);
}

/* process command line args
 */
void read_args(int argc, char **argv) {
	if ((_n = take_parm(argc, argv, "-use_names", 0)))
		remove_parms(&argc, argv, _n, 1),
		use_names = true;

	if ((_n = take_parm(argc, argv, "-segy_info", 0))) {
		print_segy_info();
		exit(0);
	}

	if (take_parm(argc, argv, "-f", 0) == 0) {
		print_usage(argc, argv);
		exit(-1);
	}

	if ((_n = take_parm(argc, argv, "-v", 1)))
	{
		verbose = atoi(argv[_n + 1]);
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-view", 0)))
	{
		enable_X11 = true;
		remove_parms(&argc, argv, _n, 1);
#ifndef WITH_SDL
		prerror_and_exit("Error, SDL2 support required for -view switch.\nPlease recompile the program with -DWITH_SDL define.\nand make sure that SDL2 and SDL2_gfx libraries and headers are properly installed.\n");
#endif
	}

	if ((_n = take_parm(argc, argv, "-scan", 0)))
		remove_parms(&argc, argv, _n, 1),
		scan = true;

	flip_endianess = false;
	if ((_n = take_parm(argc, argv, "-flip_endianess", 0))) {
		remove_parms(&argc, argv, _n, 1),
		flip_endianess = true;
	}

	if ((_n = take_parm(argc, argv, "-traces_per_record", 1)))
	{
		n_traces = atoi(argv[_n + 1]);
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-samples_per_trace", 1)))
	{
		n_samples = atoi(argv[_n + 1]);
	remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-no_header", 0)))
		remove_parms(&argc, argv, _n, 1),
		no_header = true;

	if ((_n = take_parm(argc, argv, "-x", 1)))
	{
		initial_seek = atoll(argv[_n + 1]);
		remove_parms(&argc, argv, _n, 2);
	}

	skip_ntraces = 0;
	if ((_n = take_parm(argc, argv, "-skip_n_traces", 1))) {
		skip_ntraces = atoi(argv[_n + 1]);
		remove_parms(&argc, argv, _n, 2);
	}

	only_ntraces = -1;
	if ((_n = take_parm(argc, argv, "-only_n_traces", 1))) {
		only_ntraces = atoi(argv[_n + 1]);
		remove_parms(&argc, argv, _n, 2);
	}

	skip_nsamples = -1;
	if ((_n = take_parm(argc, argv, "-skip_n_samples", 1))) {
		skip_nsamples = atoi(argv[_n + 1]);
		remove_parms(&argc, argv, _n, 2);
	}

	only_nsamples = -1;
	if ((_n = take_parm(argc, argv, "-only_n_samples", 1))) {
		only_nsamples = atoi(argv[_n + 1]);
		remove_parms(&argc, argv, _n, 2);
	}

	if (take_parm(argc, argv, "-dump_header_fields", 1)
			&& take_parm(argc, argv, "-dump_trace_fields", 1)) {
		prerror_and_exit(
				"Error: -dump_header_fields and -dump_trace_fields are not allowed together.\n");
	}

	if ((_n = take_parm(argc, argv, "-do_op", 1))) {
		apply_correction = true;
		correction_op = argv[_n + 1][0];
		correction_val = atof(argv[_n + 1] + 2);
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-dump_header_fields", 1))) {
		dump_header_fields = true;
		if (use_names) {
			from_names_to_offsets(segy_header_names, argv[_n + 1],
					my_header_fields_offsets);
			header_fields_to_dump = my_header_fields_offsets;
		} else
			header_fields_to_dump = argv[_n + 1];
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-print_rec_seq_num", 0))) {
		remove_parms(&argc, argv, _n, 1),
		print_rec_seq_num = true;
	}

	if ((_n = take_parm(argc, argv, "-change_header_fields", 1))) {
		change_header_fields = true;
		if (use_names) {
			from_names_to_offsets(segy_header_names, argv[_n + 1],
					my_header_fields_offsets);
			header_fields_to_change = my_header_fields_offsets;
		} else
			header_fields_to_change = argv[_n + 1];
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-only_traces_with", 1))) {
		only_traces_with = true;
		if (use_names) {
			from_names_to_offsets(trace_header_names, argv[_n + 1],
					my_valid_traces_fields_offsets);
			traces_fields_valid_values = my_valid_traces_fields_offsets;
		} else
			traces_fields_valid_values = argv[_n + 1];
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-dump_trace_fields", 1))) {
		dump_fields = true;
		if (use_names) {
			from_names_to_offsets(trace_header_names, argv[_n + 1],
					my_traces_fields_offsets);
			fields_to_dump = my_traces_fields_offsets;
		} else
			fields_to_dump = argv[_n + 1];
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-change_trace_fields", 1))) {
		change_fields = true;
		fields_to_change_fname = argv[_n + 1];
		file_of_fields = fopen(fields_to_change_fname, "r");
		if (file_of_fields == NULL)
			prerror_and_exit("Cannot open \"%s\" file, aborting.",
					fields_to_change_fname);
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-EBCDIC", 1))) {
		FILE *f;
		f = fopen(argv[_n + 1], "r");
		if (f == NULL) {
			prerror_and_exit("Cannot open \"%s\" file, aborting.",
					argv[_n + 1]);
		}
		int ii = 0;
		char c;
		while (!feof(f) && ii < 3200) {
			c = fgetc(f);
			if (c != '\n' && c != '\r')
				my_ebcdic[ii++] = ascii2ebcdic[(int)c];
		}
		replace_ebcdic = true;
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-dump_xy", 1))) {
		dump_xy = true;
		source_1_or_receiver_2 = 0;
		if (strcmp(argv[_n + 1], "SOURCE") == 0)
			source_1_or_receiver_2 = 1;
		if (strcmp(argv[_n + 1], "RECEIVER") == 0)
			source_1_or_receiver_2 = 2;
		if (source_1_or_receiver_2 != 1 && source_1_or_receiver_2 != 2)
			prerror_and_exit(
					"either SOURCE or RECEIVER must be given with the switch -dump_xy.\n");
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-add_xy", 1))) {
		add_xy = true;
		add_coordinates_fname = malloc(strlen(argv[_n + 1]) + 1);
		get_field(argv[_n + 1], 1, add_coordinates_fname, ',');
		char *tmp = malloc(strlen(argv[_n + 1]) + 1);
		get_field(argv[_n + 1], 2, tmp, ',');
		source_1_or_receiver_2 = 0;
		if (strcmp(tmp, "SOURCE") == 0)
			source_1_or_receiver_2 = 1;
		if (strcmp(tmp, "RECEIVER") == 0)
			source_1_or_receiver_2 = 2;
		if (source_1_or_receiver_2 != 1 && source_1_or_receiver_2 != 2)
			prerror_and_exit(
					"either SOURCE or RECEIVER must be given with the switch -add_xy.\n");
		read_xy(add_coordinates_fname);
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-irc", 1))) {
		shot_renumber = true;
		initial_record = atol(argv[_n + 1]);
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-vertical_stack", 1))) {
		vertical_stack = atol(argv[_n + 1]);
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-do_ps", 1))) {
		plot_data = true;
		char *pf = argv[_n + 1];
		char pf1[100];
		if (!get_field(pf, 1, pf1, ','))
			prerror_and_exit(
					"FATAL ERROR: -do_ps require additional arguments.\n");
		if (pf1[0] == 'A') {
			if (pf1[1] == 'C')
			{
				page_format = 5; // CUSTOM SIZE
				char pf2[100], pf3[100];
				get_field(pf1+3, 1, pf2, 'x');
				get_field(pf1+3, 2, pf3, 'x');
				if(pf2 == NULL || pf3 == NULL)
					prerror_and_exit(
							"FATAL ERROR: page dimension error.\n");
				paper_size_x[page_format] = atoi(pf2);
				paper_size_y[page_format] = atoi(pf3);
			}
			else
				page_format = atoi(pf1 + 1);
		} else
			prerror_and_exit("FATAL ERROR: Unknown page size.\n");
		if (!get_field(pf, 2, pf1, ','))
			prerror_and_exit(
					"FATAL ERROR: Number of traces to plot not given.\n");
		num_traces_per_cm = atoi(pf1);
		if (!get_field(pf, 3, pf1, ','))
			prerror_and_exit("FATAL ERROR: Scale factor not given.\n");
		trace_scale = atof(pf1);
		remove_parms(&argc, argv, _n, 2);
	}

	no_EBCDIC_stamp = false;
	if ((_n = take_parm(argc, argv, "-no_EBCDIC_stamp", 0))) {
		remove_parms(&argc, argv, _n, 1),
		no_EBCDIC_stamp = true;
	}

	if ((_n = take_parm(argc, argv, "-itc", 1))) {
		trace_renumber = true;
		initial_trace_seq = atol(argv[_n + 1]);
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-convert", 1)))
	{
		convert_to = argv[_n + 1][0];
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-dump", 0)))
	{
		dump = 1;
		remove_parms(&argc, argv, _n, 1);
	}

	if ((_n = take_parm(argc, argv, "-info", 0))) {
		if ((_n = take_parm(argc, argv, "-f", 1)))
			if (open_segy(&segy_file, argv[_n + 1],
					"rb", initial_seek))
				prerror_and_exit("Cannot open input file.\n");

		get_segy_header(&segy_file, 2);
		exit(0);
	}

	output_segy = 0;
	if ((_n = take_parm(argc, argv, "-o", 1))) {
		output_segy = 1;
		if (open_segy(&out_segy_file, argv[_n + 1],
				"wb", 0)) {
			prerror_and_exit("Cannot open output file.\n");
			exit(-1);
		}
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-trace", 2))) {
		trace_start = atoi(argv[_n + 1]);
		trace_end = atoi(argv[_n + 2]);
		remove_parms(&argc, argv, _n, 3);
	}

	trace_offset = 12; /* Trace_number_within_field_record */
	if ((_n = take_parm(argc, argv, "-num_trace_offset", 1))) {
		trace_offset = atoi(
				argv[_n + 1]);
		remove_parms(&argc, argv, _n, 2);
	}

	if ((_n = take_parm(argc, argv, "-record", 2))) {
		rec_start = atoi(argv[_n + 1]);
		rec_end = atoi(argv[_n + 2]);
		remove_parms(&argc, argv, _n, 3);
	}

	if (rec_start || rec_end || trace_start || trace_end) {
		if (trace_start == 0 && trace_end == 0) {
			trace_start = 0;
			trace_end = 99999999;
		}

		if (rec_start == 0 && rec_end == 0) {
			rec_start = 0;
			rec_end = 99999999;
		}
	}
	else
	{
		all_file = 1;
		rec_start = trace_start = 0;
		rec_end = trace_end = 99999999;
	}

	initial_seek = 0;
	if ((_n = take_parm(argc, argv, "-f", 1))) {
		if (open_segy(&segy_file, argv[_n + 1], "rb",
				initial_seek)) {
			prerror_and_exit("Cannot open input file.\n");
			exit(-1);
		}
		remove_parms(&argc, argv, _n, 2);
	}
	if(argc > 1)
	{
		printf("Error, the following given args are unknown:\n");
		for(int i = 1; i < argc; i++) printf("%s\n", argv[i]);
		exit(-1);
	}
}

void do_copy_header() {
	if (output_segy)
		copy_segy_header(&segy_file, &out_segy_file);
}

void do_write_header() {
	if (output_segy && !no_header) {
		write_segy_header(&out_segy_file);
		fflush(out_segy_file.fp);
	}
}

void do_change_trace() {
	if (output_segy) {
		copy_segy_trace_header(&segy_file, &out_segy_file);

		if (vertical_stack > 1 && output_segy) {
			set_short(
					get_short((void*) &out_segy_file.trace_header.HEADER + 116)
							* vertical_stack,
					(void*) &out_segy_file.trace_header.HEADER + 116);
			set_short(
					get_short((void*) &out_segy_file.trace_header.HEADER + 114)
							/ vertical_stack,
					(void*) &out_segy_file.trace_header.HEADER + 114);
		}

		if ((count - 1) % n_traces == 0) {
			current_record++;
			current_trace = initial_trace_seq;
		}

		if (shot_renumber)
			set_int(current_record, out_segy_file.trace_header.HEADER + 8);

		if (trace_renumber)
			set_int(current_trace, out_segy_file.trace_header.HEADER + 12);

		/* Add coordinates if needed.
		 */
		if (add_xy) {
			long rec, seq, num;
			rec = GET_SEGYTRACEH_Original_field_record_number(
					&segy_file.trace_header);
			seq = GET_SEGYTRACEH_Trace_sequence_number_within_reel(
					&segy_file.trace_header);
			num = GET_SEGYTRACEH_Trace_number_within_field_record(
					&segy_file.trace_header);

			// find the correct entry
			int jj;
			for (jj = 0; jj < num_coords; jj++) {
				if (rec == coords[jj].original_field_record
						&& seq == coords[jj].trace_seq_within_reel
						&& num == coords[jj].trace_seq_within_field_record) {
					set_short(
							coords[jj].unit_of_measure_1_feet_or_meters_2_arcsec,
							out_segy_file.trace_header.HEADER + 88); // Unit coordinate system (1-meters, 2-seconds of arc)
					set_short(coordinates_scaling_factor,
							out_segy_file.trace_header.HEADER + 70);
					switch (source_1_or_receiver_2) {
					case 1:
						set_int(round(coords[jj].lon_or_x),
								out_segy_file.trace_header.HEADER + 72);
						set_int(round(coords[jj].lat_or_y),
								out_segy_file.trace_header.HEADER + 76);
						break;
					case 2:
						set_int(round(coords[jj].lon_or_x),
								out_segy_file.trace_header.HEADER + 80);
						set_int(round(coords[jj].lat_or_y),
								out_segy_file.trace_header.HEADER + 84);
						break;
					default:
						prerror_and_exit(
								"Coordinates data corrupted, unknown type\n");
						break;
					}
				}
			}
		}

		/* Change fields if required.
		 */
		if (change_fields) {
			int field_segy_offset, rec_nr_to_change, trace_nr_to_change,
					trace_seq_to_change;
			char field_segy_type;
			char field_segy_value[1000];
			int m_field_nr = 1;
			char m_field[1000], m_fields[1000], line[100000];
			char all_fields[10000];
			char _b0[100000], _b1[100000];

			if (file_of_fields) {
				if (!feof(file_of_fields)) {
					fgets(line, 99999, file_of_fields);
					if (line[strlen(line) - 1] == '\n')
						line[strlen(line) - 1] = 0;

					/* Find rec/seq/trace number.
					 */
					if (get_field(line, 1, _b0, ':') == NULL)
						prerror_and_exit(
								"Error: file_of_fields format error.\n");
					if (get_field(_b0, 2, _b1, '=') == NULL)
						prerror_and_exit(
								"Error: file_of_fields format error.\n");
					if (get_field(_b1, 1, _b0, '/') == NULL)
						prerror_and_exit(
								"Error: file_of_fields format error.\n");
					rec_nr_to_change = atol(_b0);
					if (get_field(_b1, 2, _b0, '/') == NULL)
						prerror_and_exit(
								"Error: file_of_fields format error.\n");
					trace_seq_to_change = atol(_b0);
					if (get_field(_b1, 3, _b0, '/') == NULL)
						prerror_and_exit(
								"Error: file_of_fields format error.\n");
					trace_nr_to_change = atol(_b0);
					if (verbose)
						printf("Changing Rec/Seq/Num = %d/%d/%d :",
								rec_nr_to_change, trace_seq_to_change,
								trace_nr_to_change);

					if (rec_nr_to_change
							!= GET_SEGYTRACEH_Original_field_record_number(
									&out_segy_file.trace_header)
							|| trace_nr_to_change
									!= GET_SEGYTRACEH_Trace_number_within_field_record(
											&out_segy_file.trace_header)
							|| trace_seq_to_change
									!= GET_SEGYTRACEH_Trace_sequence_number_within_reel(
											&out_segy_file.trace_header))
						prerror_and_exit(
								"Error \"change_fields\" file got unsynchronized.");

					/* Find the fields to change and their value/type.
					 */
					if (get_field(line, 2, _b0, ':') == NULL)
						prerror_and_exit(
								"Error: file_of_fields format error.\n");
					if (get_field(_b0, 2, all_fields, '=') == NULL)
						prerror_and_exit(
								"Error: file_of_fields format error.\n");
					m_field_nr = 1;
					if (use_names) {
						while (get_field(all_fields, m_field_nr, m_fields, ';')
								!= NULL) {
							get_field(m_fields, 1, m_field, ',');
							trim(m_field);
							field_segy_offset =
									trace_header_types[get_parameter_index_by_name(
											trace_header_names, m_field) * 2];
							field_segy_type =
									trace_header_types[get_parameter_index_by_name(
											trace_header_names, m_field) * 2 + 1];
							get_field(m_fields, 2, m_field, ',');
							strcpy(field_segy_value, m_field);
							if (verbose)
								printf("%d,%c,%s; ", field_segy_offset,
										field_segy_type, field_segy_value);
							set_str_val(out_segy_file.trace_header.HEADER,
									field_segy_value, field_segy_offset,
									field_segy_type);
							m_field_nr++;
						}
					} else {
						while (get_field(all_fields, m_field_nr, m_fields, ';')
								!= NULL) {
							get_field(m_fields, 1, m_field, ',');
							field_segy_offset = atoi(m_field);
							get_field(m_fields, 2, m_field, ',');
							field_segy_type = m_field[0];
							get_field(m_fields, 3, m_field, ',');
							strcpy(field_segy_value, m_field);
							if (verbose)
								printf("%d,%c,%s; ", field_segy_offset,
										field_segy_type, field_segy_value);
							set_str_val(out_segy_file.trace_header.HEADER,
									field_segy_value, field_segy_offset,
									field_segy_type);
							m_field_nr++;
						}
					}
					if (verbose)
						printf("\n");
				}
			}
		}
		copy_segy_trace_data(&segy_file, &out_segy_file);
		int out_n_samples = GET_SEGYTRACEH_Number_of_samples_in_this_trace(
				&out_segy_file.trace_header);
		if (apply_correction) {
			for (i = 0; i < out_n_samples; i++) {
				switch (correction_op) {
				case '*':
					SET_OUT_SAMPLE(GET_OUT_SAMPLE(i) * correction_val, i);
					break;
				case '/':
					SET_OUT_SAMPLE(GET_OUT_SAMPLE(i) / correction_val, i);
					break;
				case '-':
					SET_OUT_SAMPLE(GET_OUT_SAMPLE(i) - correction_val, i);
					break;
				case '+':
					SET_OUT_SAMPLE(GET_OUT_SAMPLE(i) + correction_val, i);
					break;
				}
			}
		}
	}
}

void do_write_trace() {
	if (output_segy) {
		write_segy_trace_header(&out_segy_file);
		write_segy_trace_data(&out_segy_file);
		fflush(out_segy_file.fp);
		count++;
		current_trace++;
	}
}

int more_data() {
	if (only_ntraces == -1)
		return !feof(segy_file.fp);
	return only_ntraces >= processed_traces;
}

void do_close_files() {
	if (segy_file.fp && segy_file.fp != stdin)
		fclose(segy_file.fp);
	if (output_segy && out_segy_file.fp && out_segy_file.fp != stdout)
		fclose(out_segy_file.fp);
	if (plot_data) {
		printf("showpage\n");
	}
}

void do_plot_shots() {
	if (!plot_data)
		return;

	int i;

	double num_samples = n_samples;
	double sample_interval = GET_SEGYTRACEH_Sample_interval(&segy_file.trace_header);
	double total_pixel_length = paper_size_y[page_format];
	double total_pixel_height = paper_size_x[page_format];

	double offset_y = 12.0;
	double offset_x = 12.0;
	double num_traces = (double)num_traces_per_cm * (double)(total_pixel_length - 2 * offset_x) / 10.0;
	if(actual_row > num_traces) return;
	double scale_y = (total_pixel_height - 2 * offset_y) / num_samples;
	double scale_x = (total_pixel_length - 2 * offset_x) / num_traces;
	double _x, _y;


	if (actual_row == 1) {
		printf("%%!PS-Adobe-3.0 EPSF-3.0\n");
		printf("%%BoundingBox: %.1f %.1f %.1f %.1f\n", 0.0, 0.0,
				total_pixel_length, total_pixel_height);
		printf("%%Pages: 1\n");
		printf("%%%% Generated by segy-change\n");
		printf("/l {lineto} def\n");
		printf("/centershow\n");
		printf("{ dup stringwidth pop\n");
		printf("-2 div\n");
		printf("0 rmoveto\n");
		printf("show } def\n");
		printf("/Times-Roman findfont\n");
		printf("4 scalefont\n");
		printf("setfont\n");
		printf("newpath\n");
		printf("0.5 setlinewidth\n");
		printf("%f %f scale\n", 0.039370079 * 72.0, 0.039370079 * 72.0);
		printf("%f 0 translate\n", total_pixel_height);
		printf("%d rotate\n", 90);
		printf("%.1f %.1f moveto\n", (total_pixel_length / 2.0),
				total_pixel_height - offset_y + 5);
		printf("(%s) centershow\n", segy_file.fname);
		printf("%.1f %.1f moveto\n", (total_pixel_length / 2.0), offset_y - 10);
		printf("(%s) centershow\n", "Record number");
		printf("closepath\n");
		printf("newpath\n");

		printf("3.5 %.1f translate\n", (total_pixel_height / 2.0));
		printf("90 rotate\n");
		printf("0 0 moveto\n");
		printf("(%s) centershow\n", "Time (ms)");
		printf("-90 rotate\n");
		printf("-3.5 %.1f translate\n", -(total_pixel_height / 2.0));
		printf("closepath\n");


		printf("newpath\n");
		printf("%.2f %.2f moveto\n", offset_x, offset_y);
		printf("%.2f %.2f lineto\n", (total_pixel_length - offset_x), offset_y);
		printf("%.2f %.2f lineto\n", (total_pixel_length - offset_x),
				(total_pixel_height - offset_y));
		printf("%.2f %.2f lineto\n", offset_x, (total_pixel_height - offset_y));
		printf("%.2f %.2f lineto\n", offset_x, offset_y);
		printf("closepath\n");
		printf("stroke\n");
		printf("0.05 setlinewidth\n");
		printf("/Times-Roman findfont\n");
		printf("2 scalefont\n");
		printf("setfont\n");

		int st = num_samples / 10;
		st = (st + 5) / 10;
		st = st * 10;
		if(st < 1) st = 1;
		for(int j = 0; j <= num_samples; j += st)
		{
			double _y = j * scale_y;
			printf("newpath\n");
			printf("%.2f %.2f moveto\n", offset_x - 2.5, (total_pixel_height - _y - offset_y));
			printf("%.2f %.2f lineto\n", offset_x, (total_pixel_height - _y - offset_y));
			printf("closepath\nstroke\n");

			printf("newpath\n");
			printf("%.2f %.2f moveto\n", offset_x - 5, (total_pixel_height - _y - offset_y));
			printf("(%.2f) centershow\n", (double)delay_time + ((double)j * (double)sample_interval) / 1000.0);
		}
	}

	int num_record = GET_SEGYTRACEH_Original_field_record_number(
			&segy_file.trace_header);
	int every_one_record = ((num_traces / 10) + 5) / 10;
	every_one_record *= 10;
	if(every_one_record == 0) every_one_record = 1;

	if ((((num_record % every_one_record) == 0 || every_one_record == 0)
			&& num_record != actual_record) || actual_row == 0) {
		printf("newpath\n");
		printf("%.2f %.2f moveto\n", (offset_x + actual_row * scale_x),
				offset_y);
		printf("%.2f %.2f lineto\n", (offset_x + actual_row * scale_x),
				offset_y - 2.5);
		printf("closepath\nstroke\n");

		printf("newpath\n");
		printf("%.2f %.2f moveto\n", (offset_x + actual_row * scale_x),
				offset_y - 4);
		printf("(%d) centershow\n", num_record);
	}

	printf("newpath\n");
	printf("%.2f %.2f moveto\n", (offset_x + actual_row * scale_x),
			(total_pixel_height - offset_y));

	double y, step = 1;
	//if(num_samples > 2000) step = num_samples / 2000;
	for (y = 0; y < num_samples; y += step) {
		i = y;
		_x = segy_file.trace_data_double[i] * trace_scale;
		_y = i * scale_y;

		if (fabs(_x) > total_pixel_length / num_traces)
			_x = my_sgn(_x) * total_pixel_length / num_traces;
		if (_x >= 0)
			printf("%.2f %.2f %.2f %.2f rectfill\n", (offset_x + (actual_row - 1) * scale_x),
								(total_pixel_height - _y - offset_y),
								_x,
								scale_y * step);

//			printf("%.2f %.2f l\n", (offset_x + actual_row * scale_x),
//					(total_pixel_height - _y - offset_y));
//		printf("%.2f %.2f l\n", (offset_x + _x + actual_row * scale_x),
//				(total_pixel_height - _y - offset_y));

	}
	printf("closepath\n");
	printf("stroke\n");
	actual_row++;
	actual_record = num_record;
}

//######################### X11         E X T E N S I O N

#ifdef WITH_SDL
#define BPP 4
#define WIDTH 1400
#define HEIGHT 600
#define BORDER_Y_UP 24
#define BORDER_Y_DOWN 60
#define BORDER_X_LEFT 100
#define BORDER_X_RIGHT 32

//SDL_Surface *screen, *surface;
SDL_Window *sdlWindow;
SDL_Renderer *sdlRenderer;
SDL_Event event;
SDL_Surface *surface;
SDL_Texture *bitmapTex = NULL;
Uint8 *keys;

double max_val_wave, min_val_wave, max_val_proc_wave, min_val_proc_wave,
		num_samples;
double scale;
int zoom_wave, cur_palette = 1;
unsigned int height, width, _dummy;
char text[256];
int pos_x, pos_y;
int s;
bool exposed;
bool rescale;
int mouse_x, mouse_y, mouse_last_x, mouse_last_y;
#define MAX_WAVES 100000

double *waves[MAX_WAVES];
SEGY_trace_header trace_headers[MAX_WAVES];

Uint32 *waves_pixmap_lines = NULL;
long waves_num_samples[MAX_WAVES];
double waves_sample_interval[MAX_WAVES];

long num_waves = 0;
long max_num_samples_per_wave = 0;
long first_wave = 0;
long first_sample = 0;
double zoom_step_x = 1;
double zoom_step_y = 1;
short ct = 0xa0;

enum PALETTE_TYPE {
	REDBLUE = 0,
	REDBLACK = 1,
	REDBLACKHC = 2,
	LASTPALETTEINDEX = 3 // Keep this to last_index + 1 !!
};

Uint32 DisplayData_GetRGBAColor(SDL_Surface *surface, int palette_type,
		double value) {
	Uint8 r, g, b;
	switch (palette_type) {
	case REDBLUE:
		if (value > 0) {
			if (value > 0xff)
				value = 0xff;
			if (value == 0)
				r = g = b = 0xff;
			else if (value == 0xff)
				r = g = b = 0;
			else {
				r = 0xff - value;
				g = 0xff - value;
				b = 0xff;
			}
		} else {
			value = abs(value);
			if (value > 0xff)
				value = 0xff;
			if (value == 0)
				r = g = b = 0xff;
			else if (value == 0xff)
				r = g = b = 0;
			else
				r = 0xff;
			g = 0xff - value;
			b = 0xff - value;
		}
		return SDL_MapRGBA(surface->format, r, g, b, 255);
		break;
	case REDBLACK:
		if (value > 0) {
			if (value > 0xff)
				value = 0xff;
			if (value == 0)
				r = g = b = 0xff;
			else if (value == 0xff)
				r = g = b = 0;
			else {
				r = 0xff - value;
				g = 0xff - value;
				b = 0xff - value;
			}
		} else {
			value = abs(value);
			if (value > 0xff)
				value = 0xff;
			if (value == 0)
				r = g = b = 0xff;
			else if (value == 0xff)
				r = g = b = 0;
			else
				r = 0xff;
			g = 0xff - value;
			b = 0xff - value;
		}
		return SDL_MapRGBA(surface->format, r, g, b, 255);
		break;
	case REDBLACKHC:
		if (value > 0) {
			if (value > ct) {
				r = g = b = 0;
			}
			else {
				double v = (ct - value) * 255.0 / ct;
				r = v;
				g = v;
				b = v;
			}
		}
		else {
			double col = 0xb0;
			value = abs(value);
			if (value > ct) {
				g = b = 0;
				r = col;
			}
			else {
			    double v = (ct - value) * 255.0 / ct;
			    double v0 = (ct - value) * (255.0 - col) / ct;
			    g = b = v;
			    r = col + v0;
			}
		}
		return SDL_MapRGBA(surface->format, r, g, b, 255);
		break;
	default:
		return 0;
		break;

	}
}

// Update the pixmap vertical lines according to wave values if line == -1 updates all lines.
void DisplayData_UpdatePixmapLines(long line) {
	long i, j;
	double scale_pix = 255.0 / ((double) (max_val_wave - min_val_wave) / 2);

	waves_pixmap_lines = (Uint32*) realloc(waves_pixmap_lines,
			sizeof(Uint32) * num_waves * max_num_samples_per_wave);

	long bias = (max_val_wave + min_val_wave) / 2;
	bias = 0; // I don't like de-bias

	if (line == -1) {
#pragma omp parallel private(i,j)
#pragma omp for schedule(dynamic, 100)
		for (i = 0; i < num_waves; i++) {
			for (j = 0; j < waves_num_samples[i]; j++) {
				waves_pixmap_lines[i + j * num_waves] =
						DisplayData_GetRGBAColor(surface, cur_palette,
								(scale * (waves[i][j] - bias)) * scale_pix);
			}
		}
	} else {
		for (j = 0; j < waves_num_samples[line]; j++) {
			waves_pixmap_lines[line + j * num_waves] = DisplayData_GetRGBAColor(
					surface, cur_palette,
					(scale * (waves[line][j] - bias)) * scale_pix);
		}
	}
}

void DisplayData_replaceWave(long index, double *values, int count,
		double sample_interval) {
	int i;
	waves_num_samples[index] = count;
	waves_sample_interval[index] = sample_interval;
	if (max_num_samples_per_wave < count)
		max_num_samples_per_wave = count;
	if (waves[index] == NULL) {
		fprintf(stderr, "Error allocating memory at AddWave function.");
		exit(-1);
	}

	for (i = 0; i < count; i++) {
		waves[index][i] = values[i];
	}
}

void DisplayData_AddWave(double *values, int count, double sample_interval, SEGY_trace_header tr_header) {
	int i;
	int MAX_SAMPLE_VALUE, MIN_SAMPLE_VALUE;

	waves[num_waves] = (double*) malloc(count * sizeof(double));
	waves_num_samples[num_waves] = count;
	trace_headers[num_waves] = tr_header;
	waves_sample_interval[num_waves] = sample_interval;
	if (max_num_samples_per_wave < count)
		max_num_samples_per_wave = count;
	if (waves[num_waves] == NULL) {
		fprintf(stderr,
				"Error allocating memory at DisplayData_AddWave function.");
		exit(-1);
	}

	MAX_SAMPLE_VALUE = values[0] / zoom_wave;
	MIN_SAMPLE_VALUE = values[0] / zoom_wave;
	for (i = 0; i < count; i++) {
		waves[num_waves][i] = values[i];
		if (values[i] / zoom_wave > MAX_SAMPLE_VALUE)
			MAX_SAMPLE_VALUE = values[i] / zoom_wave;
		if (values[i] / zoom_wave < MIN_SAMPLE_VALUE)
			MIN_SAMPLE_VALUE = values[i] / zoom_wave;
		//if (MIN_SAMPLE_VALUE == MAX_SAMPLE_VALUE) return;
		if (min_val_wave > MIN_SAMPLE_VALUE)
			min_val_wave = MIN_SAMPLE_VALUE;
		if (max_val_wave < MAX_SAMPLE_VALUE)
			max_val_wave = MAX_SAMPLE_VALUE;
	}
	num_waves++;
}

void DisplayData_SetPixel(SDL_Surface *surface, long x, long y, Uint32 colour) {
	Uint32 *pixmem32;
	y = y * surface->pitch / BPP;
	pixmem32 = (Uint32*) surface->pixels + y + x;
	*pixmem32 = colour;
}

int cnt = 0;
void DisplayData_DrawWaves() {
	DisplayData_CheckZoomAndShift(false); // The check must be performed two times to take place.
	DisplayData_CheckZoomAndShift(false);

	int i, j, y;
	int _nwaves = num_waves;
	if (_nwaves > width)
		_nwaves = width;

	SDL_FillRect(surface, NULL, 0xffffffff);
	Uint32 *pixmap = (Uint32*) surface->pixels;
#pragma omp parallel private(i,j)
#pragma omp for schedule(dynamic, 1000)
	for (j = BORDER_Y_UP; j < height - BORDER_Y_DOWN; j++) {
		double hy = ((long) ((double) (j - BORDER_Y_UP) * zoom_step_y
				+ (double) first_sample)) * (double) num_waves;
		y = j * surface->pitch / BPP;
		for (i = BORDER_X_LEFT; i < width - BORDER_X_RIGHT; i++) {
			pixmap[i + y] = waves_pixmap_lines[(long) ((double) (i
					- BORDER_X_LEFT) * zoom_step_x + (double) first_wave + hy)];
		}
//		memcpy(pixmap + y, waves_pixmap_lines[j * num_waves + first_wave], width * sizeof(Uint32));
	}

	bitmapTex = SDL_CreateTextureFromSurface(sdlRenderer, surface);
	SDL_RenderClear(sdlRenderer);
	SDL_RenderCopy(sdlRenderer, bitmapTex, NULL, NULL);
	Uint32 color = 0xff000000;
	hlineColor(sdlRenderer, BORDER_X_LEFT - 1, width - BORDER_X_RIGHT,
			BORDER_Y_UP - 1, color);
	hlineColor(sdlRenderer, BORDER_X_LEFT - 1, width - BORDER_X_RIGHT,
			height - BORDER_Y_DOWN, color);
	vlineColor(sdlRenderer, BORDER_X_LEFT - 1, BORDER_Y_UP - 1,
			height - BORDER_Y_DOWN, color);
	vlineColor(sdlRenderer, width - BORDER_X_RIGHT, BORDER_Y_UP - 1,
			height - BORDER_Y_DOWN, color);

	int start = BORDER_Y_UP;
	double step = (height - BORDER_Y_DOWN - BORDER_Y_UP - 1) / 10;
	char s[1024];
	gfxPrimitivesSetFont(NULL, 0, 0);
	color = 0x80000000;
	for (j = BORDER_Y_UP + step; j < height - BORDER_Y_DOWN; j += step) {
		double hy = ((long) ((double) (j - BORDER_Y_UP) * zoom_step_y
				+ (double) first_sample));
		sprintf(s, " %8.1lfms", delay_time + 1000 * hy * waves_sample_interval[first_wave]); // The first wave is used to set sample interval
		hlineColor(sdlRenderer, BORDER_X_LEFT, width - BORDER_X_RIGHT, j, color);
		stringColor(sdlRenderer, 0, j - 4, s, 0xff000000);
	}

	step = (height - BORDER_X_RIGHT - BORDER_X_LEFT - 1) / 10;
	for (j = BORDER_X_LEFT + step; j < width - BORDER_X_RIGHT; j += step) {
		int trc = (long) ((double) (j - BORDER_X_LEFT) * zoom_step_x
				+ (double) first_wave);

		sprintf(s, "%d", GET_SEGYTRACEH_Original_field_record_number(&trace_headers[trc])); // The first wave is used to set sample interval
		vlineColor(sdlRenderer, j, BORDER_Y_UP, BORDER_Y_UP - 5, color);

		stringColor(sdlRenderer, j - (strlen(s) * 8) / 2, 10, s, 0xff000000);
	}

	SDL_DestroyTexture(bitmapTex);
	SDL_RenderPresent(sdlRenderer);

}

void DisplayData_CheckZoomAndShift(bool adapt) {
	double h, w;
	h = height - (BORDER_Y_UP + BORDER_Y_DOWN);
	w = width - (BORDER_X_LEFT + BORDER_X_RIGHT);

	if (first_sample > max_num_samples_per_wave - h * zoom_step_y)
		first_sample = max_num_samples_per_wave - h * zoom_step_y;
	if (first_sample < 0 || adapt)
		first_sample = 0;

	if (first_wave > num_waves - w * zoom_step_x)
		first_wave = num_waves - w * zoom_step_x;
	if (first_wave < 0 || adapt)
		first_wave = 0;

	if (zoom_step_y
			> (double) (max_num_samples_per_wave - first_sample) / (double) h
			|| adapt) {
		zoom_step_y = (double) (max_num_samples_per_wave - first_sample)
				/ (double) h;
		//return;
	}
	if (zoom_step_y < 0.1) {
		zoom_step_y = 0.1;
		//return;
	}

	if (zoom_step_x > (double) (num_waves - first_wave) / (double) w || adapt) {
		zoom_step_x = (double) (num_waves - first_wave) / (double) w;
		//return;
	}

	if (zoom_step_x < 0.1) {
		zoom_step_x = 0.1;
		//return;
	}
}

void DisplayData_Init() {
	width = WIDTH;
	height = HEIGHT;
	rescale = false;

	zoom_wave = 1;
	scale = 1.0;
	min_val_wave = 1e300;
	max_val_wave = -1e300;
	mouse_x = mouse_y = mouse_last_x = mouse_last_y = -1;
	num_waves = 0;
	max_num_samples_per_wave = 0;
	first_wave = 0;
	first_sample = 0;
	zoom_step_x = 1;
	zoom_step_y = 1;
}

void DisplayData_Create() {
	SDL_CreateWindowAndRenderer(width, height, SDL_WINDOW_RESIZABLE, &sdlWindow,
			&sdlRenderer);
	SDL_SetWindowTitle(sdlWindow,
			"Segy-change v1.1, active keys: +,-,*,/ & mouse drag, shift + mouse drag");
	surface = SDL_CreateRGBSurface(SDL_SWSURFACE, width, height, 32, 0x00FF0000,
			0x0000FF00, 0x000000FF, 0xFF000000);

	if (sdlWindow == NULL) {
		fprintf(stderr, "Can't init SDL: %s\n", SDL_GetError());
		exit(1);
	}
}

void DisplayData_SetTitle(const char *title) {
	SDL_SetWindowTitle(sdlWindow, title);
}

char DisplayData_GetKey() {
	SDL_PollEvent(&event);
	if (SDL_GetWindowFromID(event.window.windowID) == sdlWindow) {
		switch (event.type) {
		case SDL_KEYDOWN:
			if (event.key.state == SDL_PRESSED) {
//fprintf(stderr, "%d == %d\n", SDLK_q, event.key.keysym.sym);
				switch (event.key.keysym.sym) {
				case SDLK_q:
					return 'q';
				case SDLK_SPACE:
					return ' ';
				case SDLK_GREATER:
					return '>';
				case SDLK_KP_MINUS:
					return '-';
				case SDLK_KP_PLUS:
					return '+';
				}
			}
		}
	}
	return 0;
}

void DisplayData_GetTraceAndSampleAtMouse(long *trace, long *sample) {
	int x, y;
	SDL_GetMouseState(&x, &y);
	*trace = (long) ((double) (x - BORDER_X_LEFT) * zoom_step_x
			+ (double) first_wave);
	*sample = ((long) ((double) (y - BORDER_Y_UP) * zoom_step_y
			+ (double) first_sample));
}

void DisplayData_Yield(bool exit) {
	bool quit = false;
	int x, y, sample, wave;
	while (!quit) {
		if (exit)
			quit = true;
		usleep(1000);
		SDL_PollEvent(&event);
		if (SDL_GetWindowFromID(event.window.windowID) == sdlWindow) {
			switch (event.type) {
			case SDL_WINDOWEVENT:
				switch (event.window.event) {
				case SDL_WINDOWEVENT_CLOSE:
					quit = true;
					while (SDL_PollEvent(&event))
						;
					event.type = 0;
					return;
					break;
				case SDL_WINDOWEVENT_RESIZED:
					width = event.window.data1;
					height = event.window.data2;
					if (surface)
						free(surface);
					surface = SDL_CreateRGBSurface(SDL_SWSURFACE, width, height,
							32, 0x00FF0000, 0x0000FF00, 0x000000FF, 0xFF000000);
					DisplayData_DrawWaves();
					break;
				default:
					break;
				}
				break;
			case SDL_MOUSEWHEEL:
				SDL_GetMouseState(&x, &y);
				sample = first_sample + y * zoom_step_y;
				wave = first_wave + x * zoom_step_x;
				zoom_step_x -= (double) event.wheel.y * 0.1;
				zoom_step_y -= (double) event.wheel.y * 0.1;
				first_sample += ((double) event.wheel.y * zoom_step_y) / 2;
				first_wave += ((double) event.wheel.y * zoom_step_x) / 2;
				DisplayData_DrawWaves();
				break;
			case SDL_MOUSEMOTION:
				if (event.motion.state & SDL_BUTTON(SDL_BUTTON_LMASK)) {
					const Uint8 *state = SDL_GetKeyboardState(NULL);
					if (state[SDL_SCANCODE_LSHIFT]) {
						if (mouse_last_x != -1) {
							zoom_step_x -= (double) (event.motion.x
									- mouse_last_x) * 0.001 * zoom_step_x;
							wave = first_wave + x * zoom_step_x;
						}
						if (mouse_last_y != -1) {
							zoom_step_y -= (double) (event.motion.y
									- mouse_last_y) * 0.001 * zoom_step_y;
							sample = first_sample + y * zoom_step_y;
						}
						//first_sample += ((double) event.motion.y * zoom_step_y) / 2;
						//first_wave += ((double) event.motion.x * zoom_step_x) / 2;
					} else {
						if (mouse_last_x != -1) {
							first_wave -= (event.motion.x - mouse_last_x)
									* zoom_step_x;
						}
						if (mouse_last_y != -1) {
							first_sample -= (event.motion.y - mouse_last_y)
									* zoom_step_y;
						}
					}
					mouse_last_y = event.motion.y;
					mouse_last_x = event.motion.x;
					DisplayData_DrawWaves();

				} else {
					mouse_last_y = -1;
					mouse_last_x = -1;
					long trace, sample;
					DisplayData_GetTraceAndSampleAtMouse(&trace, &sample);
					if(sample < 0) sample = 0;
					if(trace < 0) trace = 0;
					if(trace > num_waves - 1) trace = num_waves - 1;
					if(sample > max_num_samples_per_wave - 1) sample = max_num_samples_per_wave - 1;
					printf("\033[2J\033[1;1H[Record NR, %d, Trace NR=%d, Sample NR=%d]\n", GET_SEGYTRACEH_Original_field_record_number(&trace_headers[trace]), GET_SEGYTRACEH_Trace_number_within_field_record(&trace_headers[trace]), sample);
					printf("Sample value=%lf\n", waves[trace][sample]);
				}

				break;
			case SDL_KEYDOWN:
				if (event.key.state == SDL_PRESSED) {
//fprintf(stderr, "%d == %d\n", SDLK_q, event.key.keysym.sym);
					switch (event.key.keysym.sym) {
					case SDLK_q:
						return;
						break;
					case SDLK_p:
						cur_palette = (cur_palette + 1) % (LASTPALETTEINDEX);
						DisplayData_UpdatePixmapLines(-1);
						DisplayData_DrawWaves();
						break;
					case SDLK_SPACE:
						first_wave = 0;
						first_sample = 0;
						zoom_step_x = zoom_step_y = 1;
						scale = 1;
						DisplayData_UpdatePixmapLines(-1);
						DisplayData_DrawWaves();
						break;
					case SDLK_RIGHT:
						first_wave += 10 / zoom_step_x;
						DisplayData_DrawWaves();
						break;
					case SDLK_LEFT:
						first_wave -= 10 / zoom_step_x;
						DisplayData_DrawWaves();
						break;
					case SDLK_DOWN:
						first_sample += 10 * zoom_step_y;
						DisplayData_DrawWaves();
						break;
					case SDLK_UP:
						first_sample -= 10 * zoom_step_y;
						DisplayData_DrawWaves();
						break;
					case SDLK_KP_MULTIPLY:
						scale *= 1.1;
						DisplayData_UpdatePixmapLines(-1);
						DisplayData_DrawWaves();
						break;
					case SDLK_KP_DIVIDE:
						scale /= 1.1;
						DisplayData_UpdatePixmapLines(-1);
						DisplayData_DrawWaves();
						break;
					case SDLK_PAGEUP:
						zoom_step_x += 0.1;
						DisplayData_DrawWaves();
						break;
					case SDLK_PAGEDOWN:
						zoom_step_x -= 0.1;
						DisplayData_DrawWaves();
						break;
					case SDLK_KP_MINUS:
						zoom_step_y++;
						DisplayData_DrawWaves();
						break;
					case SDLK_KP_PLUS:
						zoom_step_y--;
						DisplayData_DrawWaves();
						break;
					}
				}
			}
			event.type = 0;
			SDL_RenderCopy(sdlRenderer, bitmapTex, NULL, NULL);
			SDL_RenderPresent(sdlRenderer);
		}

	}

	return;
}

#endif

#ifndef SEGYCHANGE_LIBRARY

int main(int argc, char **argv) {
	setup();

	read_args(argc, argv);

#ifdef WITH_SDL
	if (enable_X11) {
		DisplayData_Init();
	}
#endif

	get_segy_header(&segy_file, verbose);

	do_copy_header();

	do_change_header();

	do_write_header();

	do_skip_ntraces();

	current_trace = initial_trace_seq;
	current_record = initial_record - 1;

	while (more_data()) {
		if (get_segy_trace(&segy_file, verbose))
		{
			fprintf(stderr,
					"Error while reading file at byte nr. %ld, record nr. %ld, trace nr. %ld\n",
					ftell(segy_file.fp),
					total_records, total_traces);
			break;
		}
		if (!more_data())
			break;
		if (!keep_trace(&segy_file, rec_start, rec_end, trace_start, trace_end,
				trace_offset))
			continue;
		do_dump_trace();
		do_change_trace();
		do_plot_shots();
		do_write_trace();

#ifdef WITH_SDL
		if (enable_X11) {
			DisplayData_AddWave(segy_file.trace_data_double, n_samples,
					GET_SEGYTRACEH_Sample_interval(&segy_file.trace_header)
							/ 1000000.0, segy_file.trace_header);
		}
#endif
	}

	do_close_files();
	if (scan)
	{
		fprintf(stderr, "Total records = %ld\n", total_records);
		fprintf(stderr, "Total traces = %ld\n", total_traces);
		if(min_num_samples == max_num_samples)
		{
			fprintf(stderr, "num_samples = %ld\n", min_num_samples);
		}
		else
		{
			fprintf(stderr, "min_num_samples = %ld at record nr. %ld, trace nr. %ld\n",
				min_num_samples,
				min_num_samples_rec_num,
				min_num_samples_trace_num);

			fprintf(stderr, "max_num_samples = %ld at record nr. %ld, trace nr. %ld\n",
				max_num_samples,
				max_num_samples_rec_num,
				max_num_samples_trace_num);
		}

		fprintf(stderr, "[min_sample_value, max_sample_value] = [%lf, %lf]\n", trace_min_val, trace_max_val);
	}

#ifdef WITH_SDL
	if (enable_X11) {
		DisplayData_Create();
		DisplayData_UpdatePixmapLines(-1);
		DisplayData_CheckZoomAndShift(true);
		DisplayData_DrawWaves();
		DisplayData_Yield(false);
	}
#endif

	return 0;
}
#endif
