
/**
    COPYRIGHT DISCLAIMER

    Vincent Le Guilloux, Peter Schmidtke and Pierre Tuffery, hereby
	disclaim all copyright interest in the program “fpocket” (which
	performs protein cavity detection) written by Vincent Le Guilloux and Peter
	Schmidtke.

    Vincent Le Guilloux  28 November 2008
    Peter Schmidtke      28 November 2008
    Pierre Tuffery       28 November 2008

    GNU GPL

    This file is part of the fpocket package.

    fpocket is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    fpocket is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with fpocket.  If not, see <http://www.gnu.org/licenses/>.

**/

#ifndef DH_MDPARAMS
#define DH_MDPARAMS

/* ----------------------------- INCLUDES ------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>

#include "utils.h"
#include "memhandler.h"
#include "fparams.h"

/* --------------------------- PUBLIC MACROS ---------------------------------*/


/* Options of the test program */
#define M_MDPAR_INPUT_FILE 'L'
#define M_MDPAR_INPUT_FILE2 'f'
#define M_MDPAR_OUTPUT_FILE 'o'
#define M_MDPAR_SCORING_MODE 'S'
#define M_MDPAR_OUTPUT_ALL_SNAPSHOTS 'a'


#define M_MDP_OUTPUT_FILE1_DEFAULT "mdpout_snapshots_concat.pqr"
#define M_MDP_OUTPUT_FILE2_DEFAULT "mdpout_freq_grid.dx"
#define M_MDP_OUTPUT_FILE3_DEFAULT "mdpout_freq_iso_0_5.pdb"
#define M_MDP_OUTPUT_FILE4_DEFAULT "mdpout_descriptors.txt"
#define M_MDP_OUTPUT_FILE5_DEFAULT "mdpout_mdpocket.pdb"
#define M_MDP_OUTPUT_FILE6_DEFAULT "mdpout_mdpocket_atoms.pdb"
#define M_MDP_OUTPUT_FILE7_DEFAULT "mdpout_all_atom_pdensities.pdb"
#define M_MDP_OUTPUT_FILE8_DEFAULT "mdpout_dens_grid.dx"
#define M_MDP_OUTPUT_FILE9_DEFAULT "mdpout_dens_iso_8.pdb"
#define M_MDP_DEFAULT_ISO_VALUE_FREQ 0.5
#define M_MDP_DEFAULT_ISO_VALUE_DENS 8.0

#define M_MAX_FILE_NAME_LENGTH 300
#define M_NB_MC_ITER 2500
//#define M_MIN_ASPH_RAY 3.0
//#define M_MAX_ASPH_RAY 6.0


#define M_MDP_USAGE "\
***** USAGE (mdpocket) *****\n\
\n\
1 : Pocket finding on a MD trajectory -> list of pre-aligned pdb ordered by time file:\n\
\t./bin/mdpocket -L pdb_list                                  \n\
2 : Pocket characterization on a MD trajectory -> list of pre-aligned pdb ordered by \n\
    time file and wanted pocket pdb file \n\
\t./bin/mdpocket -L pdb_list -f wanted_pocket.pdb \n\
\t an example of a wanted pocket file can be obtained by running (1) \n\
\t (mdpout_iso_8.pdb) and non wanted grid points should be deleted by hand (i.e. PyMOL).\n\
\nOPTIONS (find standard parameters in brackets)           \n\n\
\t-o (char *) : common prefix of output file (mdpout_snapshots) \n\n\
\t-S : if you put this flag, the pocket score is matched to the \n\
\t density grid \n\
\t-a : output all bfactor colored snapshots\n\n\
\t-m (float)  : Minimum radius of an alpha-sphere.      (3.0)\n\
\t-M (float)  : Maximum radius of an alpha-sphere.      (6.0)\n\
\t-A (int)    : Minimum number of apolar neighbor for        \n\
\t              an a-sphere to be considered as apolar.   (3)\n\
\t-i (int)    : Minimum number of a-sphere per pocket.   (30)\n\
\t-D (float)  : Maximum distance for first clustering        \n\
\t              algorithm.                             (1.73)\n\
\t-s (float)  : Maximum distance for single linkage          \n\
\t              clustering                              (2.5)\n\
\t-n (integer): Minimum number of neighbor close from        \n\
\t              each other (not merged otherwise).        (3)\n\
\t-r (float)  : Maximum distance between two pockets         \n\
\t              barycenter (merged otherwise).          (4.5)\n\
\t-p (float)  : Minimum proportion of apolar sphere in       \n\
\t              a pocket (remove otherwise)             (0.0)\n\
\t-v (integer): Number of Monte-Carlo iteration for the      \n\
\t              calculation of each pocket volume.     (2500)\n\
\t-b (integer): Space approximation for the basic method     \n\
\t              of the volume calculation. Not used by       \n\
\t              default (Monte Carlo approximation is)       \n\
\nSee the manual (man fpocket), or the full documentation for\n\
more information.\n\
***************************\n" /**< the usage print content*/

/* --------------------------- PUBLIC STRUCTURES -----------------------------*/

/**
 Structure containing input and output params for mdpocket
 **/
typedef struct s_mdparams
{
	char **fsnapshot;       /**< path of the snapshot form of the structure */
        char fwantedpocket[M_MAX_PDB_NAME_LEN];    /**< path of the wanted pocket file*/

	char *f_pqr,        /**< name of the pqr concatenated snapshot file*/
		*f_densdx,      /**< name of the density dx grid file*/
                *f_freqdx,      /**< name of the frequency dx grid file*/
                *f_densiso,     /**< name of the density iso pdb file*/
                *f_freqiso,     /**< name of the frequency iso pdb file*/
                *f_desc,    /**< name of the descriptor file*/
                *f_ppdb,    /**< name of the pocket pdb output file */
                *f_apdb,    /**< name of the atoms pdb output file */
                *f_appdb;    /**< name of the all atoms pdb output file */
	int nfiles;         /**< number of files to analyse*/
	s_fparams *fpar ;   /**< fparams container*/
        int flag_scoring;   /**< perform fpocket scoring on grid points instead of voronoi vertice counting*/
        int bfact_on_all;   /**< flag, if 1, output all snapshots with surface coloured by bfactors*/

} s_mdparams ;

/* ----------------------------- PROTOTYPES --------------------------------- */

s_mdparams* init_def_mdparams(void) ;
s_mdparams* get_mdpocket_args(int nargs, char **args) ;

int add_list_snapshots(char *str_list_file, s_mdparams *par) ;
int add_snapshot(char *snapbuf, s_mdparams *par) ;

void print_mdparams(s_mdparams *p, FILE *f) ;
void print_mdpocket_usage(FILE *f) ;
void free_mdparams(s_mdparams *p) ;

#endif
