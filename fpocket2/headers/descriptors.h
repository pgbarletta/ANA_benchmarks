
/*
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

#ifndef DH_DESCR
#define DH_DESCR

/* ---------------------------------INCLUDES--------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <string.h>

#include "voronoi.h"
#include "voronoi_lst.h"
#include "atom.h"
#include "aa.h"
#include "utils.h"


/* --------------------------------STRUCTURES---------------------------------*/
/**
        Structure containing descriptors for binding pockets
*/
typedef struct s_desc 
{
	float hydrophobicity_score, /**< Hydropathie score - for each aa */
              volume_score,         /**< Volume score - for each aa */
              volume,               /**< Pocket volume */
              prop_polar_atm,       /**< Proportion of polar atoms */
              mean_asph_ray,        /**< Mean alpha sphere radius */
              masph_sacc,           /**< Mean alpha sphere solvent accessibility */
              apolar_asphere_prop,  /**< Proportion of apolar alpha spheres */
              mean_loc_hyd_dens,    /**< Mean local hydrophobic density (from alpha spheres) */
              as_density,           /**< Pocket density, defined as mean distance between alpha spheres*/
              as_max_dst,           /**< Maximum distance between two alpha spheres */
              /**< The following descriptors are all normalized using observed
                 values among all pocket found by the algorithm. These
                 are not set in descriptor.c, but in pocket.c as we have to check
                 all pocket first to store boundaries of the descriptor to
                 normalize. */

              flex,                  /**< Normalized flexibility - based on B factors - ABUSIVE */
              nas_norm,              /**< Normalized number of alpha sphere */
              polarity_score_norm,   /**< Normalized polarity score*/
              mean_loc_hyd_dens_norm,/**< Normalized mean local hydrophobic density */
              prop_asapol_norm,      /**< Normalized proportion of apolar alphasphere */
              as_density_norm,       /**< Normalized alpha sphere density */
              as_max_dst_norm,       /**< normalized maximum distance between alpha sphere centers*/
                
                /**< The following descriptors are various surface calculations*/
              surf_vdw,              /**< Van der Waals surface of the pocket*/
              surf_vdw14,            /**< Van der Waals surface + 1.4 A probe*/
              surf_vdw22,            /**< Van der Waals surface + 2.2 A probe*/
              surf_pol_vdw,          /**< polar van der Waals surface of the pocket*/
              surf_pol_vdw14,        /**< polar van der Waals surface + 1.4 A probe*/
              surf_pol_vdw22,        /**< polar an der Waals surface + 2.2 A probe*/
              surf_apol_vdw,         /**< polar van der Waals surface of the pocket*/
              surf_apol_vdw14,       /**< polar van der Waals surface + 1.4 A probe*/
              surf_apol_vdw22        /**< polar van der Waals surface + 2.2 A probe*/
        ;
	
	int aa_compo[20] ;	/**< Absolute amino acid composition */
	int nb_asph,		/**< Number of alpha spheres */
            polarity_score,	/**< Polarity score (based on amino acids properties ; see aa.c & aa.h) */
            charge_score ;	/**< Sum of all net charges at pH = 7 (see aa.c & aa.h) */
        float as_max_r ;        /**< Alpha sphere maximum radius*/
        float drug_score;       /**< Drug score of the binding site*/
        int interChain,                             /**< 0 if pocket in single chain, 1 if between 2 chains*/
                characterChain1,                         /**< 0 if protein, 1 if nucl acid pocket, 2 if HETATM pocket*/
                characterChain2,                         /**< 0 if protein, 1 if nucl acid pocket, 2 if HETATM pocket*/
                numResChain1,                            /**<number of resdiues on chain 1*/
                numResChain2;                           /**<number of res on chain 2*/
        char nameChain1[2],        /**<name of the first chain in contact with the pocket*/
                nameChain2[2];     /**<name of the second chain in contact with the pocket, if there is*/
        char ligTag[8];           /**<het atom tag of ligands situated in the pocket*/
        int n_abpa;               /**<number of abpas in the binding site*/
} s_desc ;

/* ------------------------------PROTOTYPES---------------------------------- */

s_desc* allocate_s_desc(void) ;
void reset_desc(s_desc *desc) ;


void set_descriptors(s_atm **tatoms, int natoms, s_vvertice **tvert, int nvert, s_desc *desc, int niter,s_pdb *pdb,int flag_do_expensive_calculations) ;

int get_vert_apolar_density(s_vvertice **tvert, int nvert, s_vvertice *vert) ;
void set_atom_based_descriptors(s_atm **atoms, int natoms, s_desc *desc,s_atm *all_atoms, int all_natoms);
void set_aa_desc(s_desc *desc, const char *aa_name) ;
int countResidues(s_atm *atoms, int natoms, char chain[2]);

#endif
