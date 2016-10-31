#include "../headers/mdpout.h"

/*

## GENERAL INFORMATION
##
## FILE 					mdpout.c.
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			28-11-10
##
## SPECIFICATIONS
##
##	This file implements a memory handler. Whenever you call
##	my_bloc_malloc function, the allocated pointer will be
##	stored in a simple chained list. Then, if a malloc fails,
##	or if you call the my_exit function, all allocated variable
##	will be freed if not NULL of course. Therefore, a bloc
##	allocated by my_... functions MUST be freed by the my_free,
##	or a double free error should appears.
##
##	WARNING
##
##	This is an easy way to handle memory. However, if many bloc
##	are allocated and removed during the programm, it might be
##	solwed down as for each free, as we have to look for the bloc
##	to free in the list, and remove it.
##
##	REMEMBER: if you use my_malloc, you MUST use my_free. If you
##	free a bloc allocated by my_(..), and if you call free_all at
##	the end, a double free coprruption will occure.
##
## MODIFICATIONS HISTORY
##
##	01-01-08	(vp) Created (random date...)
##
## TODO or SUGGESTIONS
##

*/


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



/**
   ## FUNCTION:
        write_md_grid

   ## SPECIFICATION:
	Write the md grid to a file.

   ## PARAMETRES:
	@ s_mdgrid g : structure containing the grid
        @ FILE *f : file handle for output file
        @ FILE *fiso : file handle for iso pdb output file
        @ s_mdparams *par : parameters for mdpocket
        @ float isovalue : isovalue at which one wants to extract the iso PDB

   ## RETURN:
	void :

*/


void write_md_grid(s_mdgrid *g, FILE *f, FILE *fiso,s_mdparams *par,float isovalue)
{
    int cx,cy,cz;
    float cv;
    float rx,ry,rz;
    size_t cnt=0;
    /*write the header of the dx file*/
    fprintf(f,"# Data calculated by mdpocket, part of the fpocket package\n");
    fprintf(f,"# This is a standard DX file of occurences of cavities within MD trajectories.\n");
    fprintf(f,"# The file can be visualised using the freely available VMD software\n");
    fprintf(f,"# fpocket parameters used to create this dx file : \n");
    fprintf(f,"# \t-m %2.f (min alpha sphere size) -M %.2f (max alpha sphere size)\n",par->fpar->asph_min_size, par->fpar->asph_max_size);
    fprintf(f,"# \t-i %d (min number of alpha spheres per pocket)\n",par->fpar->min_pock_nb_asph);
    fprintf(f,"# \t-D %.2f (Max distance for 1st clustering algo)\n#\t-r %.2f (Max dist for 2nd clustering algo)\n#\t-s %.2f (Max dist for third clustering algo)\n",par->fpar->clust_max_dist, par->fpar->refine_clust_dist, par->fpar->sl_clust_max_dist);
    fprintf(f,"# \t-n %d (Min neighbour atoms for multi linkage clustering)\n",par->fpar->sl_clust_min_nneigh);
    if(par->flag_scoring) fprintf(f,"# \t-S (Map drug score to density map!)\n");
    fprintf(f,"object 1 class gridpositions counts %d %d %d\n",g->nx,g->ny,g->nz);
    fprintf(f,"origin %.2f %.2f %.2f\n",g->origin[0],g->origin[1],g->origin[2]);
    fprintf(f,"delta %.2f 0 0\n",g->resolution);
    fprintf(f,"delta 0 %.2f 0\n",g->resolution);
    fprintf(f,"delta 0 0 %.2f\n",g->resolution);
    fprintf(f,"object 2 class gridconnections counts %d %d %d\n",g->nx,g->ny,g->nz);
    fprintf(f,"object 3 class array type double rank 0 items %d data follows\n",g->nx*g->ny*g->nz);
    int i=0;
    for(cx=0;cx<g->nx;cx++){
        for(cy=0;cy<g->ny;cy++){
            for(cz=0;cz<g->nz;cz++){
                if(i==3) {
                    i=0;
                    fprintf(f,"\n");
                }
                cv=g->gridvalues[cx][cy][cz];
                fprintf(f,"%.3f ",cv);
                if(cv>=isovalue){
                    cnt++;
                    rx=g->origin[0]+cx*g->resolution;
                    ry=g->origin[1]+cy*g->resolution;
                    rz=g->origin[2]+cz*g->resolution;
                    fprintf(fiso,"ATOM  %5d  C   PTH     1    %8.3f%8.3f%8.3f%6.2f%6.2f\n",(int)cnt,rx,ry,rz,0.0,0.0);
                }
                i++;
            }
        }
    }
}



/**
   ## FUNCTION:
        write_md_pocket_atoms

   ## SPECIFICATION:
        writes a pdb file containing all pocket atoms at each snapshot in
        distinct models

   ## PARAMETRES:
	@ FILE *f : output file handle
        @ int *ids : list of atom identifiers
        @ s_pdb *prot : the protein handle
        @ int nids : number of ids in ids list
        @ int sn : number of snapshots

   ## RETURN:
	void

*/
void write_md_pocket_atoms(FILE *f,int *ids,s_pdb *prot, int nids, int sn){
    s_atm *cura;
    int i,j,flag;
    i=0;
    j=0;
    fprintf(f,"MODEL        %d\n",sn);
    for(i=0;i<nids;i++){
        flag=0;
        while(flag==0 && j<prot->natoms){
            cura=prot->latoms_p[j];
            if(cura->id==ids[i]){
                flag=1;
                write_pdb_atom_line(f, "ATOM", cura->id, cura->name,
						 cura->pdb_aloc, cura->res_name, cura->chain,
						 cura->res_id, cura->pdb_insert, cura->x, cura->y, cura->z, cura->occupancy,
						 cura->bfactor,	cura->symbol, cura->charge);
            }
            j++;
        }
    }
    fprintf(f,"ENDMDL\n");
    
}



/**
   ## FUNCTION:
        write_first_bfactor_density

   ## SPECIFICATION:
	writes the protein structure with pocket densities in the bfactor column.

   ## PARAMETRES:
	@ FILE *f : output file handle
        @ s_pdb *prot : protein handle

   ## RETURN:
	void

*/
void write_first_bfactor_density(FILE *f,s_pdb *prot){
    s_atm *cura;
    int i;
    //fprintf(f,"MODEL        %d\n",sn);
    for(i=0;i<prot->natoms;i++){
        cura=prot->latoms_p[i];
        write_pdb_atom_line(f, "ATOM", cura->id, cura->name,
						 cura->pdb_aloc, cura->res_name, cura->chain,
						 cura->res_id, cura->pdb_insert, cura->x, cura->y, cura->z, cura->occupancy,
						 cura->bfactor,	cura->symbol, cura->charge);
    }
    fprintf(f,"TER\n");
    fprintf(f,"END\n");

}
