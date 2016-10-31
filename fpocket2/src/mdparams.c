
#include "../headers/mdparams.h"

/**

## GENERAL INFORMATION
##
## FILE 					mdparams.c
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			28-11-08
##
## SPECIFICATIONS
##
##	Handle parameters (parse the command line and sore values)
##	for the mdpocket programm.
##
## MODIFICATIONS HISTORY
##
##	28-11-08	(v)  Comments UTD + relooking
##	27-11-08	(v)  Minor Relooking
##	01-04-08	(v)  Added comments and creation of history
##	01-01-08	(vp) Created (random date...)
##
## TODO or SUGGESTIONS
##
##	(v) Check and update if necessary comments of each function!!
##	(v) Review the main function and handle all possible crashes.
##

*/

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

/**
   ## FUNCTION:
	init_def_mdparams

   ## SPECIFICATION:
	Initialisation of default parameters.

   ## PARAMETRES: void

   ## RETURN:
	s_mdparams*: Pointer to allocated paramers.

*/
s_mdparams *init_def_mdparams(void)
{
	s_mdparams *par = (s_mdparams*) my_malloc(sizeof(s_mdparams)) ;

	par->f_pqr = (char *)my_malloc(M_MAX_FILE_NAME_LENGTH*sizeof(char)) ;
	par->f_freqdx = (char *)my_malloc(M_MAX_FILE_NAME_LENGTH*sizeof(char)) ;
        par->f_densdx = (char *)my_malloc(M_MAX_FILE_NAME_LENGTH*sizeof(char)) ;
        par->f_freqiso = (char *)my_malloc(M_MAX_FILE_NAME_LENGTH*sizeof(char)) ;
        par->f_densiso = (char *)my_malloc(M_MAX_FILE_NAME_LENGTH*sizeof(char)) ;
        par->f_desc = (char *)my_malloc(M_MAX_FILE_NAME_LENGTH*sizeof(char)) ;
        par->f_ppdb = (char *)my_malloc(M_MAX_FILE_NAME_LENGTH*sizeof(char)) ;
        par->f_apdb = (char *)my_malloc(M_MAX_FILE_NAME_LENGTH*sizeof(char)) ;
        par->f_appdb = (char *)my_malloc(M_MAX_FILE_NAME_LENGTH*sizeof(char));
	strcpy(par->f_pqr, M_MDP_OUTPUT_FILE1_DEFAULT) ;
	strcpy(par->f_freqdx, M_MDP_OUTPUT_FILE2_DEFAULT) ;
        strcpy(par->f_freqiso, M_MDP_OUTPUT_FILE3_DEFAULT) ;
        strcpy(par->f_desc, M_MDP_OUTPUT_FILE4_DEFAULT) ;
        strcpy(par->f_ppdb, M_MDP_OUTPUT_FILE5_DEFAULT) ;
        strcpy(par->f_apdb, M_MDP_OUTPUT_FILE6_DEFAULT) ;
        strcpy(par->f_appdb, M_MDP_OUTPUT_FILE7_DEFAULT);
        strcpy(par->f_densdx, M_MDP_OUTPUT_FILE8_DEFAULT) ;
        strcpy(par->f_densiso, M_MDP_OUTPUT_FILE9_DEFAULT) ;
	par->fsnapshot = NULL ;
        par->fwantedpocket[0] = 0 ;
	par->nfiles = 0 ;
        par->flag_scoring=0;
        par->bfact_on_all=0;

	return par ;
}

/**
   ## FUNCTION:
	get_dpocket_args

   ## SPECIFICATION:
	This function analyse the user's command line and parse it to store
	parameters for the descriptor calculator programm.

   ## PARAMETRES:
	@ int nargs   : Number of arguments
	@ char **args : Arguments of main program

   ## RETURN:
	s_mdparams*: Pointer to parameters

*/
s_mdparams* get_mdpocket_args(int nargs, char **args)
{
	int i,
		status = 0,
		nstats = 0,
                npdb=0;

	char *str_list_file = NULL ;

	s_mdparams *par = init_def_mdparams() ;
	par->fpar = get_fpocket_args(nargs, args) ;

	/* Read arguments by flags */
	for (i = 1; i < nargs; i++) {
		if (strlen(args[i]) == 2 && args[i][0] == '-') {
			switch (args[i][1]) {
				case M_MDPAR_OUTPUT_FILE :
						if(nstats >= 1) fprintf(stdout, "! More than one single file for the stats output file has been given. Ignoring this one.\n") ;
						else {
							if(i < nargs-1) {
								if(strlen(args[++i]) < M_MAX_FILE_NAME_LENGTH) {
									remove_ext(args[i]) ;
									sprintf(par->f_pqr, "%s.pqr", args[i]) ;
									sprintf(par->f_freqdx, "%s_freq.dx", args[i]) ;
                                                                        sprintf(par->f_densdx, "%s_dens.dx", args[i]) ;
                                                                        sprintf(par->f_freqiso, "%s_freq_iso_0_5.pdb", args[i]) ;
                                                                        sprintf(par->f_densiso, "%s_dens_iso_8.pdb", args[i]) ;
                                                                        sprintf(par->f_desc, "%s_descriptors.txt", args[i]) ;
                                                                        sprintf(par->f_ppdb, "%s_mdpocket.pdb", args[i]) ;
                                                                        sprintf(par->f_apdb, "%s_mdpocket_atoms.pdb", args[i]) ;
                                                                        sprintf(par->f_appdb, "%s_all_atom_pdensities.pdb", args[i]) ;
								}
								else fprintf(stdout, "! Output file name is too long... Keeping default.") ;
							}
							else {
								fprintf(stdout, "! Invalid output file name argument missing.\n") ;
								status += 1 ;
							}
						}
						break ;

				case M_MDPAR_INPUT_FILE :
						if(i < nargs-1) str_list_file = args[++i] ;
						else {
							fprintf(stdout, "! Input file name argument missing.\n") ;
							status += 1 ;
						}
						 break ;
                                case M_MDPAR_INPUT_FILE2 :
						if(npdb >= 1) fprintf(stderr,
							"! Only first input pdb will be used.\n") ;
						else {
							strcpy(par->fwantedpocket, args[++i]) ; npdb++ ;
						}
						break ;
                                case M_MDPAR_SCORING_MODE :
                                                par->flag_scoring=1;
                                                break;
                                case M_MDPAR_OUTPUT_ALL_SNAPSHOTS :
                                                par->bfact_on_all=1;
                                                break;
				default:
					//  Check fpocket parameters!
					if(!is_fpocket_opt(args[i][1])) {
						fprintf(stdout, "> Unknown option '%s'. Ignoring it.\n",
								args[i]) ;
					}
					break ;
			}
		}
	}

	if(status > 0) {
		free_mdparams(par) ;
		par = NULL ;
 		print_mdpocket_usage(stdout);
	}
	else {
		if(str_list_file) {
			int res = add_list_snapshots(str_list_file, par) ;
			if(res <= 0) {
				fprintf(stdout, "! No data has been read.\n") ;
				free_mdparams(par) ;
				par = NULL ;
 				print_mdpocket_usage(stdout);
			}
		}
		else {
			fprintf(stdout, "! No input file given... Try again :).\n") ;
			free_mdparams(par) ;
			par = NULL ;
	 		print_mdpocket_usage(stdout);
		}
	}

	return par;
}

/**
   ## FUNCTION:
	add_list_snapshots

   ## SPECIFICATION:
	Load a list of snapshot pdb file path. This file should have the
	following format:

	snapshot_pdb_file
	snapshot_pdb_file2
	snapshot_pdb_file3
	(...)

 	Each snapshot file will be stored in the parameters structure.

   ## PARAMETRES:
	@ char *str_list_file : Path of the file containing all data
	@ s_mdparams *par      : Structures that stores all thoses files

   ## RETURN:
	int: Number of file read.

*/
int add_list_snapshots(char *str_list_file, s_mdparams *par)
{
	FILE *f;
	int n,
		nread = 0,
		status ;

	char buf[M_MAX_PDB_NAME_LEN*2 + 6],
		 snapbuf[M_MAX_PDB_NAME_LEN];

	/* Loading data. */
	f = fopen(str_list_file, "r") ;
/*
	printf(str_list_file);
*/
	if(f) {
		while(fgets(buf, 210, f)) {
/*
			printf("B: %s\n" , buf);
*/
			n = par->nfiles ;
			status = sscanf(buf, "%s", snapbuf) ;
			if(status < 1) {

				fprintf(stderr, "! Skipping row '%s' with bad format (status %d).\n",
								buf, status) ;
			}
			else {
				nread += add_snapshot(snapbuf, par) ;
			}

		}
	}
	else {
		fprintf(stderr, "! File %s doesn't exists\n", str_list_file) ;
	}
        fclose(f);
        return nread ;
}

/**
   ## FUNCTION:
	add_snapshot

   ## SPECIFICATION:
	Add a set of data to the list of set of data in the parameters. this function
	is used for the mdpocket program only.

	The function will try to open the file, and data will be stored only if the
	file exists, and if the name of the ligand is valid.

   ## PARAMETERS:
	@ char *snapbuf     : The snapshots path
	@ s_mdparams *par: The structure than contains parameters.

   ## RETURN:
	int: 1 if everything is OK, 0 if not.

*/
int add_snapshot(char *snapbuf, s_mdparams *par)
{
	int nm1 ;

	FILE *f = fopen_pdb_check_case(snapbuf, "r") ;
	if(f) {
                nm1 = par->nfiles ;
                par->nfiles += 1 ;

                
                par->fsnapshot = (char**) my_realloc(par->fsnapshot, (par->nfiles)*sizeof(char*)) ;

                par->fsnapshot[nm1] = (char *)my_malloc((strlen(snapbuf)+1)*sizeof(char)) ;

                strcpy(par->fsnapshot[nm1], snapbuf) ;

                fclose(f) ;

	}
	else {
		fprintf(stdout, "! The pdb file '%s' doesn't exists.\n", snapbuf) ;
		return 0 ;
	}

	return 1 ;
}

/**
   ## FUNCTION:
	print_mdparams

   ## SPECIFICATION:
	Print function, usefull to debug

   ## PARAMETRES:
	@ s_mdparams *p : The structure than will contain the parsed parameter
	@ FILE *f      : The file to write in

   ## RETURN:
    void

*/
void print_mdparams(s_mdparams *p, FILE *f)
{

	if(p) {
		fprintf(f, "==============\nParameters of the program: \n");
		int i ;
		for(i = 0 ; i < p->nfiles ; i++) {
			fprintf(f, "> Snaphot %d: '%s'\n", i+1, p->fsnapshot[i]) ;
		}
		fprintf(f, "==============\n");
                if(p->fwantedpocket[0]!=0){
                    fprintf(f,"Wanted pocket given in file : %s\n",p->fwantedpocket);
                    fprintf(f, "==============\n");
                }
	}
	else fprintf(f, "> No parameters detected\n");
}

/**
   ## FUNCTION:
	print_mdparams_usage

   ## SPECIFICATION:
	Displaying usage of the programm in the given buffer

   ## PARAMETRES:
	@ FILE *f: buffer to print in

   ## RETURN:

*/
void print_mdpocket_usage(FILE *f)
{
	f = (f == NULL) ? stdout:f ;

	fprintf(f, M_MDP_USAGE) ;
}

/**
   ## FUNCTION:
	free_params

   ## SPECIFICATION:
	Free parameters

   ## PARAMETRES:
	@ s_mdparams *p: Pointer to the structure to free

   ## RETURN:
	void

*/
void free_mdparams(s_mdparams *p)
{
	if(p) {
		
		if(p->fsnapshot) {
			my_free(p->fsnapshot);
			p->fsnapshot = NULL;
		}
		if(p->f_pqr) {
			my_free(p->f_pqr);
			p->f_pqr = NULL;
		}
		if(p->f_densdx) {
			my_free(p->f_densdx);
			p->f_densdx = NULL ;
		}
                if(p->f_freqdx) {
			my_free(p->f_freqdx);
			p->f_freqdx = NULL ;
		}
		if(p->f_desc) {
			my_free(p->f_desc);
			p->f_desc = NULL;
		}
		free_fparams(p->fpar);
 		my_free(p) ;
	}
}
