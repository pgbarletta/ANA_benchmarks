
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

#ifndef DH_MDPMEM
#define DH_MDPMEM

/* ----------------------------- INCLUDES ------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>

#include "fpocket.h"
#include "fpout.h"
#include "writepocket.h"
#include "tpocket.h"
#include "dparams.h"
#include "descriptors.h"
#include "neighbor.h"
#include "pocket.h"
#include "cluster.h"
#include "refine.h"
#include "aa.h"
#include "utils.h"
#include "mdparams.h"
#include "mdpbase.h"
#include "memhandler.h"
/* ---------------------------------MACROS----------------------------------*/

/* -------------------------------PROTOTYPES--------------------------------*/


s_mdconcat *init_md_concat(void);
void alloc_first_md_concat(s_mdconcat *m,size_t n);
void realloc_md_concat(s_mdconcat *m,size_t n);
void free_mdconcat(s_mdconcat *m);

#endif
