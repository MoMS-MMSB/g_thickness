#include <stdlib.h>
#include <stdio.h>

#include <ctype.h> // Needed for toupper

#include <gromacs/copyrite.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/pbc.h>
#include <gromacs/rmpbc.h>
#include <gromacs/smalloc.h>
#include <gromacs/statutil.h>
#include <gromacs/vec.h>
#include <gromacs/xvgr.h>

#include "grid_mode.h"
#include "dist_mode.h"

static const char *authors[] = {
    "Written by Jonathan Barnoud (jonathan.barnoud@inserm.fr)",
    "Copyright (c) 2012  Jonathan Barnoud, Luca Monticelli"
};
static const char *gpl[] = {
    "This program is free software; you can redistribute it and/or",
    "modify it under the terms of the GNU General Public License",
    "as published by the Free Software Foundation; either version 2",
    "of the License, or (at your option) any later version."
};
static const char *version[] = {
    ":-) g_thickness - version 1.0 (-:"
};
void sp_print(FILE *out, const char *txt)
{
    int i, s;
    s = (int) (80 - strlen(txt)) / 2.0;
    for (i=0; i < s; i++) 
        fprintf(out, " ");
    fprintf(out, "%s\n", txt);
}


typedef struct GeneralData {
    int ngrps;
    atom_id **index;
    int *isize;
    char **grpnames;
    const char *traj_fn;
    int adt;
} GeneralData;

typedef struct t_modes {
    GridHeight *grid_store;
    DistMode *dist_store;
    GeneralData *general;
} t_modes;

/*****************************************************************************
 *                               I/O stuff                                   *
 *****************************************************************************/
void clean_modes(t_modes *modes) {
    int group;
    clean_grids(modes->grid_store);
    clean_dist(modes->dist_store);
    for (group = 0; group < modes->general->ngrps; ++group) {
        sfree(modes->general->index[group]);
        sfree(modes->general->grpnames[group]);
    }
    sfree(modes->general->index);
    sfree(modes->general->isize);
    sfree(modes->general->grpnames);
}

/** Read user choices and prepare the run
 */
t_modes handle_user(int argc, char **argv,
        output_env_t *oenv, t_topology **top, int *ePBC) {
    int i;
    /* Output variable */
    t_modes modes;
    /* Variables for the reading of the arguments */
    static const char *axtitle[] = { NULL, "z", "x", "y", NULL };
    static const char *prof_axtitle[] = { NULL, "d", "z", "x", "y", NULL };
    int axis = 0;
    int axis_prof = 0;
    int sl = 100;
    int sl2 = -1;
    int adt = -1;
    gmx_bool bGrid = TRUE;
    gmx_bool bDist = TRUE;
    gmx_bool bCOM = TRUE;
    /* Variables for the reading of the common index file */
    atom_id **index = NULL;
    int *isize = NULL;
    char **grpnames = NULL;
    static const int ngrps = 2;
    
    const char *desc[] = {
        "Calculate the local thickness of a membrane.",
        "[PAR]",
        "The program can calculate the thickness landscape of the membrane",
        "using the [TT]-og[tt] option and the thickness profile as a function",
        "to the distance of a group using the [TT]-od[tt] option. At least",
        "one of these two options have to be used.",
        "[PAR]",
        "Thickness is calculated as the distance between its two leaflets",
        "along the normal axis. This axis have to be a unit axis; it can be",
        "chosen using the [TT]-d[tt] option.",
        "[PAR]",
        "When calculating a landscape, [TT]-sl[tt] and [TT]-sl2[tt] correspond",
        "to the number of cells in each dimension. If [TT]-sl2[tt] is negative",
        "then it takes the value of [TT]-sl[tt]. When calculating a landscape,",
        "[TT]-sl[tt] corresponds to the number of bins; [TT]-sl2[tt] is",
        "ignored.",
        "[PAR]",
        "The distance to a reference group is calculated, by default, as the",
        "distance to the center of mass of the reference group. It can be",
        "calculated as the minimum distance using [TT]-nocom[tt].",
        "[PAR]",
        "The [TT]-adt[tt] option allows to chose how frequently the thickness",
        "is calculated. When this option is set to a value greater than one,",
        "the center of mass of each leaflet for a cell or a bin is averaged",
        "over several frames.",
        "[PAR]",
        "See the README for more details."
    };

    t_pargs pa[] = {
        { "-d",    FALSE, etENUM,  {axtitle}, "Membrane normal dimension."},
        { "-sl", FALSE, etINT, {&sl}, "Number of grid cells per side."},
        { "-sl2", FALSE, etINT, {&sl2}, "Number of grid cells on the second "
            "dimension. If lesser or equal 0 the value of -sl is used."},
        { "-adt", FALSE, etINT, {&adt},
            "Thickness will be averaged when nsteps \% adt will be null or at "
                "the end if adt is lesser than 0 or bigger than the simulation "
                "length."},
        { "-com", FALSE, etBOOL, {&bCOM},
            "If true center of mass distance, else use minimum distance."},
    };
    #define NPA asize(pa)
    t_filenm fnm[] = {
        { efTPX, "-s", NULL, ffREAD},  /* this is for the topology   */
        { efTRX, "-f", NULL, ffREAD},  /* this is for the trajectory */
        { efNDX, "-n", NULL, ffREAD},  /* this is for the index file */
        /* output for the grid mode data and sampling */
        { efDAT, "-og", "thickness_grid", ffOPTWR }, 
        { efDAT, "-ogs", "thickness_grid_sampling", ffOPTWR }, 
        /* output for the dist mode data and sampling */
        { efXVG, "-od", "thickness_dist", ffOPTWR }, 
        { efXVG, "-ods", "thickness_dist_sampling", ffOPTWR }, 
    };
    #define NFILE asize(fnm)

    /* Display g_thickness credit */
    for (i=0; i < (int)asize(version); i++)
        sp_print(stderr, version[i]);
    fprintf(stderr, "\n");
    for (i=0; i < (int)asize(authors); i++)
        sp_print(stderr, authors[i]);
    fprintf(stderr, "\n");
    for (i=0; i < (int)asize(gpl); i++)
        sp_print(stderr, gpl[i]);
    fprintf(stderr, "\n");
    /* Display GROMACS credit */
    CopyRight(stderr,argv[0]);
    /* Parse the command line arguments */
    parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,
	    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL,oenv);
	bGrid = opt2bSet("-og",NFILE,fnm);
	bDist = opt2bSet("-od",NFILE,fnm);

	if (! (bDist || bGrid)) {
	    gmx_fatal(FARGS, "You need to choose at least one output"
	                     "(see -og and -od options)");
	}

	/* Convert axis in int */
    axis = toupper(axtitle[0][0]) - 'X';

    /* Look at -sl2 */
    if (sl2 <= 0) {
        sl2 = sl;
    }

	/* Read topology */
    (*top)=read_top(ftp2fn(efTPX,NFILE,fnm),ePBC);

	/* Read index */
    snew(index, ngrps);
    snew(isize, ngrps);
    snew(grpnames, ngrps);
    printf("Select groups for the leaflets:\n");
    get_index(&((*top)->atoms),ftp2fn(efNDX,NFILE,fnm),ngrps,
            isize,index,grpnames);

	/* Create the mode objects */
    snew(modes.general, 1);
	modes.general->ngrps = ngrps;
	modes.general->index = index;
	modes.general->isize = isize;
	modes.general->grpnames = grpnames;
	modes.general->traj_fn = ftp2fn(efTRX,NFILE,fnm);
	modes.general->adt = adt;

	modes.grid_store = NULL;
	modes.dist_store = NULL;
	if (bGrid) {
	    modes.grid_store = build_grids((int [2]){sl, sl2}, axis,
	            opt2fn("-og",NFILE,fnm), opt2fn("-ogs",NFILE,fnm));
	}
	if (bDist) {
        modes.dist_store = build_dist(sl, axis,
                opt2fn("-od",NFILE,fnm), opt2fn("-ods",NFILE,fnm), *oenv,
                ftp2fn(efNDX,NFILE,fnm), *top, bCOM);
	}
	
	return modes;
}

/*****************************************************************************
 *                            Trajectory reading                             *
 *****************************************************************************/
void do_frame(t_modes modes, t_pbc *pbc, int ePBC, matrix box, rvec *x,
        gmx_rmpbc_t gpbc, int natoms, t_topology *top) {
    int leaflet = 0;
    int atom = 0;
    if (pbc) {
        set_pbc(pbc,ePBC,box);
        /* make molecules whole again */
        gmx_rmpbc(gpbc,natoms,box,x);
    }
    grid_start_frame(modes.grid_store, box);
    dist_start_frame(modes.dist_store, box, top, x);
    for (leaflet = 0; leaflet < modes.general->ngrps; ++leaflet) {
        for (atom = 0; atom < modes.general->isize[leaflet]; ++atom) {
            grid_store(modes.grid_store, leaflet, 
                    x[modes.general->index[leaflet][atom]], pbc);
            dist_store(modes.dist_store, leaflet,
                    modes.general->index[leaflet][atom], x, pbc);
        }
    }
    grid_end_frame(modes.grid_store, modes.general->adt);
    dist_end_frame(modes.dist_store, modes.general->adt);
}

void read_traj(t_modes modes, output_env_t oenv, t_topology *top, int ePBC) {
    int natoms;
    real time;
    rvec *x;
    matrix box;
    t_pbc *pbc;
    t_trxstatus *status;
    gmx_rmpbc_t gpbc=NULL;

    /* Read the first frame to get basic informations about the system */
    natoms=read_first_x(oenv,&status,modes.general->traj_fn,&time,&x,box);
    /* Set PBC stiff */
    if (ePBC != epbcNONE)
        snew(pbc,1);
    else
        pbc = NULL;
    gpbc = gmx_rmpbc_init(&top->idef,ePBC,natoms,box);
    /* Read the trajectory */
    do {
        do_frame(modes, pbc, ePBC, box, x, gpbc, natoms, top);
    } while(read_next_x(oenv,status,&time,natoms,x,box));
}

int main(int argc, char **argv) {
    output_env_t oenv;
    t_topology *top;
    int ePBC;
    t_modes modes;

    /* Read user input */
    modes = handle_user(argc, argv, &oenv, &top, &ePBC);
    /* Read the trajectory */
    read_traj(modes, oenv, top, ePBC);
    /* Write results */
    grid_end(modes.grid_store, modes.general->adt);
    dist_end(modes.dist_store, modes.general->adt);
    /* Clean everything */
    clean_modes(&modes);
    return 0;
}
