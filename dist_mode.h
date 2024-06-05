#ifndef _dist_mode_h
#define _dist_mode_h

#include <math.h>

#include <gromacs/statutil.h>
#include <gromacs/macros.h>
#include <gromacs/smalloc.h>
#include <gromacs/typedefs.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/pbc.h>
#include <gromacs/xvgr.h>
#include <gromacs/vec.h>

#include "distances.h"

typedef struct DistMode {
    real *height[3];    
    int  *sampling[3];
    int  length;
    FILE *out_dist;
    FILE *out_sampling;
    real width;
    int axis[2];
    real box_width;
    int nframes;
    atom_id *ref_index;
    int ref_size;
    real mass;
    gmx_bool bCOM;
    rvec *com;
} DistMode; 

DistMode *build_dist(int length, int normal_axis,
        const char *dist_fn, const char *sampling_fn, output_env_t oenv,
        const char *index_fn, t_topology *top, gmx_bool bCOM);

void clean_dist(DistMode *dist_store);

void clean_dist(DistMode *dist_store);

void dist_start_frame(DistMode *dist_store, matrix box, t_topology *top,
                      rvec *x);

void dist_end_frame(DistMode *dist_store, int adt);

void dist_store(DistMode *dist, int leaflet, int atom, rvec *x, t_pbc *pbc);

void dist_end(DistMode *dist_store, int adt);

#endif
