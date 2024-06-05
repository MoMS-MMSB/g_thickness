#include "distances.h"
#include <stdio.h>
#include <ctype.h>
#include <gromacs/sysstuff.h>
#include <gromacs/futil.h>
#include <gromacs/string2.h>
#include <gromacs/macros.h>
#include <gromacs/smalloc.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/statutil.h>
#include <gromacs/gmxfio.h>
#include <gromacs/maths.h>
#include <gromacs/vec.h>
#include <gromacs/pbc.h>

void make_2D(rvec vector, int axis, rvec result) {
    int i;
    for (i=0; i<DIM; ++i) {
        if (i != axis) {
            result[i] = vector[i];
        }
        else {
            result[i] = 0;
        }
    }
}

void pbc_rvec_sub(const t_pbc *pbc,const rvec xi,const rvec xj,rvec dx) {
    if (pbc)
        pbc_dx(pbc,xi,xj,dx);
    else
        rvec_sub(xi, xj, dx);
}

real get_distance(rvec pointA, rvec pointB, t_pbc *pbc) {
    rvec dx;
    real dist;
    pbc_rvec_sub(pbc, pointA, pointB, dx);
    dist = sqrtf(norm2(dx));
    return dist;
}

real min_dist(rvec pointA, atom_id *group, int grp_size, rvec *x, t_pbc *pbc, int axis) {
    int i;
    real dist;
    real min_dist = GMX_REAL_MAX;
    rvec pointA_mod, pointB;
    make_2D(pointA, axis, pointA_mod);
    for (i=0; i<grp_size; ++i) {
        make_2D(x[group[i]], axis, pointB);
        dist = get_distance(pointA_mod, pointB, pbc);
        if (dist < min_dist) {
            min_dist = dist;
        }
    }
    return min_dist;
}

real get_mass(atom_id *group, int grp_size, t_topology *top) {
    real mass = 0;
    int i = 0;
    for (i=0; i<grp_size; ++i) {
        mass += top->atoms.atom[group[i]].m;
    }
    return mass;
}

rvec *center_of_mass(atom_id *group, int grp_size, rvec *x,
        t_topology *top, real mass) {
    rvec *com = NULL;
    int i = 0, dim=0;
    snew(com, 1);
    for (i=0; i<grp_size; ++i) {
        for (dim=0; dim<DIM; ++dim) {
            (*com)[dim] += x[group[i]][dim] * top->atoms.atom[group[i]].m;
        }
    }
    for (dim=0; dim<DIM; ++dim) {
        (*com)[dim] /= mass;
    }
    return com;
}

real dist_2D(rvec pointA, rvec pointB, t_pbc *pbc, int axis) {
    rvec pointA_mod = {0,0,0}, pointB_mod = {0,0,0};
    make_2D(pointA, axis, pointA_mod);
    make_2D(pointB, axis, pointB_mod);
    return get_distance(pointA_mod, pointB_mod, pbc);
}
