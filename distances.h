#ifndef _distances_h
#define _distances_h

#include <gromacs/typedefs.h>

void make_2D(rvec vector, int axis, rvec result);

/** Substract 2 vector of coordinates and take PBC into account if needed
 *
 * :Parameters:
 *     - pbc : the periodic box description
 *     - xi, xj : xj will be substracted to xi
 *     - dx  : destination vector
 **/
void pbc_rvec_sub(const t_pbc *pbc,const rvec xi,const rvec xj,rvec dx);


/** Get the distance between two points with or without periodic conditions
 */
real get_distance(rvec pointA, rvec pointB, t_pbc *pbc);

/** Get the minimum distance between a point and a group of points
 */
real min_dist(rvec pointA, atom_id *group, int grp_size, rvec *x, t_pbc *pbc,
        int axis);

/** Get mass of a group
 */
real get_mass(atom_id *group, int grp_size, t_topology *top);

/** Get the center of mass of a group of atoms
 */
rvec *center_of_mass(atom_id *group, int grp_size, rvec *x,
        t_topology *top, real mass);

/** Get the distance between two point in 2D
 */
real dist_2D(rvec pointA, rvec pointB, t_pbc *pbc, int axis);

#endif	/* _distances_h */
