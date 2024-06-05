#include "dist_mode.h"

void _average_field_dist(DistMode *dist_store) {
    int leaflet, i;
    for (leaflet = 0; leaflet < 2; ++leaflet) {
        for (i=0; i < dist_store->length; ++i) {
            dist_store->height[leaflet][i] /=
                dist_store->sampling[leaflet][i];
        }
    }
}

void _calculate_thickness_dist(DistMode *dist_store) {
    int i;
    int minsamp = 0;
    for (i=0; i < dist_store->length; ++i) {
        minsamp = min(dist_store->sampling[0][i],
                      dist_store->sampling[1][i]);
        if (minsamp > 0) {
            dist_store->height[2][i] +=
                (real)fabs(dist_store->height[0][i] - 
                        dist_store->height[1][i]) * minsamp;
            dist_store->sampling[2][i] += minsamp;
        }
    }
}

DistMode *build_dist(int length, int normal_axis,
        const char *dist_fn, const char *sampling_fn, output_env_t oenv,
        const char *index_fn, t_topology *top, gmx_bool bCOM) {
    DistMode *dist_store;
    int prof, i;
    atom_id **index;
    int *isize;
    char **grpnames;

    /* Check dimensions */
    if (length <= 0) {
        fprintf(stderr,
                "I can not build a distance profile with this length: %d\n",
                length);
        exit(1);
    }

    /* Allocate empty structure */
    snew(dist_store, 1);
    
    dist_store->nframes = 0;
    /* Store the shape */
    dist_store->length = length;
    /* Set the slice width to 0, just in case */
    dist_store->width = 0;
    /* Set box_width sommation to 0 */
    dist_store->box_width = 0.0;
    /* Define the axis */
    dist_store->axis[0] = normal_axis;
    dist_store->axis[1] = 0;

    /* Allocate the profiles */
    for (prof = 0; prof < 3; ++prof) {
        snew(dist_store->height[prof], length);
        snew(dist_store->sampling[prof], length);
        for (i=0; i<length; ++i) {
            dist_store->height[prof][i] = 0;
            dist_store->sampling[prof][i] = 0;
        }
    }

    /* Get the reference group index */
    snew(index, 1);
    snew(isize, 1);
    snew(grpnames, 1);
    printf("Select reference group for distance calcultation:\n");
    get_index(&(top->atoms), index_fn, 1, isize,index,grpnames);
    dist_store->ref_index = index[0];
    dist_store->ref_size = isize[0];

    /* Calculate the reference group mass if needed */
    dist_store->bCOM = bCOM;
    if (bCOM) {
        dist_store->mass = get_mass(index[0], isize[0], top);
    }
    else {
        dist_store->mass = 0;
    }

    /* Open the files */
    dist_store->out_dist = xvgropen(dist_fn,"Thickness",
            "Distance from Protein (nm)","z coordinate (nm)",oenv);
    dist_store->out_sampling = xvgropen(sampling_fn,"Sampling",
            "Distance from Protein (nm)","Average number of hit",oenv);
    return dist_store;
}

void clean_dist(DistMode *dist_store) {
    int prof = 0;
    if (dist_store) {
        for (prof = 0; prof < 3; ++prof) {
            sfree(dist_store->height[prof]);
            sfree(dist_store->sampling[prof]);
        }
        sfree(dist_store->ref_index);
        fclose(dist_store->out_dist);
        fclose(dist_store->out_sampling);
        sfree(dist_store);
    }
}

void dist_start_frame(DistMode *dist_store, matrix box, t_topology *top,
                      rvec *x) {
    int i = 0;
    real max_box_size = 0;
    if (dist_store) {
        /* Find what the maximum distance is */
        for (i=0; i<DIM; ++i) {
            if (i != dist_store->axis[0]) {
                max_box_size += box[i][i] * box[i][i];
            }
        }
        max_box_size = sqrt(max_box_size)/2;
        dist_store->nframes += 1;
        dist_store->width = max_box_size/dist_store->length;
        dist_store->box_width += max_box_size;
        /* Get reference group center of mass if needed */
        if (dist_store->bCOM) {
            dist_store->com = center_of_mass(dist_store->ref_index,
                    dist_store->ref_size, x, top, dist_store->mass);
        }
    }
}

void dist_end_frame(DistMode *dist_store, int adt) {
    int i, leaflet;
    if (dist_store && adt > 0 && dist_store->nframes % adt == 0) {
        /* Average the fields */
        _average_field_dist(dist_store);
        /* Calculate the thickness */
        _calculate_thickness_dist(dist_store);
        /* Empty the fields */
        for (leaflet = 0; leaflet < 2; ++leaflet) {
            for (i=0; i < dist_store->length; ++i) {
                dist_store->height[leaflet][i] = 0.0;
                dist_store->sampling[leaflet][i] = 0;
            }
        }
        /* Free the center of mass if needed */
        if (dist_store->bCOM) {
            sfree(dist_store->com);
        }
    }
}

void dist_store(DistMode *dist, int leaflet, int atom, rvec *x, t_pbc *pbc) {
    int slice = 0;
    int i = 0;
    real distance = 0;
    rvec *com = NULL;
    /*real distance2 = 0;*/
    if (dist) {
        if (dist->bCOM) {
            distance = dist_2D(x[atom], *(dist->com), pbc, dist->axis[0]);
        }
        else {
            distance = min_dist(x[atom], dist->ref_index, dist->ref_size,
                    x, pbc, dist->axis[0]);
        }

        slice = distance/dist->width;
        if (slice >= dist->length) {
            printf("arrggg: %d >= %d\n", slice, dist->length);
        }
        dist->height[leaflet][slice] += x[atom][dist->axis[0]];
        dist->sampling[leaflet][slice] += 1;
    }
}

void dist_end(DistMode *dist_store, int adt) {
    if (dist_store) {
        int i, leaflet;
        real bin_size = 0;
        dist_store->box_width /= dist_store->nframes;
        bin_size = dist_store->box_width/dist_store->length;
        if (adt < 0 || adt > dist_store->nframes) {
            _average_field_dist(dist_store);
            _calculate_thickness_dist(dist_store);
        }
        for (i=0; i < dist_store->length; ++i) {
            dist_store->height[2][i] /= dist_store->sampling[2][i];
        }
        /* Write the output */
        for (i=0; i<dist_store->length; ++i) {
            if (dist_store->sampling[2][i] > 0) {
                fprintf(dist_store->out_dist, "%7.3f %7.3f\n",
                        i*bin_size, dist_store->height[2][i]);
                fprintf(dist_store->out_sampling, "%7.3f %7d\n",
                        i*bin_size, dist_store->sampling[2][i]);
            }
        }
    }
}
