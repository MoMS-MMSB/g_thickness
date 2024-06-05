#include "grid_mode.h"

void _average_field(GridHeight *grid_store) {
    int leaflet, i, j;
    for (leaflet = 0; leaflet < 2; ++leaflet) {
        for (i=0; i < grid_store->shape[0]; ++i) {
            for (j=0; j < grid_store->shape[1]; ++j) {
                grid_store->grids[leaflet][i][j] /=
                    grid_store->sampling[leaflet][i][j];
            }
        }
    }
}

void _calculate_thickness(GridHeight *grid_store) {
    int i, j;
    int minsamp = 0;
    for (i=0; i < grid_store->shape[0]; ++i) {
        for (j=0; j < grid_store->shape[1]; ++j) {
            minsamp = min(grid_store->sampling[0][i][j],
                    grid_store->sampling[1][i][j]);
            if (minsamp > 0) {
                grid_store->grids[2][i][j] +=
                    (real)fabs(grid_store->grids[0][i][j] - 
                            grid_store->grids[1][i][j]) * minsamp;
                grid_store->sampling[2][i][j] += minsamp;
            }
        }
    }
}

/** Contruct an instance of GridHeight
 *
 * All the dimensions described in the "shape" array have to be greater than 0.
 */
GridHeight *build_grids(int shape[2], int normal_axis,
        const char *grid_fn, const char *sampling_fn) {
    GridHeight *grid_store;
    int grid, i;

    /* Check dimensions */
    if (shape[0] <= 0 || shape[1] <= 0) {
        fprintf(stderr,
                "I can not build a grid with this dimensions: (%d, %d)\n",
                shape[0], shape[1]);
        exit(1);
    }

    /* Allocate empty structure */
    snew(grid_store, 1);
    
    grid_store->nframes = 0;
    for (i=0; i<2; ++i) {
        /* Store the shape */
        grid_store->shape[i] = shape[i];
        /* Set the slice width to 0, just in case */
        grid_store->width[i] = 0;
        /* Set box_width sommation to 0 */
        grid_store->box_width[i] = 0.0;
    }
    /* Define the axis */
    grid_store->axis[0] = normal_axis;
    switch (normal_axis) {
        case 0:
            grid_store->axis[1] = 1; grid_store->axis[2] = 2;
            break;
        case 1:
            grid_store->axis[1] = 0; grid_store->axis[2] = 2;
            break;
        case 2:
            grid_store->axis[1] = 0; grid_store->axis[2] = 1;
            break;
        default:
            gmx_fatal(FARGS,"Invalid axes. Terminating. \n");
    }

    /* Allocate the grids */
    for (grid = 0; grid < 3; ++grid) {
        (grid_store->grids)[grid] = realMatrix(shape[0], shape[1], 0.0);
        (grid_store->sampling)[grid] = intMatrix(shape[0], shape[1], 0);
    }

    /* Open the files */
    grid_store->out_grid = ffopen(grid_fn, "w");
    if (grid_store->out_grid == NULL) {
        fprintf(stderr, "Error oppenning %s for grid mode\n", grid_fn);
        exit(1);
    }
    grid_store->out_sampling = ffopen(sampling_fn, "w");
    if (grid_store->out_sampling == NULL && sampling_fn != NULL) {
        fprintf(stderr, "Error oppenning %s for grid mode\n", sampling_fn);
        exit(1);
    }

    return grid_store;
}

/** Clean an instance of GridHeight
 */
void clean_grids(GridHeight *grid_store) {
    int grid = 0;
    if (grid_store) {
        for (grid = 0; grid < 3; ++grid) {
            deleteRealMat(grid_store->grids[grid], grid_store->shape[0]);
            deleteIntMat(grid_store->sampling[grid], grid_store->shape[0]);
        }
        ffclose(grid_store->out_grid);
        ffclose(grid_store->out_sampling);
        sfree(grid_store);
    }
}

void grid_start_frame(GridHeight *grid_store, matrix box) {
    int i = 0;
    int axis = 0;
    if (grid_store) {
        grid_store->nframes += 1;
        for (i=0; i<2; ++i) {
            axis = grid_store->axis[i+1];
            grid_store->width[i] = box[axis][axis]/grid_store->shape[i];
            grid_store->box_width[i] += box[axis][axis];
        }
    }
}

void grid_end_frame(GridHeight *grid_store, int adt) {
    int i, j, leaflet;
    int minsamp=0;
    if (grid_store && adt > 0 && grid_store->nframes % adt == 0) {
        /* Average the fields */
        _average_field(grid_store);
        /* Calculate the thickness */
        _calculate_thickness(grid_store);
        /* Empty the fields */
        for (leaflet = 0; leaflet < 2; ++leaflet) {
            for (i=0; i < grid_store->shape[0]; ++i) {
                for (j=0; j < grid_store->shape[1]; ++j) {
                    grid_store->grids[leaflet][i][j] = 0.0;
                    grid_store->sampling[leaflet][i][j] = 0;
                }
            }
        }
    }
}

void grid_store(GridHeight *grid, int leaflet, rvec atom, t_pbc *pbc) {
    int axis = 0;
    int slice[2] = {0, 0};
    int i = 0;
    if (grid) {
        put_atom_in_box((real (*)[3])pbc->box,atom);
        axis = grid->axis[0];
        for (i =0; i<2; ++i) {
            slice[i] = atom[grid->axis[i+1]]/grid->width[i];
        }
        grid->grids[leaflet][slice[0]][slice[1]] += atom[axis];
        grid->sampling[leaflet][slice[0]][slice[1]] += 1;
    }
}

void grid_end(GridHeight *grid_store, int adt) {
    char labels[] = "XYZ";
    if (grid_store) {
        int i, j, leaflet;
        if (adt < 0 || adt > grid_store->nframes) {
            _average_field(grid_store);
            _calculate_thickness(grid_store);
        }
        for (i=0; i < grid_store->shape[0]; ++i) {
            for (j=0; j < grid_store->shape[1]; ++j) {
                grid_store->grids[2][i][j] /= grid_store->sampling[2][i][j];
            }
        }
        /* Write the output */
        fprintf(grid_store->out_grid, "@xwidth %7.3f\n",
                grid_store->box_width[0]/grid_store->nframes);
        fprintf(grid_store->out_grid, "@ywidth %7.3f\n",
                grid_store->box_width[1]/grid_store->nframes);
        fprintf(grid_store->out_sampling, "@xwidth %7.3f\n",
                grid_store->box_width[0]/grid_store->nframes);
        fprintf(grid_store->out_sampling, "@ywidth %7.3f\n",
                grid_store->box_width[1]/grid_store->nframes);
        fprintf(grid_store->out_grid, "@xlabel %c (nm)\n",
                labels[grid_store->axis[1]]);
        fprintf(grid_store->out_grid, "@ylabel %c (nm)\n",
                labels[grid_store->axis[2]]);
        fprintf(grid_store->out_sampling, "@xlabel %c (nm)\n",
                labels[grid_store->axis[1]]);
        fprintf(grid_store->out_sampling, "@ylabel %c (nm)\n",
                labels[grid_store->axis[2]]);
        fprintf(grid_store->out_grid, "@legend Thickness (nm)\n");
        fprintf(grid_store->out_sampling, "@legend Thickness (nm)\n");
        for (leaflet = 2; leaflet < 3; ++leaflet) {
            for (i=0; i < grid_store->shape[0]; ++i) {
                for (j=0; j < grid_store->shape[1]; ++j) {
                    if (j > 0) {
                        fprintf(grid_store->out_grid, "\t");
                        fprintf(grid_store->out_sampling, "\t");
                    }
                    fprintf(grid_store->out_grid, "%7.3f",
                            grid_store->grids[leaflet][i][j]);
                    fprintf(grid_store->out_sampling, "%d",
                            grid_store->sampling[leaflet][i][j]);
                }
                fprintf(grid_store->out_grid, "\n");
                fprintf(grid_store->out_sampling, "\n");
            }
        }
    }
}


