# g_thickness: Calculate bilayer local thickness

![maintenance-status](https://img.shields.io/badge/maintenance-deprecated-red.svg)

## Installation
``g_thickness`` depends on the GROMACS software package, which needs to be
installed.  Only versions 4.5.x have been tested, but ``g_thickness`` might be
compatible with other versions of GROMACS.

To install ``g_thickness``, GROMACS needs to be loaded. You can load
it using:

    source /path_to_gromacs/bin/GMXRC

Go into the source directory of the program, then run ``make``. The
``g_thickness`` executable should be created.  Make sure that
this executable is in the research path of your shell.

## Usage
Here we assume that ``g_thickness`` is in the research path of your shell. To
get some help just run ``g_thickness -h``. All available options will be
listed.

A classical use would be:

    g_thickness -f traj.xtc -s topol.tpr -n index.ndx -og thickness_grid.dat

GROMACS needs to be loaded for ``g_thickness`` to work.

### Required arguments
Three arguments are required for any use of the program. They are:

* ``-f``: the path to the trajectory to read;
* ``-s``: the path to the topology (tpr file);
* ``-n``: the path to an index file that describe the group of atoms you are
  interested in;

The ``-d`` option defines the reference axis (i.e. the axis normal to the
membrane). The axis is set at Z by default.

### Output control
The following arguments control the output. You can get either one or both of
the possible outputs but you need to select at least one of them.

* ``-og``: produce the thickness landscape; a file path can be given as
  argument. The sampling for each cell of the grid is written in the file
  specified with the ``-osg`` option. The produced files can be converted into
  a picture. See the `Generate picture from landscapes`_ section to know more
  about that.
* ``-od``: produce the thickness profile as a function of distance to a group.
  The sampling is written in the file given with the ``-osd`` option. Distance
  is calculated as a function of the center of mass of a reference group; the
  minimum distance can be used instead using the ``-nocom`` option.

### Sampling control
There is two ways to adjust the sampling. The ``-sl`` option corresponds to the
number of bins in the profile or to the number of cell per side in the
landscape. The ``-sl2`` option adjusts the number of cell per side on the
second dimension of a profile; if it is not used then the value of ``-sl`` is
used. See the following table to know which axis are controlled by ``-sl`` and
``-sl2``.
```
====== ======= ========
``-d`` ``-sl`` ``-sl2``
====== ======= ========
x      y       z
y      x       z
z      x       y
====== ======= ========
```
The thickness averaging in a given cell or bin is done as follow: at each frame
the atoms of interest in one leaflet in the considered volume are stored for
this volume; the same is done for the other leaflet;  every N frames, the
center of mass of the stored atoms is calculated for each leaflet, then the
distance between the two center of mass along the membrane normal axis is
calculated and stored; at the end of the calculation, the distances are
averaged. The number of frames N can be chosen using the ``-adt`` option. The
bigger this value will be the bigger will be the sampling in each cell or bin,
but the less accurate will be the result. If the value of ``-adt`` is negative
or bigger than the number of frame in the trajectory, then distance between
leaflets is calculated only once at the end.

### Generate pictures from landscapes
The landscape output is a text file describing the order parameter values on a
grid. The file format is not XPM like most grid outputs produced by GROMACS
tools so the ``xpm2ps`` utility can not be used to produce usable pictures. The
``dispgrid`` python script aims to exploit the data and to produce pictures
from them.

You need python 2.7, with the numpy and matplotlib modules to run dispgrid.

Basic usage of dispgrid is:

    dispgrid input.dat output.png

See the help available by typing ``dispgrid -h`` for more features.
