# This is a Gromacs 4.5 template makefile for your own utility programs using pkg-config.
#
# Copy this file to whatever directory you are using for your own software
#
# Usage:
# $ source /path/to/GMXRC
# $ make -f Makefile.pkg
#
#change the name of the program here
NAME=g_thickness

#add extra c file to compile here
EXTRA_SRC=matrix.c distances.c dist_mode.c grid_mode.c

###############################################################3
#below only boring default stuff
#only change it if you know what you are doing ;-)

#what should be done by default
all: $(NAME)

#if GMXLDLIB is defined we add it to PKG_CONFIG_PATH
ifeq "$(origin GMXLDLIB)" "undefined"
  $(error "GMXLDLIB not found, please source GMXRC")
else
  export PKG_CONFIG_PATH:=${PKG_CONFIG_PATH}:${GMXLDLIB}/pkgconfig
endif

#get CPPFLAGS and LDFLAGS from pkg-config
CPPFLAGS=`pkg-config --cflags libgmx`
LDFLAGS=`pkg-config --libs libgmx`

#generate a list of object (.o) files
OBJS=$(patsubst %.c,%.o,$(NAME).c $(EXTRA_SRC))

#main program depend on all objects, rest is done by implicit rules
$(NAME): $(OBJS)

%.o: %.c
	cc  `pkg-config --cflags libgmx`  -c -o $@ $<

g_thickness: distances.o dist_mode.o grid_mode.o matrix.o g_thickness.o
	cc $^ -o $@ `pkg-config --libs libgmx` -lmd -lgmx -lm -ldl -lpthread


#clean up rule
clean:
	rm -f $(NAME) $(OBJS)

#all, clean are phony rules, e.g. they are always run
.PHONY: all clean
