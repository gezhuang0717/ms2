TARGET = ../lib/liboptics.so
OBJS = encoptics.o loadopticsprm.o

include ../common.mk

loadopticsprm.o: optics.fh
encoptics.o: optics.fh ../traj/trajectory.fh
