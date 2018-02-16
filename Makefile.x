#
# Makefile for XVMC
#
CC=g++

# debug and development options
#  CFLAGS=-g -Wall -I. -c -DALL_INLINE  -DUSE_SOBOL
#CFLAGS=-g -Wall -I. -c -DDEBUG -DDOUBLE_REAL -DUSE_SOBOL
#CFLAGS=-g -Wall -I. -c -DDEBUG
#CFLAGS=-g -Wall -I. -c -DDEBUG -DUSE_SOBOL
#CFLAGS=-g -Wall -I. -c
#LFLAGS=-g -Wall -lpthread -lmpatrol -lbfd -liberty
#  LFLAGS=-g -Wall -lpthread

# optimization options
#CFLAGS=-O2 -I. -c -DALL_INLINE -DUSE_SOBOL -DONE_PROCESS -DICCR_TEST
#CFLAGS=-O2 -I. -c -DALL_INLINE -DDOUBLE_REAL -DUSE_SOBOL
#CFLAGS=-O2 -I. -c -DALL_INLINE -DCHECK_ENERGY
#CFLAGS=-O2 -I. -c -DALL_INLINE -DCHECK_ENERGY -DDOUBLE_REAL -DUSE_SOBOL
#CFLAGS=-O2 -I. -c -DALL_INLINE -DCHECK_ENERGY -DUSE_SOBOL -DONE_PROCESS
#CFLAGS=-O2 -ffast-math -I. -c -DALL_INLINE -DCHECK_ENERGY -DUSE_SOBOL -DICCR_TEST
CFLAGS=-O2 -ffast-math -I. -c -DALL_INLINE -DUSE_SOBOL # -DCMS_MOD
LFLAGS=-g -O2 -ffast-math -lpthread

CFLAGS= -ffast-math -I. -c -DALL_INLINE -DREF_POINT # -DUSE_SOBOL -DCMS_MOD
LFLAGS= -ffast-math -pthread 

OBJECTS= xvmc.o xvmc_util.o xvmc_aux.o xvmc_messages.o ranmar.o sobseq.o \
         contour.o multi_leaf.o init_beam.o io_routines.o \
         array_2d.o array_3d.o init_xvmc.o InitCrossSectionTables.o \
         beam_core.o treatment_plan.o \
         MC_plane.o MC_halfcyl.o MC_material.o \
         MC_region.o MC_slab.o MC_volume_6p.o \
         MC_object.o MC_doubleslab.o MC_jaws_focus.o MC_mlc.o \
         MC_mlc_2focus.o MC_mlc_rfocus.o MC_mlc_elekta.o MC_mlc_varian.o \
         beam_modifier.o irregular_modifier.o simple_mlc.o real_mlc.o \
         beam_model.o p_beam_model.o \
         e_beam_1point_mono.o e_beam_1point_poly.o \
         e_beam_triple_mono.o e_beam_triple_poly.o \
         e_beam_bumpy_mono.o  e_beam_bumpy_poly.o \
         p_beam_1point_mono.o p_beam_1point_spec.o \
         p_beam_2gauss_mono.o p_beam_2gauss_poly.o \
         electron_data_inp.o electron_data.o electron_transport_data.o \
         photon_data_inp.o photon_data.o \
         moller_XS.o bhabha_XS.o brems_XS.o mscat.o \
         compton_XS_diff.o compton_XS_total.o \
         pair_XS_diff.o pair_XS_total.o photo_XS_total.o \
         calc_dose.o portal_dose.o evaluate.o \
         multi_electron.o one_electron.o multi_photon.o kerma_photon.o \
         write_bmp.o

HEADERS= definitions.h global.h xvmc.h xvmc_util.h ranmar.h sobseq.h \
         contour.h multi_leaf.h array_2d.h array_3d.h \
         beam_core.h treatment_plan.h MCBitSet.h \
         MC_plane.h MC_halfcyl.h MC_material.h \
         MC_region.h MC_slab.h MC_volume_6p.h \
         MC_object.h MC_doubleslab.h MC_jaws_focus.h MC_mlc.h \
         MC_mlc_2focus.h MC_mlc_rfocus.h MC_mlc_elekta.h MC_mlc_varian.h \
         beam_modifier.h irregular_modifier.h simple_mlc.h real_mlc.h \
         beam_model.h e_beam_model.h p_beam_model.h \
         e_beam_1point_mono.h e_beam_1point_poly.h \
         e_beam_triple_mono.h e_beam_triple_poly.h \
         e_beam_bumpy_mono.h  e_beam_bumpy_poly.h \
         p_beam_1point_mono.h p_beam_1point_spec.h \
         p_beam_2gauss_mono.h p_beam_2gauss_poly.h \
         electron_data_inp.h electron_data.h electron_transport_data.h \
         photon_data_inp.h photon_data.h \
         moller_XS.h bhabha_XS.h brems_XS.h \
         compton_XS_diff.h compton_XS_total.h \
         pair_XS_diff.h pair_XS_total.h photo_XS_total.h \
         multi_electron.h portal_dose.h

all: xvmc

xvmc: $(OBJECTS) 
	$(CC) $(LFLAGS) -o xvmc $(OBJECTS)

xvmc_orig: $(OBJECTS) p_beam_2gauss_mono.o calc_dose.o init_xvmc.o io_routines.o
	$(CC) $(LFLAGS) -o xvmc $(OBJECTS)

xvmc_moving_jaw: $(OBJECTS) p_beam_2gauss_mono_moving_jaw.o calc_dose.o init_xvmc.o io_routines.o
	$(CC) $(LFLAGS) -o xvmc_moving_jaw $(OBJECTS)
	
xvmc_nu_value: $(OBJECTS) p_beam_2gauss_mono_nu_value.o io_routines.o
	$(CC) $(LFLAGS) -o xvmc_nu $(OBJECTS)

xvmc_phps_read: $(OBJECTS) calc_dose_phsp_read.o init_xvmc_phsp_read.o io_routines.o
	$(CC) $(LFLAGS) -o xvmc_phps_read $(OBJECTS)

xvmc.o: xvmc.cpp $(HEADERS)
	$(CC) $(CFLAGS) xvmc.cpp

xvmc_util.o: xvmc_util.cpp $(HEADERS)
	$(CC) $(CFLAGS) xvmc_util.cpp

xvmc_aux.o: xvmc_aux.cpp $(HEADERS)
	$(CC) $(CFLAGS) xvmc_aux.cpp

xvmc_messages.o: xvmc_messages.cpp $(HEADERS)
	$(CC) $(CFLAGS) xvmc_messages.cpp

ranmar.o: ranmar.cpp $(HEADERS)
	$(CC) $(CFLAGS) ranmar.cpp

sobseq.o: sobseq.cpp $(HEADERS)
	$(CC) $(CFLAGS) sobseq.cpp

contour.o: contour.cpp $(HEADERS)
	$(CC) $(CFLAGS) contour.cpp

multi_leaf.o: multi_leaf.cpp $(HEADERS)
	$(CC) $(CFLAGS) multi_leaf.cpp

io_routines.o: io_routines.cpp $(HEADERS)
	$(CC) $(CFLAGS) io_routines.cpp

array_2d.o: array_2d.cpp $(HEADERS)
	$(CC) $(CFLAGS) array_2d.cpp

array_3d.o: array_3d.cpp $(HEADERS)
	$(CC) $(CFLAGS) array_3d.cpp

init_xvmc.o: init_xvmc.cpp $(HEADERS)
	$(CC) $(CFLAGS) init_xvmc.cpp

InitCrossSectionTables.o: InitCrossSectionTables.cpp $(HEADERS)
	$(CC) $(CFLAGS) InitCrossSectionTables.cpp

beam_core.o: beam_core.cpp $(HEADERS)
	$(CC) $(CFLAGS) beam_core.cpp

treatment_plan.o: treatment_plan.cpp $(HEADERS)
	$(CC) $(CFLAGS) treatment_plan.cpp

init_beam.o: init_beam.cpp $(HEADERS)
	$(CC) $(CFLAGS) init_beam.cpp

MC_plane.o: MC_plane.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_plane.cpp

MC_halfcyl.o: MC_halfcyl.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_halfcyl.cpp

MC_material.o: MC_material.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_material.cpp

MC_region.o: MC_region.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_region.cpp

MC_slab.o: MC_slab.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_slab.cpp

MC_volume_6p.o: MC_volume_6p.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_volume_6p.cpp

MC_object.o: MC_object.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_object.cpp

MC_doubleslab.o: MC_doubleslab.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_doubleslab.cpp

MC_jaws_focus.o: MC_jaws_focus.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_jaws_focus.cpp

MC_mlc.o: MC_mlc.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_mlc.cpp

MC_mlc_2focus.o: MC_mlc_2focus.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_mlc_2focus.cpp

MC_mlc_rfocus.o: MC_mlc_rfocus.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_mlc_rfocus.cpp

MC_mlc_elekta.o: MC_mlc_elekta.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_mlc_elekta.cpp

MC_mlc_varian.o: MC_mlc_varian.cpp $(HEADERS)
	$(CC) $(CFLAGS) MC_mlc_varian.cpp

beam_modifier.o: beam_modifier.cpp $(HEADERS)
	$(CC) $(CFLAGS) beam_modifier.cpp

irregular_modifier.o: irregular_modifier.cpp $(HEADERS)
	$(CC) $(CFLAGS) irregular_modifier.cpp

simple_mlc.o: simple_mlc.cpp $(HEADERS)
	$(CC) $(CFLAGS) simple_mlc.cpp

real_mlc.o: real_mlc.cpp $(HEADERS)
	$(CC) $(CFLAGS) real_mlc.cpp

beam_model.o: beam_model.cpp $(HEADERS)
	$(CC) $(CFLAGS) beam_model.cpp

e_beam_1point_mono.o: e_beam_1point_mono.cpp $(HEADERS)
	$(CC) $(CFLAGS) e_beam_1point_mono.cpp

e_beam_1point_poly.o: e_beam_1point_poly.cpp $(HEADERS)
	$(CC) $(CFLAGS) e_beam_1point_poly.cpp

e_beam_triple_mono.o: e_beam_triple_mono.cpp $(HEADERS)
	$(CC) $(CFLAGS) e_beam_triple_mono.cpp

e_beam_triple_poly.o: e_beam_triple_poly.cpp $(HEADERS)
	$(CC) $(CFLAGS) e_beam_triple_poly.cpp

e_beam_bumpy_mono.o: e_beam_bumpy_mono.cpp $(HEADERS)
	$(CC) $(CFLAGS) e_beam_bumpy_mono.cpp

e_beam_bumpy_poly.o: e_beam_bumpy_poly.cpp $(HEADERS)
	$(CC) $(CFLAGS) e_beam_bumpy_poly.cpp

p_beam_model.o: p_beam_model.cpp $(HEADERS)
	$(CC) $(CFLAGS) p_beam_model.cpp

p_beam_1point_mono.o: p_beam_1point_mono.cpp $(HEADERS)
	$(CC) $(CFLAGS) p_beam_1point_mono.cpp

p_beam_1point_spec.o: p_beam_1point_spec.cpp $(HEADERS)
	$(CC) $(CFLAGS) p_beam_1point_spec.cpp

p_beam_2gauss_mono.o: p_beam_2gauss_mono.cpp $(HEADERS)
	$(CC) $(CFLAGS) p_beam_2gauss_mono.cpp

p_beam_2gauss_poly.o: p_beam_2gauss_poly.cpp $(HEADERS)
	$(CC) $(CFLAGS) p_beam_2gauss_poly.cpp
	
electron_data_inp.o: electron_data_inp.cpp $(HEADERS)
	$(CC) $(CFLAGS) electron_data_inp.cpp

electron_data.o: electron_data.cpp $(HEADERS)
	$(CC) $(CFLAGS) electron_data.cpp

photon_data_inp.o: photon_data_inp.cpp $(HEADERS)
	$(CC) $(CFLAGS) photon_data_inp.cpp

photon_data.o: photon_data.cpp $(HEADERS)
	$(CC) $(CFLAGS) photon_data.cpp

electron_transport_data.o: electron_transport_data.cpp $(HEADERS)
	$(CC) $(CFLAGS) electron_transport_data.cpp

moller_XS.o: moller_XS.cpp $(HEADERS)
	$(CC) $(CFLAGS) moller_XS.cpp

bhabha_XS.o: bhabha_XS.cpp $(HEADERS)
	$(CC) $(CFLAGS) bhabha_XS.cpp

brems_XS.o: brems_XS.cpp $(HEADERS)
	$(CC) $(CFLAGS) brems_XS.cpp

mscat.o: mscat.cpp $(HEADERS)
	$(CC) $(CFLAGS) mscat.cpp

compton_XS_diff.o: compton_XS_diff.cpp $(HEADERS)
	$(CC) $(CFLAGS) compton_XS_diff.cpp

compton_XS_total.o: compton_XS_total.cpp $(HEADERS)
	$(CC) $(CFLAGS) compton_XS_total.cpp

pair_XS_diff.o: pair_XS_diff.cpp $(HEADERS)
	$(CC) $(CFLAGS) pair_XS_diff.cpp

pair_XS_total.o: pair_XS_total.cpp $(HEADERS)
	$(CC) $(CFLAGS) pair_XS_total.cpp

photo_XS_total.o: photo_XS_total.cpp $(HEADERS)
	$(CC) $(CFLAGS) photo_XS_total.cpp

calc_dose.o: calc_dose.cpp $(HEADERS)
	$(CC) $(CFLAGS) -D_REENTRANT calc_dose.cpp

portal_dose.o: portal_dose.cpp $(HEADERS)
	$(CC) $(CFLAGS) portal_dose.cpp

multi_electron.o: multi_electron.cpp $(HEADERS)
	$(CC) $(CFLAGS) multi_electron.cpp

one_electron.o: one_electron.cpp $(HEADERS)
	$(CC) $(CFLAGS) one_electron.cpp

multi_photon.o: multi_photon.cpp $(HEADERS)
	$(CC) $(CFLAGS) multi_photon.cpp

kerma_photon.o: kerma_photon.cpp $(HEADERS)
	$(CC) $(CFLAGS) kerma_photon.cpp

evaluate.o: evaluate.cpp $(HEADERS)
	$(CC) $(CFLAGS) evaluate.cpp

write_bmp.o: write_bmp.c
	gcc -c -Wall write_bmp.c

p_beam_2gauss_mono_orig.o: p_beam_2gauss_mono.cpp $(HEADERS)
	$(CC) $(CFLAGS) p_beam_2gauss_mono.cpp

p_beam_2gauss_mono_nu_value.o: p_beam_2gauss_mono.cpp $(HEADERS)
	$(CC) $(CFLAGS) p_beam_2gauss_mono.cpp -DUSE_NU

p_beam_2gauss_mono_moving_jaw.o: p_beam_2gauss_mono.cpp $(HEADERS)
	$(CC) $(CFLAGS) p_beam_2gauss_mono.cpp -DMOVING_JAW

init_xvmc_phsp_read.o: init_xvmc.cpp $(HEADERS)
	$(CC) $(CFLAGS) init_xvmc.cpp -DPHSP_READ

calc_dose_phsp_read.o: calc_dose.cpp $(HEADERS)
	$(CC) $(CFLAGS) -D_REENTRANT calc_dose.cpp -DPHSP_READ

clean:
	rm -f *.o
	