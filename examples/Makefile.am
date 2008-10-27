INCLUDES= -I. -I.. 

AM_LDFLAGS = ../lib/libmind.la $(LIBS)

EXTRA_DIST = run_examples.sh 

bin_PROGRAMS =  root_2_dst fit_tracks simple_smear

root_2_dst_SOURCES = root_2_dst.cpp 

fit_tracks_SOURCES = fit_tracks.cpp

simple_smear_SOURCES = simple_smear.cpp

examples: $(bin_PROGRAMS)

dst: root_2_dst

fit: fit_tracks

smear: simple_smear

CLEANFILES = *~* 




