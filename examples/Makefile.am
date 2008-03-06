INCLUDES= -I. -I.. 

AM_LDFLAGS = ../lib/libmind.la $(LIBS)

EXTRA_DIST = run_examples.sh 

bin_PROGRAMS =  root_2_dst fit_tracks

root_2_dst_SOURCES = root_2_dst.cpp 

fit_tracks_SOURCES = fit_tracks.cpp

examples: $(bin_PROGRAMS)

dst: root_2_dst

fit: fit_tracks

CLEANFILES = *~* 




