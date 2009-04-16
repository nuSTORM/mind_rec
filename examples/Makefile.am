INCLUDES= -I. -I.. 

AM_LDFLAGS = ../lib/libmind.la $(LIBS)

EXTRA_DIST = run_examples.sh 

bin_PROGRAMS =  fit_tracks

fit_tracks_SOURCES = fit_tracks.cpp

examples: $(bin_PROGRAMS)

fit: fit_tracks

CLEANFILES = *~* 




