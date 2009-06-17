INCLUDES= -I. -I.. 

AM_LDFLAGS = ../lib/libmind.la $(LIBS)

EXTRA_DIST = run_examples.sh 

bin_PROGRAMS =  fit_tracks rec_had_show

fit_tracks_SOURCES = fit_tracks.cpp

rec_had_show_SOURCES =  rec_had_show.cpp

examples: $(bin_PROGRAMS)

fit: fit_tracks

hadrec: rec_had_show

CLEANFILES = *~* 




