INCLUDES= -I. -I.. 

AM_LDFLAGS = ../lib/libmind.la $(LIBS)

EXTRA_DIST = run_examples.sh 

bin_PROGRAMS =  root_2_dst

root_2_dst_SOURCES = root_2_dst.cpp 

examples: $(bin_PROGRAMS)

dst: root_2_dst

CLEANFILES = *~* 




