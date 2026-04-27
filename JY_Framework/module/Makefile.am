AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(OFFLINE_MAIN)/include/phool \
  -I$(OFFLINE_MAIN)/include/fun4all \
  -I$(OFFLINE_MAIN)/include/trackbase \
  -I$(OFFLINE_MAIN)/include/trackbase_historic \
  -I$(OFFLINE_MAIN)/include/calobase \
  -I$(OFFLINE_MAIN)/include/globalvertex \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  SiCaloMatcher.h \
  SiCaloMatchedTrack.h \
  SiCaloMatchedTrackContainer.h

lib_LTLIBRARIES = \
  libtrackcalomatch.la

libtrackcalomatch_la_SOURCES = \
  SiCaloMatcher.cc \
  SiCaloMatchedTrack.cc \
  SiCaloMatchedTrackContainer.cc

libtrackcalomatch_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -ltrack_io \
  -ltrack \
  -lcalo_io \
  -lg4detectors \
  -lg4detectors_io \
  -lcdbobjects \
  -lphg4hit \
  -lg4eval \
  -lCLHEP

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libtrackcalomatch.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)