AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(top_srcdir) \
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
  $(top_srcdir)/SiCaloMatcher.h \
  $(top_srcdir)/SiCaloMatchedTrack.h \
  $(top_srcdir)/SiCaloMatchedTrackContainer.h

lib_LTLIBRARIES = \
  libtrackcalomatch.la

libtrackcalomatch_la_SOURCES = \
  $(top_srcdir)/SiCaloMatcher.cc \
  $(top_srcdir)/SiCaloMatchedTrack.cc \
  $(top_srcdir)/SiCaloMatchedTrackContainer.cc \
  trackcalomatch_Dict.cc

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

trackcalomatch_Dict.cc: $(top_srcdir)/SiCaloMatchedTrack.h $(top_srcdir)/SiCaloMatchedTrackContainer.h $(top_srcdir)/LinkDef.h
	rootcling -f $@ $(CINTDEFS) $(AM_CPPFLAGS) \
	  $(top_srcdir)/SiCaloMatchedTrack.h \
	  $(top_srcdir)/SiCaloMatchedTrackContainer.h \
	  $(top_srcdir)/LinkDef.h

BUILT_SOURCES = \
  testexternals.cc \
  trackcalomatch_Dict.cc

CLEANFILES = \
  trackcalomatch_Dict.cc \
  trackcalomatch_Dict.h \
  trackcalomatch_Dict_rdict.pcm \
  libtrackcalomatch_rdict.pcm

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
	rm -f $(BUILT_SOURCES) $(CLEANFILES)