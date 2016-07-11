AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include

lib_LTLIBRARIES = \
libDPRun16dAu.la

noinst_HEADERS = \
DPRun16dAuLinkDef.h

pkginclude_HEADERS = \
PhotonEvent.h \
EMCC.h \
CNTE.h \
CNTDE.h \
Reconstruction.h \
RDanalyzer.h

libDPRun16dAu_la_SOURCES = \
PhotonEvent.cxx \
EMCC.cxx \
CNTE.cxx \
CNTDE.cxx \
Reconstruction.cpp \
DPRun16dAu_Dict.C 

libDPRun16dAu_la_LIBADD = \
-L$(libdir) \
-L$(OFFLINE_MAIN)/lib \
-lSubsysReco \
-luti \
-lHardpDST \
-lCNT \
-lcgl \
-lpad \
-lemc \
-llpc \
-lTHmul \
-lboost_thread \
-lboost_date_time \
-lTOAD \
-lFROG \
-lfun4allfuncs

# Rule for generating CINT dictionaries from class headers.
DPRun16dAu_Dict.C: \
  PhotonEvent.h \
  EMCC.h \
  CNTE.h \
  CNTDE.h \
  Reconstruction.h \
  DPRun16dAuLinkDef.h 
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
