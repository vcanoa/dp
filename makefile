all:	libDPLocal.so pairs_FG pairs_BG

clean:
	rm libDPLocal.so pairs_FG

pairs_FG: pairs_FG.cc
	`root-config --cxx` `root-config --cflags --libs` -L. -lDPLocal $^ -o pairs_FG

pairs_BG: pairs_BG.cc
	`root-config --cxx` `root-config --cflags --libs` -L. -lDPLocal $^ -o pairs_BG

libDPLocal.so: PhotonEvent.o EMCC.o CNTE.o CNTDE.o dpUtil.o dpReco.o
	 `root-config --cxx` -shared -WL,-soname,libDPLocal.so `root-config --cflags --libs` $^ -o libDPLocal.so 
	rm $^

PhotonEvent.o: PhotonEvent.cxx
	`root-config --cxx` `root-config --cflags --libs` -c $^

EMCC.o:EMCC.cxx
	`root-config --cxx` `root-config --cflags --libs` -c $^
CNTE.o:CNTE.cxx 	
	`root-config --cxx` `root-config --cflags --libs` -c $^
CNTDE.o:CNTDE.cxx
	`root-config --cxx` `root-config --cflags --libs` -c $^
dpUtil.o:dpUtil.cxx
	`root-config --cxx` `root-config --cflags --libs` -c $^
dpReco.o:dpReco.cxx
	`root-config --cxx` `root-config --cflags --libs` -c $^
