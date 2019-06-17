art_make_exec( NAME testEarthModel
               SOURCE testEarthModel.cc
               LIBRARIES OscLibFunc
	     )
art_make_exec( NAME testOsc
               SOURCE testOsc.cc
               LIBRARIES OscLibFunc
	     )
art_make_exec( NAME testOscRndm
               SOURCE testOscRndm.cc
               LIBRARIES OscLibFunc
	     )

install_headers()
install_source()
