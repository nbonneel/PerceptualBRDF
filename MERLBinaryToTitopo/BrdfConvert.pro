TEMPLATE = app

message("QMake rules: Application")

contains( CONFIG, "debug" ) {
	contains( CONFIG, "release" ) {
		error("CONFIG cannot contains both release and debug flags.") 
	}
}

#CONFIG *= openmp
LIBS *= -lgomp

INCLUDEPATH += . .. ../..

DEFINES *= PAS_NICO

macx{
    INCLUDEPATH += /opt/local/include ../local/include
}

QMAKE_CXXFLAGS *= -frounding-math -Wall -fpermissive -fPIE

# Modif JPF - SSE fp instructions on osX - Needed for exceptions
macx{
    QMAKE_CXXFLAGS += -mfpmath=sse
	QMAKE_CXXFLAGS -= -fopenmp
}

debug {
        DEFINES *= DEBUG
		OBJECTS_DIR = .obj.debug
		DESTDIR = ../bin.debug
# commented out because it kills performance drastically
#		QMAKE_CXXFLAGS += -fsanitize=address
#		LIBS *= -fsanitize=address
		LIBS *= -L../lib.debug 
		PRE_TARGETDEPS *= ../lib.debug/libBrdf.a 
}
release {
		OBJECTS_DIR = .obj.release
		DESTDIR = ../bin.release
		LIBS *= -L../lib.release 
		DEFINES *= NDEBUG
		PRE_TARGETDEPS *= ../lib.release/libBrdf.a 
		QMAKE_CXXFLAGS *= -g
}

macx{
	LIBS -= -fopenmp
}

openmp {
	QMAKE_CXXFLAGS *= -fopenmp
	LIBS *= -lgomp
}



SOURCES = IsotropicMERLBRDF.cpp  main.cpp  BrdfIO.cpp
HEADERS = IsotropicMERLBRDF.h MERLBRDFRead.h BrdfIO.h BRDF.h argstream.h Spectrum.h


