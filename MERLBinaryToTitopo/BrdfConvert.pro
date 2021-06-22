TEMPLATE = lib
CONFIG *= staticlib 
#CONFIG *= openmp

QMAKE_CXXFLAGS *= -fPIC -g -frounding-math -Wall

# Modif JPF - SSE fp instructions on osX - Needed for exceptions
macx{
    QMAKE_CXXFLAGS += -mfpmath=sse
}

contains( CONFIG, "debug" ) {
	contains( CONFIG, "release" ) {
		error("CONFIG cannot contains both release and debug flags. Please use 'qmake CONFIG=debug' or 'qmake CONFIG=release'") 
	}
}

DEFINES *= PAS_NICO
INCLUDEPATH += . .. ../..

# Modif JPF - added include directories for osX with macports and local ylm 
macx{
    INCLUDEPATH += /opt/local/include ../local/include
}

debug {
		DEFINES *= DEBUG
		OBJECTS_DIR = .obj.debug
		DESTDIR = ../lib.debug
}
release {
		OBJECTS_DIR = .obj.release
		DESTDIR = ../lib.release
		DEFINES *= NDEBUG
}


SOURCES = IsotropicMERLBRDF.cpp  main.cpp  BrdfIO.cpp
HEADERS = IsotropicMERLBRDF.h MERLBRDFRead.h BrdfIO.h BRDF.h argstream.h Spectrum.h


