# BRDF Conversion tool

Use this tool to convert BRDFs from/to MERL binary format to theta/phi format (a.k.a. titopo):

      ./MERLBinaryToTitopo -i gold-paint.binary -o gold-paint.titopo

Compilation:
      qmake
      make

QMake is a tool to generate Makefiles, provided by Qt. On Ubuntu you get qmake by installing the package qt5-qmake.

The code should compile/run on all plateforms but has only been tested on Ubuntu so far.
