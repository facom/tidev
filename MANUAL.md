tidev
=====

A general framework to calculate the evolution of rotation for tidally
interacting bodies using the formalism by Efroimsky.

Installing libconfig and gsl from sources
-----------------------------------------

If you are working on a non-debian environment, Mac or a distributed
environment we recommend to download and compile dependencies from the
sources.  Follow the next procedure to do it:

1. Untar the sources provided with the package in the util/src
   directory:

	util/src $ tar zxvf libconfig-1.4.9.tar.gz

	util/src $ tar zxvf gsl-1.16.tar.gz

2. Configure and compile:

      	util/src $ cd gsl-1.16 && ./configure --prefix=$(pwd) && make && make install

   	util/src $ cd libconfig-1.4.9 && ./configure --prefix=$(pwd) && make && make install

3. Copy binary library and header files:

        util/src $ cp gsl-1.16/lib/*.{a,la} ../lib

        util/src $ cp gsl-1.16/include/gsl ../include

        util/src $ cp libconfig-1.4.9/lib/*.{a,la} ../lib

        util/src $ cp libconfig-1.4.9/include/*.{h,h++} ../include

