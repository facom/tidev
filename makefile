#######################################################################
#   _______     ________     
#  /_  __(_)___/ / ____/   __
#   / / / / __  / __/ | | / /
#  / / / / /_/ / /___ | |/ / 
# /_/ /_/\__,_/_____/ |___/  
# Tidal Spin Evolution
#######################################################################
# Copyright (C) 2012 
# Mario Melita (melita@iafe.uba.ar), Pablo Cuartas (quarktas@gmail.com)
# Jorge Zuluaga (zuluagajorge@gmail.com
#######################################################################
# MASTER MAKEFILE
#######################################################################
CC=g++
OPTIM=-O4
#OPTIM=-g
CFLAGS=$(OPTIM) -c -I. $(OPTIONS)
LFLAGS=$(OPTIM) -lm -lgsl -lgslcblas -lconfig++

%.out:%.o
	$(CC) $(LFLAGS) $^ -o $@

%.o:%.cpp tidev.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o *.exe *.out *~
	find . -name *~ -exec rm -rf {} \;

cleanall:clean
	rm -rf *.dat


#CFLAGS=$(OPTIM) -c -I. -Iutil/include $(OPTIONS)
#LFLAGS=$(OPTIM) -lm util/lib/libgsl.a util/lib/libgslcblas.a -lconfig++
