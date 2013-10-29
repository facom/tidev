#######################################################################
#   _______     ________     
#  /_  __(_)___/ / ____/   __
#   / / / / __  / __/ | | / /
#  / / / / /_/ / /___ | |/ / 
# /_/ /_/\__,_/_____/ |___/  
# Tidal Spin Evolution
#######################################################################
# Copyright (C) 2013 
# Jorge Zuluaga (zuluagajorge@gmail.com, Mario Melita (melita@iafe.uba.ar)
# Pablo Cuartas (quarktas@gmail.com), Bayron Portilla (bayron@gmail.com)
#######################################################################
# MASTER MAKEFILE
#######################################################################
CC=g++
OPTIM=-O4
#OPTIM=-g
CFLAGS=$(OPTIM) -w -c -I. -Iutil/include $(OPTIONS)
LFLAGS=$(OPTIM) -lm util/lib/libgsl.a util/lib/libgslcblas.a util/lib/libconfig++.a
EDITOR=emacs -nw

%.out:%.o
	$(CC) $^ $(LFLAGS) -o $@

%.o:%.cpp tidev.cpp
	$(CC) $(CFLAGS) $< -o $@

run:
	./tidev-full.out 

cleanout:
	rm -rf *.out *.o *.exe

clean:cleanout
	rm -rf *.log *~ *.dump
	find . -name *~ -exec rm -rf {} \;

cleandump:
	rm -rf *.dump

cleanall:clean cleanout cleandump
	rm -rf *.dat

edit: 
	emacs -nw *.md TODO makefile* *.cpp *.cfg

commit:
	git commit -am "Commit"
	git push origin master
