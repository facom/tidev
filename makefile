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
CFLAGS=$(OPTIM) -w -c -I. $(OPTIONS)
LFLAGS=$(OPTIM) -lm -lgsl -lgslcblas -lconfig++
EDITOR=emacs -nw

%.out:%.o
	$(CC) $(LFLAGS) $^ -o $@

%.o:%.cpp tidev.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o *.exe *.out *.log *~ \#*
	find . -name *~ -exec rm -rf {} \;
	rm -rf *.png

cleanout:
	rm -rf *.out *.o

cleanall:clean
	rm -rf *.dat

edit: 
	emacs -nw *.md TODO makefile* *.cpp *.cfg

commit:
	git commit -am "Commit"
	git push
