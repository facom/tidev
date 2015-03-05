import os,commands,time
import numpy as np
system=os.system
System=commands.getoutput

#READ RUN CONFIGURATION
conf=dict()
execfile("runs.cfg",{},conf)
verbose=False

#LAUNCH RUNS
nrun=conf["nrun"]
while True:

    for i in xrange(nrun):
        suffix="%03d"%i

        if verbose:print "Checking run %d..."%i
        
        #CHECK IF RUN IS RUNNIG
        if os.path.lexists("runs/run_%s/.start"%suffix):
            if verbose:print "\tRunning..."
            continue

        #CHECK NUMBER OF RUNNING PROCESSES
        nproc=System("at -l |wc -l")
        if verbose:print "Number of running processes %d..."%int(nproc)
        
        if int(nproc)>conf["runsim"]:
            if verbose:print "\tMaximum number of processes running..."
            break
        
        print "Launching run %d..."%(i+1)
        system("cd runs/run_%s;at now < run.sh"%suffix)

    print "%d process checked..."%i
    
    if i==nrun-1:break

    print "Sleeping..."
    time.sleep(1)
    
print "Done."

