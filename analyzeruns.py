import os,commands
import numpy as np
system=os.system
System=commands.getoutput

#READ RUN CONFIGURATION
conf=dict()
execfile("runs.cfg",{},conf)

#LAUNCH RUNS
ne=0
ress=dict()
for i in xrange(conf["nrun"]):
    suffix="%03d"%i
    print "Analyzing run %d..."%(i+1)
    if os.path.lexists("runs/run_%s/.end"%suffix):
        ne+=1
        print "\tEnded."
        resonance=System("tail -n 1 runs/run_%s/*.dat |cut -f 8 -d ' '"%suffix)
        resval=np.round(float(resonance),1)
        resstr="%.1f"%resval
        if resstr in ress.keys():
            ress[resstr]+=1
        else:
            ress[resstr]=1
totruns=ne

for res in ress.keys():
    ress[res]/=(1.0*totruns);

print ress
        
