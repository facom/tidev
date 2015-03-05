from os import system
import numpy as np

#READ RUN CONFIGURATION
conf=dict()
execfile("runs.cfg",{},conf)

#RANGE
param_vec=np.linspace(conf["param_ini"],conf["param_end"],
                      conf["nrun"])

#GENERATE RUN FILES
system("rm -r runs/run*")
for i in xrange(conf["nrun"]):
    paramval=param_vec[i]

    suffix="%03d"%i
    system("rm -rf runs/run_%s"%suffix)
    system("cp -rf runs/template runs/run_%s"%suffix)
    
    #CREATE CONFIGURATION FILE
    filecfg=open("runs/run_%s/tidev.cfg"%suffix,"w")
    fc=open("tidev.cfg")
    for line in fc:
        exec("cond='%s' in line"%conf["param"])
        if cond:
            line="%s = %.17e;"%(conf["param"],paramval)
        filecfg.write(line)
    filecfg.close()

    #CREATE LAUNCHING SCRIPT

