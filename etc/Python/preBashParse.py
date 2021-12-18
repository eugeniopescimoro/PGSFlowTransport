#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 09:09:44 2021
@author: pmxep5
It instructs bash to grep values from an OpenFOAM log file following word patterns
"""
#!/usr/bin/env python3
#import numpy as np
import os
import subprocess

sim = 1
FS = 5
os.chdir("/data/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_3")
for i in range(0, sim):
    b = i+FS
    os.chdir("TS%d" % b)
    current_dir = os.getcwd()
    os.makedirs("./LOGs", exist_ok = True)
    subprocess.run(['/bin/bash', '-c', 'cat log | grep \'Adaptive time =\' | cut -d\' \' -f4 > LOGs/logTime'])
    subprocess.run(['/bin/bash', '-c', 'cat log | grep \'Total mass =\' | cut -d\' \' -f4 > LOGs/logMass'])
    subprocess.run(['/bin/bash', '-c', 'cat log | grep \'Mean vel =\' | cut -d\' \' -f4 | tr -d \'(\' > LOGs/logVelx'])
    subprocess.run(['/bin/bash', '-c', 'cat log | grep \'Flux out =\' | cut -d\' \' -f4 > LOGs/logFlux'])    
    os.chdir("..")