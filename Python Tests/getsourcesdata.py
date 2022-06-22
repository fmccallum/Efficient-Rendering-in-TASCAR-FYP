from cpu_tools import cmd_time, spk_circle,sources_point,sources_circle,logspace
import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle


os.system("export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH")
cmd = ["tascar_renderfile","-o outfile.wav","sources.tsc"]
plugin=sys.argv[1]
fname = sys.argv[2]
arrangement = sys.argv[3]
smode =0
if arrangement == "circle":
	smode = 1
	

#source_n = list(range(0,50,10))
source_n = logspace(0,3,7,zero=False)
repeats =[2*i for i in range(len(source_n),0,-1)]

#source_n.append(150)
#source_n.append(180)
#source_n.append(220)

spk_circle(46)
results = np.zeros((len(source_n),repeats[0]))

for n in range(len(source_n)):
	print(n)
	if smode == 1:
		sources_circle(source_n[n],plugin)
	else:
		sources_point(source_n[n],plugin)
	for r in range(repeats[n]):
		print("Repeat "+str(r+1)+ " out of "+str(repeats[n]))
		use = cmd_time(cmd)
		results[n,r]= use

			
		
avg =np.sum(results,axis=1)
avg/=repeats

suffix = "_p"
if smode == 1:
	suffix = "_c"
with open("data/"+fname+suffix, 'wb') as f:
	pickle.dump([avg,source_n], f)
