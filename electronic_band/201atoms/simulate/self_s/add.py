import numpy as np
import array
import math

d1=-0.35818333333333324*1.0
d2=0.36479413159048546*1.0

#d1=-0.34031383
#d2=0.34034672

with open('see', 'r') as f1:
    row = f1.readlines()
i1,i2,i3 = [],[],[]
for i in row:
   j=i.split()
   i1.append(j[0])
   i2.append(j[1])
   i3.append(j[2])

i1 = np.array(i1).reshape(-1, 1)
i1 = i1.astype(np.float64)
i2 = np.array(i2).reshape(-1, 1)
i2 = i2.astype(np.float64)
i3 = np.array(i3).reshape(-1, 1)
i3 = i3.astype(np.float64)


o2 = open('pred_target.dat', 'w');
d=0
for i in range(len(i1)):
    if(int(i3[i])==1):
        i1[i]=i1[i]+d1
    if(int(i3[i])==2):
        i1[i]=i1[i]+d2
    o2.write("%9.6f %9.6f \n" % (float(i1[i]),float(i2[i])))
