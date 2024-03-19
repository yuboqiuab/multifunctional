import numpy as np
import array
import math
import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
import time

natom=67+134
nmo=int(natom/3)
ns=int(natom/3*2)
nw=nmo*5+ns*3
tw=nw*nw
lines=nw*nw*9

i1=np.zeros(lines,dtype='int')
i2=np.zeros(lines,dtype='int')
i3=np.zeros(lines,dtype='int')
i4=np.zeros(lines,dtype='int')
i5=np.zeros(lines,dtype='int')
i6=np.zeros(lines,dtype='float')
i7=np.zeros(lines,dtype='float')

for i in range(lines):
   if ((i>=0*tw) and (i<1*tw)):
      i1[i]=-1
      i2[i]=-1
      i3[i]= 0
   if ((i>=1*tw) and (i<2*tw)):
      i1[i]=-1
      i2[i]= 0
      i3[i]= 0
   if ((i>=2*tw) and (i<3*tw)):
      i1[i]=-1
      i2[i]= 1
      i3[i]= 0
   if ((i>=3*tw) and (i<4*tw)):
      i1[i]= 0
      i2[i]=-1
      i3[i]= 0
   if ((i>=4*tw) and (i<5*tw)):
      i1[i]= 0
      i2[i]= 0
      i3[i]= 0
   if ((i>=5*tw) and (i<6*tw)):
      i1[i]= 0
      i2[i]= 1
      i3[i]= 0
   if ((i>=6*tw) and (i<7*tw)):
      i1[i]= 1
      i2[i]=-1
      i3[i]= 0
   if ((i>=7*tw) and (i<8*tw)):
      i1[i]= 1
      i2[i]= 0
      i3[i]= 0
   if ((i>=8*tw) and (i<9*tw)):
      i1[i]= 1
      i2[i]= 1
      i3[i]= 0
   i4[i]=i%tw%nw+1;
   i5[i]=i%tw/nw+1;
   i6[i]=-0
   i7[i]=-0


with open('merge', 'r') as f1:
    row = f1.readlines()
l1,l2,l3 = [],[],[]
for i in row:
   j=i.split()
   l1.append(j[0])
   l2.append(j[1])
   l3.append(j[2])
l1 = np.array(l1).reshape(-1, 1)
l1 = l1.astype(np.float64)
l2 = np.array(l2).reshape(-1, 1)
l2 = l2.astype(np.float64)
l3 = np.array(l3).reshape(-1, 1)
l3 = l3.astype(np.float64)

for i in range(len(l1)):
    j=int(l1[i])
    i6[j]=l2[i]

o2 = open('sh2', 'w');
for i in range(lines):
   o2.write(" %4d %4d %4d %4d %4d   %9.6f   %9.6f\n" % (i1[i],i2[i],i3[i],i4[i],i5[i],i6[i],i7[i]))


# IMAGE_IDX = {(-1, -1, 0): 0, (-1, 0, 0): 1, (0, -1, 0): 2, (0, 0, 0): 3, (0, 1, 0): 4, (1, 0, 0): 5, (1, 1, 0): 6}
