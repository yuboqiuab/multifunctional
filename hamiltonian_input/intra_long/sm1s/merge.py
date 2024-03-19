import numpy as np

with open('reduce_inputsm1s', 'r') as f1:
    row1 = f1.readlines()
with open('../sm2m/reduce_inputsm2m', 'r') as f2:
    row2 = f2.readlines()


f3 = open('reduce_sm', 'w');
for i in range(len(row1)):
   k1=row1[i].split()
   for j1 in range(len(k1)):
      f3.write("%9.6f \n" % (float(k1[j1])))
   k2=row2[i].split()
   for j2 in range(len(k2)):
      f3.write("%9.6f \n" % (float(k2[j2])))
