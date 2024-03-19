import numpy as np
import array
import math
import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
import time

f1 = open('all_output_new', 'w');
with open('all_output', 'r') as f2:
    row = f2.readlines()
for i in row:
   j=float(i)
   f1.write("%9.6f\n" % (j*100.0))
