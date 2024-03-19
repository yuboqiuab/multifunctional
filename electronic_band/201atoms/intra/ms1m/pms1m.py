import numpy as np

with open('../inputms1m_1', 'r') as f1:
    row = f1.readlines()
a=len(row)
f2 = open('output', 'w');
for i in range(int(a/61)):
    f2.write("0.0000000   0.0000000   0.0000000\n")
