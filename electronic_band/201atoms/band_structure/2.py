import numpy as np
import matplotlib.pyplot as plt
import time



time1 = time.perf_counter();
num_wan = 737
k_mesh = 50
E_fermi =  -0.3783
lr=7;
hr=3;

basis_vector=[]
with open('lattice', 'r') as f1:
  for row in f1:
     basis_vector.append(row.split()[0:len(row)])
basis_vector = np.array(basis_vector).reshape(-1, 3)
basis_vector = basis_vector.astype(np.float64)

K_point_path=[]
with open('kpath', 'r') as f2:
  for row in f2:
     K_point_path.append(row.split()[0:len(row)])
K_point_path = np.array(K_point_path).reshape(-1, 3)
K_point_path = K_point_path.astype(np.float64)



V = np.dot(basis_vector[0], np.cross(basis_vector[1], basis_vector[2]) )
rec = [np.cross(basis_vector[1], basis_vector[2]) * 2 * np.pi / V,
       np.cross(basis_vector[2], basis_vector[0]) * 2 * np.pi / V,
       np.cross(basis_vector[0], basis_vector[1]) * 2 * np.pi / V]

k_point=[]
for i in range(len(K_point_path)):
    K_point_path[i] = K_point_path[i][0] * rec[0] + K_point_path[i][1] * rec[1] + K_point_path[i][2] * rec[2]
    if(i!=0):
       interval = np.array(K_point_path[i]) - np.array(K_point_path[i-1])
       interval = interval / k_mesh
       for k in range(k_mesh):
          a=K_point_path[i-1] + k * interval
          k_point.append(a)


lk=len(k_point);
phase = np.zeros((lr, lr, lr, lk),dtype='complex')
for r1 in range(lr):
    for r2 in range(lr):
        for r3 in range(lr):
            for k1 in range(lk):
                R1_vector = (r1-hr) * np.array(basis_vector[0])
                R2_vector = (r2-hr) * np.array(basis_vector[1])
                R3_vector = (r3-hr) * np.array(basis_vector[2])
                R_vec = R1_vector + R2_vector + R3_vector
                inner_product = np.dot(R_vec, k_point[k1])
                phase[r1][r2][r3][k1] = np.exp(1j * inner_product)


time2 = time.perf_counter();
print("Time for Initializing := ", time2 - time1)

lc=[]
with open('sh2', 'r') as f:
  for row in f:
     lc.append(row.split()[0:7])
lc = np.array(lc).reshape(-1,7)
ll = len(lc);

time3 = time.perf_counter();
print("Time for Reading := ", time3 - time2)

H = np.zeros((lk, num_wan, num_wan),dtype='complex')
for l in range(ll):
    for k1 in range(lk):
       i0=int(lc[l][0])+hr;
       i1=int(lc[l][1])+hr;
       i2=int(lc[l][2])+hr;
       i3=int(lc[l][3])-1;
       i4=int(lc[l][4])-1;
       f5=float(lc[l][5]);
       f6=float(lc[l][6]);
       H[k1][i3][i4] = H[k1][i3][i4] + (f5+1j*f6)*phase[i0][i1][i2][k1];

time4 = time.perf_counter();
print("Time for Building H := ", time4 - time3)

ke = np.zeros((lk, num_wan),dtype='complex')
for k1 in range(lk):
    eig = np.linalg.eigvals(H[k1])
    idx = np.argsort(eig)
    eig = eig[idx]
    ke[k1]=eig


time5 = time.perf_counter();
print("Time for Solving H := ", time5 - time4)

f3 = open('out.dat', 'w');
for i in range(num_wan):
    for k1 in range(lk):
        a=np.real(ke[k1][i])-E_fermi
        f3.write("%12.10f\n" % a)
    a=np.real(ke[0][i])-E_fermi
    f3.write("%12.10f\n" % a)
    f3.write("\n")


time6 = time.perf_counter();
print("Total time count is := ", time6 - time1)
