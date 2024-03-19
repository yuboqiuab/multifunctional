import csv
import functools
import json, gzip
import torch, os
import zlib
import math
from pathlib import Path
import os.path as osp

from itertools import product

import numpy as np
from numpy.linalg import norm
from pymatgen.core.structure import Structure
from pymatgen.core.bonds import CovalentBond
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import torch.nn.functional as F
from torch.utils.data import Dataset
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader as GraphDataLoader
import matplotlib.pyplot as plt
from torch.utils.data.dataset import random_split, Subset

from multiprocessing import Pool
from tqdm import tqdm


PATH = Path(__file__).parent





class CrystalGraphDataset(Dataset):
    def __init__(self, dataset_dir, dataset_ratio=1.0, radius=8.0, n_nbr=18):
        super().__init__()

        self.radius = radius
        self.n_nbr = n_nbr
        self.struc_dir = osp.join(dataset_dir, 'structures')
        self.ham_dir = osp.join(dataset_dir, 'forces')
        self.ids = os.listdir(self.struc_dir)
        self.ids = self.ids[:int(len(self.ids) * dataset_ratio)]


    def __len__(self):
        return len(self.ids)

    @functools.lru_cache(maxsize=None)
    def __getitem__(self, index):

        file = osp.join(self.struc_dir, self.ids[index])
        fid = self.ids[index].split('_')[1]
        with open(file, 'r') as f:
            cont = f.read()
        struc = Structure.from_str(cont, fmt='poscar')
        f_file = osp.join(self.ham_dir, 'forces_{}'.format(fid))
        fx,fy,fz = [],[],[]
        with open(f_file, 'r') as f:
            for line in f:
                fx.append(float(line.split()[-3]))
                fy.append(float(line.split()[-2]))
                fz.append(float(line.split()[-1]))
        fx = np.array(fx)
        fy = np.array(fy)
        fz = np.array(fz)
        self._get_edge_and_face(struc,fx,fy,fz,fid)


    def _get_edge_and_face(self, struc, fx,fy,fz,fid):
        vect1=[]
        with open('part1', 'r') as f:
             for row in f:
                vect1.append(row.split())
        vect1 = (np.array(vect1)).astype(np.float64)

        vect2=[]
        with open('part2', 'r') as f:
             for row in f:
                vect2.append(row.split())
        vect2 = (np.array(vect2)).astype(np.float64)

        vect3=[]
        with open('part3', 'r') as f:
             for row in f:
                vect3.append(row.split())
        vect3 = (np.array(vect3)).astype(np.float64)

        vect4=-vect1
        vect5=-vect2
        vect6=-vect3

        m=struc.lattice.matrix
        lattice=Lattice(struc.lattice.matrix)
        num=struc.num_sites
        nw=int(num/3*5+num/3*2*3)
        atom=[]
        ave1=0
        for i in range(num):
            atom.append(struc.frac_coords[i])
            ave1=ave1+struc.frac_coords[i][2]
        ave1=ave1/(num*1.0)
        ave2=0
        ave3=0
        num2=0
        num3=0
        for i in range(num):
            if(atom[i][2]<ave1):
                num2=num2+1
                ave2=ave2+atom[i][2]
            if(atom[i][2]>ave1):
                num3=num3+1
                ave3=ave3+atom[i][2]
        ave2=ave2/(num2*1.0)
        ave3=ave3/(num3*1.0)

        c,d,bc = [],[],[]
        for i in range(num):
            if((i<int(num/3)) and (atom[i][2])<ave1):
               a,b = self._get_list1(struc, vect1, atom[i])
            if((i<int(num/3)) and (atom[i][2])>ave1):
               a,b = self._get_list1(struc, vect4, atom[i])
            if((i>=int(num/3)) and (atom[i][2])<ave1 and (atom[i][2])>ave2):
               a,b = self._get_list2(struc, vect2, atom[i])
            if((i>=int(num/3)) and (atom[i][2])<ave1 and (atom[i][2])<ave2):
               a,b = self._get_list2(struc, vect3, atom[i])
            if((i>=int(num/3)) and (atom[i][2])>ave1 and (atom[i][2])>ave3):
               a,b = self._get_list2(struc, vect6, atom[i])
            if((i>=int(num/3)) and (atom[i][2])>ave1 and (atom[i][2])<ave3):
               a,b = self._get_list2(struc, vect5, atom[i])
            c.append(a)
            bc.append(a)
            d.append(b)


        for i in range(num):
            if((i<int(num/3)) and (atom[i][2])<ave1):
                temp_r = 1000.0
                for j in range(int(num/3)):
                    if(atom[j][2]>ave1):
                        a,b = lattice.get_distance_and_image(atom[i],atom[j])
                    else:
                        a = 1000.0
                    if(a<temp_r):
                        temp_r=a
                        temp_image=b
                        temp_j=j
                a,b = self._get_list1(struc, vect4, atom[temp_j])
                dis=lattice.get_cartesian_coords(atom[temp_j]+temp_image-atom[i])
                add=np.zeros(4,dtype='float')
                new=np.zeros(4,dtype='float')
                rad=dis[0]*dis[0]+dis[1]*dis[1]+dis[2]*dis[2];
                add[0]=1.0/np.sqrt(rad)
                add[1]=dis[0]/rad
                add[2]=dis[1]/rad
                add[3]=dis[2]/rad
                bc[i].append(add)
                for k in range(len(a)):
                    add=np.zeros(4,dtype='float')
                    new=np.zeros(4,dtype='float')
                    add[0]=(a[k][1]/a[k][0]/a[k][0]+dis[0])
                    add[1]=(a[k][2]/a[k][0]/a[k][0]+dis[1])
                    add[2]=(a[k][3]/a[k][0]/a[k][0]+dis[2])
                    rad=add[0]*add[0]+add[1]*add[1]+add[2]*add[2];
                    new[0]=1.0/np.sqrt(rad)
                    new[1]=add[0]/rad
                    new[2]=add[1]/rad
                    new[3]=add[2]/rad
                    bc[i].append(new)

            if((i<int(num/3)) and (atom[i][2])>ave1):
                temp_r = 1000.0
                for j in range(int(num/3)):
                    if(atom[j][2]<ave1):
                        a,b = lattice.get_distance_and_image(atom[i],atom[j])
                    else:
                        a = 1000.0
                    if(a<temp_r):
                        temp_r=a
                        temp_image=b
                        temp_j=j
                a,b = self._get_list1(struc, vect1, atom[temp_j])
                dis=lattice.get_cartesian_coords(atom[temp_j]+temp_image-atom[i])
                add=np.zeros(4,dtype='float')
                new=np.zeros(4,dtype='float')
                rad=dis[0]*dis[0]+dis[1]*dis[1]+dis[2]*dis[2];
                add[0]=1.0/np.sqrt(rad)
                add[1]=dis[0]/rad
                add[2]=dis[1]/rad
                add[3]=dis[2]/rad
                bc[i].append(add)
                for k in range(len(a)):
                    add=np.zeros(4,dtype='float')
                    new=np.zeros(4,dtype='float')
                    add[0]=(a[k][1]/a[k][0]/a[k][0]+dis[0])
                    add[1]=(a[k][2]/a[k][0]/a[k][0]+dis[1])
                    add[2]=(a[k][3]/a[k][0]/a[k][0]+dis[2])
                    rad=add[0]*add[0]+add[1]*add[1]+add[2]*add[2];
                    new[0]=1.0/np.sqrt(rad)
                    new[1]=add[0]/rad
                    new[2]=add[1]/rad
                    new[3]=add[2]/rad
                    bc[i].append(new)


            if((i>=int(num/3)) and (atom[i][2])<ave1 and (atom[i][2])>ave2):
                temp_r = 1000.0
                for j in range(int(num/3),num):
                    if(atom[j][2]>ave1 and atom[j][2]<ave3):
                        a,b = lattice.get_distance_and_image(atom[i],atom[j])
                    else:
                        a = 1000.0
                    if(a<temp_r):
                        temp_r=a
                        temp_image=b
                        temp_j=j
                a,b = self._get_list2(struc, vect5, atom[temp_j])
                dis=lattice.get_cartesian_coords(atom[temp_j]+temp_image-atom[i])
                add=np.zeros(4,dtype='float')
                new=np.zeros(4,dtype='float')
                rad=dis[0]*dis[0]+dis[1]*dis[1]+dis[2]*dis[2];
                add[0]=1.0/np.sqrt(rad)
                add[1]=dis[0]/rad
                add[2]=dis[1]/rad
                add[3]=dis[2]/rad
                bc[i].append(add)
                for k in range(len(a)):
                    add=np.zeros(4,dtype='float')
                    new=np.zeros(4,dtype='float')
                    add[0]=(a[k][1]/a[k][0]/a[k][0]+dis[0])
                    add[1]=(a[k][2]/a[k][0]/a[k][0]+dis[1])
                    add[2]=(a[k][3]/a[k][0]/a[k][0]+dis[2])
                    rad=add[0]*add[0]+add[1]*add[1]+add[2]*add[2];
                    new[0]=1.0/np.sqrt(rad)
                    new[1]=add[0]/rad
                    new[2]=add[1]/rad
                    new[3]=add[2]/rad
                    bc[i].append(new)
     #           for k in range(5):
     #               new=np.zeros(4,dtype='float')
     #               bc[i].append(new)

            if((i>=int(num/3)) and (atom[i][2])<ave1 and (atom[i][2])<ave2):
                temp_r = 1000.0
                for j in range(int(num/3),num):
                    if(atom[j][2]>ave1 and atom[j][2]<ave3):
                        a,b = lattice.get_distance_and_image(atom[i],atom[j])
                    else:
                        a = 1000.0
                    if(a<temp_r):
                        temp_r=a
                        temp_image=b
                        temp_j=j
                a,b = self._get_list2(struc, vect5, atom[temp_j])
                dis=lattice.get_cartesian_coords(atom[temp_j]+temp_image-atom[i])
                add=np.zeros(4,dtype='float')
                new=np.zeros(4,dtype='float')
                rad=dis[0]*dis[0]+dis[1]*dis[1]+dis[2]*dis[2];
                add[0]=1.0/np.sqrt(rad)
                add[1]=dis[0]/rad
                add[2]=dis[1]/rad
                add[3]=dis[2]/rad
                bc[i].append(add)
                for k in range(len(a)):
                    add=np.zeros(4,dtype='float')
                    new=np.zeros(4,dtype='float')
                    add[0]=(a[k][1]/a[k][0]/a[k][0]+dis[0])
                    add[1]=(a[k][2]/a[k][0]/a[k][0]+dis[1])
                    add[2]=(a[k][3]/a[k][0]/a[k][0]+dis[2])
                    rad=add[0]*add[0]+add[1]*add[1]+add[2]*add[2];
                    new[0]=1.0/np.sqrt(rad)
                    new[1]=add[0]/rad
                    new[2]=add[1]/rad
                    new[3]=add[2]/rad
                    bc[i].append(new)
    #            for k in range(5):
    #                new=np.zeros(4,dtype='float')
    #                bc[i].append(new)


            if((i>=int(num/3)) and (atom[i][2])>ave1 and (atom[i][2])>ave3):
                temp_r = 1000.0
                for j in range(int(num/3),num):
                    if(atom[j][2]<ave1 and atom[j][2]>ave2):
                        a,b = lattice.get_distance_and_image(atom[i],atom[j])
                    else:
                        a = 1000.0
                    if(a<temp_r):
                        temp_r=a
                        temp_image=b
                        temp_j=j
                a,b = self._get_list2(struc, vect2, atom[temp_j])
                dis=lattice.get_cartesian_coords(atom[temp_j]+temp_image-atom[i])
                add=np.zeros(4,dtype='float')
                new=np.zeros(4,dtype='float')
                rad=dis[0]*dis[0]+dis[1]*dis[1]+dis[2]*dis[2];
                add[0]=1.0/np.sqrt(rad)
                add[1]=dis[0]/rad
                add[2]=dis[1]/rad
                add[3]=dis[2]/rad
                bc[i].append(add)
                for k in range(len(a)):
                    add=np.zeros(4,dtype='float')
                    new=np.zeros(4,dtype='float')
                    add[0]=(a[k][1]/a[k][0]/a[k][0]+dis[0])
                    add[1]=(a[k][2]/a[k][0]/a[k][0]+dis[1])
                    add[2]=(a[k][3]/a[k][0]/a[k][0]+dis[2])
                    rad=add[0]*add[0]+add[1]*add[1]+add[2]*add[2];
                    new[0]=1.0/np.sqrt(rad)
                    new[1]=add[0]/rad
                    new[2]=add[1]/rad
                    new[3]=add[2]/rad
                    bc[i].append(new)
      #          for k in range(5):
      #              new=np.zeros(4,dtype='float')
      #              bc[i].append(new)

            if((i>=int(num/3)) and (atom[i][2])>ave1 and (atom[i][2])<ave3):
                temp_r = 1000.0
                for j in range(int(num/3),num):
                    if(atom[j][2]<ave1 and atom[j][2]>ave2):
                        a,b = lattice.get_distance_and_image(atom[i],atom[j])
                    else:
                        a = 1000.0
                    if(a<temp_r):
                        temp_r=a
                        temp_image=b
                        temp_j=j
                a,b = self._get_list2(struc, vect2, atom[temp_j])
                dis=lattice.get_cartesian_coords(atom[temp_j]+temp_image-atom[i])
                add=np.zeros(4,dtype='float')
                new=np.zeros(4,dtype='float')
                rad=dis[0]*dis[0]+dis[1]*dis[1]+dis[2]*dis[2];
                add[0]=1.0/np.sqrt(rad)
                add[1]=dis[0]/rad
                add[2]=dis[1]/rad
                add[3]=dis[2]/rad
                bc[i].append(add)
                for k in range(len(a)):
                    add=np.zeros(4,dtype='float')
                    new=np.zeros(4,dtype='float')
                    add[0]=(a[k][1]/a[k][0]/a[k][0]+dis[0])
                    add[1]=(a[k][2]/a[k][0]/a[k][0]+dis[1])
                    add[2]=(a[k][3]/a[k][0]/a[k][0]+dis[2])
                    rad=add[0]*add[0]+add[1]*add[1]+add[2]*add[2];
                    new[0]=1.0/np.sqrt(rad)
                    new[1]=add[0]/rad
                    new[2]=add[1]/rad
                    new[3]=add[2]/rad
                    bc[i].append(new)
   #             for k in range(5):
   #                 new=np.zeros(4,dtype='float')
   #                 bc[i].append(new)

        f3 = open('input_'+fid, 'w');
        f4 = open('output_'+fid, 'w');

   #    print(d[num-1])
   #     for i in range(len(c[33])):
   #         f2.write("%12.10f  %12.10f  %12.10f  %12.10f\n" % (c[33][i][0],c[33][i][1],c[33][i][2],c[33][i][3]))

        for ai in range(int(num/3),int(num)):
            mu=100.0
            zero=0.0
            input=[]
            for ti in range(len(bc[ai])):
                input.append(bc[ai][ti].tolist())
            for ti in range(len(input)):
                f3.write("%12.10f   %12.10f   %12.10f   %12.10f\n" % (input[ti][0],input[ti][1],input[ti][2],input[ti][3]))
            f4.write("%12.10f   %12.10f   %12.10f\n" % (fx[ai]*mu,fy[ai]*mu,fz[ai]*mu))


        f3.close()
        f4.close()

    def _get_list1(self, struc, vect, ctr):
        nf=[]
        nv=[]
        nf.append(ctr)
        m=struc.lattice.matrix
        lattice=Lattice(struc.lattice.matrix)
        num=struc.num_sites
        atom=[]
        for i in range(num):
            atom.append(struc.frac_coords[i])
        for i1 in range(len(vect)):
           fv=lattice.get_fractional_coords(vect[i1])
           if(i1<6):
               new_ctr=ctr
               element='S'
           if((i1>=6) and (i1<8)):
               new_ctr=nf[1]
               element='Mo'
           if((i1>=8) and (i1<10)):
               new_ctr=nf[3]
               element='Mo'
           if((i1>=10) and (i1<12)):
               new_ctr=nf[5]
               element='Mo'
           if((i1>=12) and (i1<14)):
               new_ctr=nf[7]
               element='S'
           if((i1>=14) and (i1<16)):
               new_ctr=nf[8]
               element='S'
           if((i1>=16) and (i1<18)):
               new_ctr=nf[8]
               element='S'
           if((i1>=18) and (i1<20)):
               new_ctr=nf[9]
               element='S'
           if((i1>=20) and (i1<22)):
               new_ctr=nf[10]
               element='S'
           if((i1>=22) and (i1<24)):
               new_ctr=nf[10]
               element='S'
           if((i1>=24) and (i1<26)):
               new_ctr=nf[11]
               element='S'
           if((i1>=26) and (i1<28)):
               new_ctr=nf[12]
               element='S'
           if((i1>=28) and (i1<30)):
               new_ctr=nf[12]
               element='S'

           fatom = new_ctr + fv
           t_r=1000.0
           for i2 in range(num):
               if(struc[i2].specie.name == element):
                   dis, image = lattice.get_distance_and_image(atom[i2],fatom)
                   if(dis<t_r):
                     t_r=dis
                     t_i=i2
                     t_image=image
           a=atom[t_i]-t_image
           nf.append(a)
           a=a-ctr
           a=lattice.get_cartesian_coords(a)
           rad=a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
           b= np.zeros(4,dtype='float')
           b[0]=1.0/np.sqrt(rad)
           b[1]=a[0]/rad
           b[2]=a[1]/rad
           b[3]=a[2]/rad
           nv.append(b)
        return nv, nf

    def _get_list2(self, struc, vect, ctr):
        nf=[]
        nv=[]
        nf.append(ctr)
        m=struc.lattice.matrix
        lattice=Lattice(struc.lattice.matrix)
        num=struc.num_sites
        atom=[]
        for i in range(num):
            atom.append(struc.frac_coords[i])
        for i1 in range(len(vect)):
           fv=lattice.get_fractional_coords(vect[i1])
           if(i1==0):
               new_ctr=ctr
               element='S'
           if((i1>=1) and (i1<4)):
               new_ctr=ctr
               element='Mo'
           if((i1>=4) and (i1<8)):
               new_ctr=nf[2]
               element='S'
           if((i1>=8) and (i1<12)):
               new_ctr=nf[3]
               element='S'
           if((i1>=12) and (i1<16)):
               new_ctr=nf[4]
               element='S'
           if((i1>=16) and (i1<17)):
               new_ctr=nf[5]
               element='Mo'
           if((i1>=17) and (i1<19)):
               new_ctr=nf[7]
               element='Mo'
           if((i1>=19) and (i1<20)):
               new_ctr=nf[9]
               element='Mo'
           if((i1>=20) and (i1<22)):
               new_ctr=nf[11]
               element='Mo'
           if((i1>=22) and (i1<23)):
               new_ctr=nf[13]
               element='Mo'
           if((i1>=23) and (i1<25)):
               new_ctr=nf[15]
               element='Mo'
           fatom = new_ctr + fv
           t_r=1000.0
           for i2 in range(num):
               if(struc[i2].specie.name == element):
                   dis, image = lattice.get_distance_and_image(atom[i2],fatom)
                   if(dis<t_r):
                     t_r=dis
                     t_i=i2
                     t_image=image
           a=atom[t_i]-t_image
           nf.append(a)
           a=a-ctr
           a=lattice.get_cartesian_coords(a) 
           rad=a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
           b= np.zeros(4,dtype='float')
           b[0]=1.0/np.sqrt(rad)
           b[1]=a[0]/rad
           b[2]=a[1]/rad
           b[3]=a[2]/rad
           nv.append(b)
        b= np.zeros(4,dtype='float')
      #  nv.append(b)
      #  nv.append(b)
      #  nv.append(b)
      #  nv.append(b)
      #  nv.append(b)
        return nv, nf

def get_orb_idx(i):
    if i < 32:
        return list(range(i * 5, i * 5 + 5))
    else:
        return list(range(160 + (i - 32) * 3, 160 + (i - 32) * 3 + 3))


IMAGE_IDX = {(-2, -2, 0): 0, (-2, -1, 0): 1, (-2, 0, 0): 2, (-1, -2, 0): 3, (-1, -1, 0): 4, (-1, 0, 0): 5, (-1, 1, 0): 6, 
             ( 0, -2, 0): 7, ( 0, -1, 0): 8, ( 0, 0, 0): 9, ( 0,  1, 0): 10, (0, 2, 0): 11, (1, -1, 0): 12, (1, 0, 0): 13,
             ( 1,  1, 0): 14, ( 1, 2, 0): 15, ( 2, 0, 0): 16, ( 2,  1, 0): 17, (2, 2, 0): 18}


def get_orb_rows(ctr_idx, nbr_idx, image):
    indices = np.zeros(64)
    mask = np.zeros(64)
    if image in IMAGE_IDX:
        image_idx = IMAGE_IDX[image] * 352 * 352
        ctr_indices = get_orb_idx(ctr_idx)
        nbr_indices = get_orb_idx(nbr_idx)
        if len(ctr_indices) == 3:
            if len(nbr_indices) == 3:
                start = 0
            else:
                start = 9
        else:
            if len(nbr_indices) == 3:
                start = 24
            else:
                start = 39
        for i, (c, n) in enumerate(product(ctr_indices, nbr_indices)):
            idx = image_idx + c + n * 352
            indices[start + i] = idx
            mask[start + i] = 1

    return indices, mask



if __name__ == '__main__':
    didx = '96'
    dataset = CrystalGraphDataset('../dataset/MoS2_{}'.format(didx))
#    a=math.acos(0.5)
#    print(a)
    def process_data(idx):
        dataset[idx]
    N = len(dataset)
    pool = Pool(processes=20)

    for _ in tqdm(pool.imap_unordered(process_data, range(N)), total=N):
        pass
