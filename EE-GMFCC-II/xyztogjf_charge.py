import os
from glob import glob
#from multiprocessing import Pool
import os
import re

ELEMENTS = ['X',  # Ghost
    'H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 'K' , 'Ca',
    'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' , 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I' , 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
CHARGES = dict(((str(i),x) for i,x in enumerate(ELEMENTS)))


def checkcharge_monomer(aa):
    if aa in ['ARG', 'LYS']:
        return 1
    elif aa in ['GLU', 'ASP']:
        return -1
    else:
        return 0

def split_file_by_empty_lines(input_filename):
    index=[]
    for line in range(len(input_filename)):
        line1 = input_filename[line].strip()
        if not line1:
            index.append(line)
    return index

def ChooseTest(log):
    T=open(log)
    name=log.split("/")[-1].split(".")[0]
    charge=[]
    if "-" in name:
        mono1=name.split("-")[0].split("_")[-1]
        mono2=name.split("-")[1].split("_")[-1]
        charge=checkcharge_monomer(mono1)+checkcharge_monomer(mono2)
    else:
        mono1=name.split("-")[0].split("_")[-1]
        charge=checkcharge_monomer(mono1)
    lines = T.readlines()
    result=lines[2:]
    index=split_file_by_empty_lines(result)
    xyz = result[:index[0]]
    xyz = [i.strip() for i in xyz]
    xyz_charge=result[index[0]+1:index[1]]
    xyz_charge=[" ".join(i.strip().split()[:4]) for i in xyz_charge]
    xyz0=["%nproc=16","%mem=20GB","#p wb97xd/6-31g* nosymm charge scf(conver=7,xqc,MaxConventionalCycles=90)","  ","test"," "]
    xyz1=[str(charge)+" 1"]
    T=xyz0+xyz1+xyz+[" "]+xyz_charge+[" "]+[" "]
    return T


R=open("xyzfile.raw")
E=R.readlines()
for p in E:
    path=p.strip()
    name=p.split("/")[-1]
    inpu=ChooseTest(path)
    #print(inpu)
    f=open(path+".gjf","w")
    f.write("\n".join(inpu))

