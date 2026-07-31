import os
import glob 
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--dir", type=str, help="pdbfile")
args = parser.parse_args()

capped=glob.glob(args.dir+"/charge_cappednew/*log")
capped = sorted(capped,key=lambda x:int(x.split("/")[-1].split("_")[0]))
capped = capped[2:-2]

capped_xyz=glob.glob(args.dir+"/charge_cappednew/*xyz")
capped_xyz = sorted(capped_xyz,key=lambda x:int(x.split("/")[-1].split("_")[0]))
capped_xyz1 = capped_xyz[2:-2]
#capped_xyz = capped_xyz[1:-3]

pair=glob.glob(args.dir+"/charge_pairamino/*log")
pair= sorted(pair,key=lambda x:int(x.split("/")[-1].split("_")[0]))

pair_xyz = glob.glob(args.dir+"/charge_pairamino/*xyz")
pair_xyz = sorted(pair_xyz,key=lambda x:int(x.split("/")[-1].split("_")[0]))

charge_txt=np.loadtxt(args.dir+"/charge.txt")

def count_duplicates(two_d_list):
    count_dict = {}
    for sublist in two_d_list:
        # compute the number of same iteraction
        sublist_tuple = tuple(sublist)
        if sublist_tuple in count_dict:
            count_dict[sublist_tuple] += 1
        else:
            count_dict[sublist_tuple] = 1
    
    return count_dict

def merge_dicts(dictA, dictB):
    dictC = {}
    # merge interaction
    for key in dictA:
        if key in dictB:
            dictC[key] = dictA[key] + dictB[key]
        else:
            dictC[key] = dictA[key]
    for key in dictB:
        if key not in dictA:
            dictC[key] = dictB[key]
    return dictC

def same_dicts(dictA, dictB):
    dictC = {}
    # 
    for key in dictA:
        if key in dictB:
            # minus the interaction number cap*Acap and cap*cap
            dictC[key] = dictA[key]
    return dictC

def filter_dict_by_value(original_dict, target_value):
    # filter the rest interaction
    filtered_dict = {key: value for key, value in original_dict.items() if value == target_value}
    return filtered_dict


def get_interaction_pair(xyz):
    capped_xyz=open(xyz).readlines()[2:]
    index1=split_file_by_empty_lines(capped_xyz)
    capped_xyz_charge=capped_xyz[index1[1]+1:]
    bc_xyz_charge=capped_xyz[index1[0]+1:index1[1]]
    #print(len(capped_xyz_charge),len(bc_xyz_charge))
    pair=[]
    for h in capped_xyz_charge:
        for j in bc_xyz_charge:
            inter=[int(h.split()[-1]),int(j.split()[-1])]
            inter.sort() # sort the interaction
            #inter=check(inter)
            #print(inter)
            pair.append(inter)
    return pair

def get_energy(log):
    log=open(log).readlines()
    ene=[]
    for i in log:
        if "SCF Done" in i:
            energy=i.split()[4]
            ene.append(energy)
        elif "Self energy of the charges =" in i:
            energy=i.split()[6]
            ene.append(energy)
        elif "Nuclei-charges interaction =" in i:
            energy=i.split()[3]
            ene.append(energy)
    ene_tot,self_ene,ene_nc=ene[-1],ene[0],ene[1]
    return float(ene_tot),float(self_ene),float(ene_nc)

def split_file_by_empty_lines(input_filename):
    index=[]
    for line in range(len(input_filename)):
        line1 = input_filename[line].strip()
        if not line1:
            index.append(line)
    return index

def compute_ele(m,n):
    ene=[]
    for i in m:
        for j in n:
            Rij=np.linalg.norm(i[:3] - j[:3])
            ele_ene_unit = 333.05*((i[3] * j[3])/Rij)
            ene.append(ele_ene_unit)
    return sum(ene)

def compute_ele_two(m,n):
    Rij=np.linalg.norm(m[:3] - n[:3])
    ele_ene_unit = 333.05*((m[3] * n[3])/Rij)
    return ele_ene_unit


def get_charge_electron(capped_xyz, pair_xyz):
    capped_xyz=open(capped_xyz).readlines()[2:]
    index1=split_file_by_empty_lines(capped_xyz)
    capped_xyz_charge=capped_xyz[index1[1]+1:]
    capped_xyz_index=[int(i.strip().split()[4]) for i in capped_xyz_charge]
    capped_xyz_charge=[[float(i.split()[0]),float(i.split()[1]),float(i.split()[2]),float(i.split()[3])] for i in capped_xyz_charge]
    charge_m=np.array(capped_xyz_charge)
    pair_xyz=open(pair_xyz).readlines()[2:]
    index2=split_file_by_empty_lines(pair_xyz)
    pair_xyz_charge=pair_xyz[index2[1]+1:]
    max_index=max([int(ii.strip().split()[-1]) for ii in pair_xyz_charge])
    charge_n = charge_txt[max_index:]
    ele=compute_ele(charge_m,charge_n)
    return ele


def split_file_by_empty_lines(input_filename):
    index=[]
    for line in range(len(input_filename)):  
        line1 = input_filename[line].strip()
        if not line1:
            index.append(line)
    return index

def round_z(number, decimals):
    return f"{number:.{decimals}f}"

cap_ene=[get_energy(i)[0]-get_energy(i)[1] for i in capped]
pair_ene=[get_energy(i)[0]-get_energy(i)[1] for i in pair]

capped_inter_num=[]
for jj in range(len(capped_xyz1)):
    num=get_interaction_pair(capped_xyz1[jj])
    for kk in num:
        capped_inter_num.append(kk)
cap_dict=count_duplicates(capped_inter_num)
cap_dict={key: -value for key,value in cap_dict.items()}
#print(cap_dict)
#for i,j in cap_dict.items():
#    print("cap",i,j)
pair_inter_num=[]
for jjj in range(len(pair_xyz)):
    num=get_interaction_pair(pair_xyz[jjj])
    for kkk in num:
        pair_inter_num.append(kkk)
        
pair_dict=count_duplicates(pair_inter_num)
result = merge_dicts(pair_dict,cap_dict)
#print(result)

#same_dict=same_dicts(pair_dict,cap_dict)
#print(result)
#for i,j in pair_dict.items():
#    print("pair",i,j)

filter_dict=filter_dict_by_value(result,2)
double_ene=[]
for h in filter_dict:
    m=charge_txt[h[0]-1]
    n=charge_txt[h[1]-1]
    double_ene.append(compute_ele_two(m,n))

#print("pair energy:",len(pair_ene),sum(pair_ene))
#print("cap energy:",len(cap_ene),sum(cap_ene))
#print("double counting:",len(double_ene),sum(double_ene)/627.503)
print("1B_energy",sum(pair_ene)-sum(cap_ene)-(sum(double_ene)/627.503))

