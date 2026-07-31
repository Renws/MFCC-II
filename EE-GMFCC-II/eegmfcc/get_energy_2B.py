import os
import glob 
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--dir", type=str, help="pdbfile")
args = parser.parse_args()


#print(args.dir+"/dimertmpfile/charge_dimer/")

dimer=glob.glob(args.dir+"/charge_dimer/*log")
#print(dimer)
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
    ene_tot,self_ene=ene[-1],ene[0]
    return float(ene_tot)-float(self_ene)

def split_file_by_empty_lines(input_filename):
    index=[]
    for line in range(len(input_filename)):
        line1 = input_filename[line].strip()
        if not line1:
            index.append(line)
    return index

def get_k_inter(mono1,mono2):
    mono1_xyz=open(mono1).readlines()[2:]
    mono2_xyz=open(mono2).readlines()[2:]
    mono1_index=split_file_by_empty_lines(mono1_xyz)
    mono2_index=split_file_by_empty_lines(mono2_xyz)
    #print(mono1_index)
    mono1_xyz_charge=mono1_xyz[mono1_index[-1]+1:]
    mono2_xyz_charge=mono2_xyz[mono2_index[-1]+1:]
    #print(len(mono1_xyz_charge))
    energy=[]
    for ii in mono1_xyz_charge:
        xyz1=np.array([float(coord) for coord in ii.strip().split()[:3]])
        charge1=float(ii.strip().split()[3])
        #print(xyz1,charge1)
        for jj in mono2_xyz_charge:
            xyz2=np.array([float(coord) for coord in jj.strip().split()[:3]])
            charge2=float(jj.strip().split()[3])
            R12=np.linalg.norm(xyz1 - xyz2)/0.5291772083
            ele_unit=(charge1*charge2)/R12
            energy.append(ele_unit)
    #print(sum(energy))
    return sum(energy)

energy_sum=[]
K_sum=[]
for h in dimer:
    dimer_ene=get_energy(h)
    h=h.split("/")[-1]
    #print(h)
    mono1=args.dir+"/charge_monomer/"+h.split(".")[0].split("-")[0]+".xyz.log"
    mono2=args.dir+"/charge_monomer/"+h.split(".")[0].split("-")[1]+".xyz.log"    
    #print(mono1,mono2)
    mono1_xyz=args.dir+"/charge_monomer/"+h.split(".")[0].split("-")[0]+".xyz"
    mono2_xyz=args.dir+"/charge_monomer/"+h.split(".")[0].split("-")[1]+".xyz"
    K_energy = get_k_inter(mono1_xyz,mono2_xyz)
    K_sum.append(K_energy)
    mono1_ene,mono2_ene =  get_energy(mono1),get_energy(mono2)
    energy_sum.append(dimer_ene-mono1_ene-mono2_ene)
    #print(h,dimer_ene,mono1_ene,mono2_ene)
    #print(h,dimer_ene-mono1_ene-mono2_ene,mono1_ene,mono2_ene)
#print("2B energy")
print("2B_energy",len(energy_sum),sum(energy_sum)+sum(K_sum))
