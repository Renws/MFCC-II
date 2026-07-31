import os
import glob 
import numpy as np
import argparse
import time

_t0 = time.time()
parser = argparse.ArgumentParser()
parser.add_argument("--dir", type=str, help="pdbfile")
args = parser.parse_args()

capped=glob.glob(args.dir+"tmpfile/charge_cappednew/*xyz.log")
capped = sorted(capped,key=lambda x:int(x.split("/")[-1].split("_")[0]))
capped = capped[2:-2]
capped_xyz=glob.glob(args.dir+"tmpfile/charge_cappednew/*xyz")
capped_xyz_raw = sorted(capped_xyz,key=lambda x:int(x.split("/")[-1].split("_")[0]))
capped_xyz = capped_xyz_raw[2:-2]

pair=glob.glob(args.dir+"tmpfile/charge_pairamino/*xyz.log")
pair= sorted(pair,key=lambda x:int(x.split("/")[-1].split("_")[0]))
pair_xyz = glob.glob(args.dir+"tmpfile/charge_pairamino/*xyz")
pair_xyz = sorted(pair_xyz,key=lambda x:int(x.split("/")[-1].split("_")[0]))

dimer=glob.glob(args.dir+"dimertmpfile/charge_dimer/*.xyz.log")
dimer.sort()
monomer=glob.glob(args.dir+"dimertmpfile/charge_monomer/*.xyz.log")
monomer= sorted(monomer,key=lambda x:int(x.split("/")[-1].split("_")[0]))
#print(monomer)

def get_energy(log):
    log=open(log).readlines()
    ene=[]
    for i in log:
        if "SCF Done:" in i:
            energy=i.strip().split()[4]
            ene.append(energy)
    ene_tot = ene[0]
    return float(ene_tot)




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
    return sum(ene)/627.503

def get_inter_frag_bc(xyz_raw):
    xyz=open(xyz_raw).readlines()[2:]
    index1=split_file_by_empty_lines(xyz)
    #print(index1)
    frag_xyz_charge=xyz[index1[1]+1:]
    frag_xyz_charge_new=[]
    for i in frag_xyz_charge:
        line=[float(j) for j in i.strip().split()[:4]]
        frag_xyz_charge_new.append(line)
    bc_xyz_charge=xyz[index1[0]+1:index1[1]]
    #print(bc_xyz_charge) 
    bc_xyz_charge_new=[]
    for k in bc_xyz_charge:
        line1=[float(l) for l in k.strip().split()[:4]]
        bc_xyz_charge_new.append(line1)
    #print(bc_xyz_charge_new) 
    return np.array(frag_xyz_charge_new), np.array(bc_xyz_charge_new)

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



pair_ene=[]
pair_K_ene=[]

for i in range(len(pair)):
    qm_ene=get_energy(pair[i])
    k_ene=compute_ele(get_inter_frag_bc(pair_xyz[i])[0],get_inter_frag_bc(pair_xyz[i])[1])
    pair_ene.append(qm_ene)
    pair_K_ene.append(k_ene) 
print(sum(pair_ene),sum(pair_K_ene),len(pair_ene))

capped_ene=[]
capped_K_ene=[]
for i in range(len(capped)):
    qm_ene=get_energy(capped[i])
    k_ene=compute_ele(get_inter_frag_bc(capped_xyz[i])[0],get_inter_frag_bc(capped_xyz[i])[1])
    capped_ene.append(qm_ene)
    capped_K_ene.append(k_ene) 
print(sum(capped_ene),sum(capped_K_ene),len(capped_ene))

dimer_ene=[]
dimer_K_ene=[]

M1_ene=[]
M1_K_ene=[]
M2_ene=[]
M2_K_ene=[]


tot_int=[]
for hh in range(len(monomer)-2):
    for  gg in range(hh+2,len(monomer)):
        mono1=monomer[hh]
        mono2=monomer[gg]
        mono1_xyz=mono1.split(".log")[0]
        mono2_xyz=mono2.split(".log")[0]
        K_energy = get_k_inter(mono1_xyz,mono2_xyz)
        tot_int.append(K_energy)
        #print(monomer[hh],monomer[gg])



inter_K=[]

for h in dimer:
    dimer_qm_ene=get_energy(h)
    dimer_xyz=h.split(".log")[0]
    #dimer_k_ene=compute_ele(get_inter_frag_bc(dimer_xyz)[0],get_inter_frag_bc(dimer_xyz)[1])
    dimer_ene.append(dimer_qm_ene)
    #dimer_K_ene.append(dimer_k_ene)

    h=h.split("/")[-1]    
    mono1=args.dir+"dimertmpfile/charge_monomer/"+h.split(".")[0].split("-")[0]+".xyz.log"
    mono2=args.dir+"dimertmpfile/charge_monomer/"+h.split(".")[0].split("-")[1]+".xyz.log"
    
    mono1_xyz=args.dir+"dimertmpfile/charge_monomer/"+h.split(".")[0].split("-")[0]+".xyz"
    mono2_xyz=args.dir+"dimertmpfile/charge_monomer/"+h.split(".")[0].split("-")[1]+".xyz"
    
    mono1_ene=get_energy(mono1)
    mono1_k_ene=compute_ele(get_inter_frag_bc(mono1_xyz)[0],get_inter_frag_bc(mono1_xyz)[1])
    M1_ene.append(mono1_ene)
    M1_K_ene.append(mono1_k_ene)

    mono2_ene=get_energy(mono2)
    mono2_k_ene=compute_ele(get_inter_frag_bc(mono2_xyz)[0],get_inter_frag_bc(mono2_xyz)[1])
    M2_ene.append(mono2_ene)
    M2_K_ene.append(mono2_k_ene)
    
    inter_k=compute_ele(get_inter_frag_bc(mono1_xyz)[0],get_inter_frag_bc(mono2_xyz)[0])
    inter_K.append(inter_k)

#print(sum(dimer_ene),sum(dimer_K_ene),len(dimer_ene))
#print(sum(M1_ene),sum(M1_K_ene),len(M1_ene))
#print(sum(M2_ene),sum(M2_K_ene),len(M2_ene))

pair_ene_tot=sum(pair_ene)
capped_ene_tot=sum(capped_ene)
dimer_ene_tot=sum(dimer_ene)
M1_ene_tot=sum(M1_ene)
M2_ene_tot=sum(M2_ene)

E_ass=pair_ene_tot-capped_ene_tot+dimer_ene_tot-M1_ene_tot-M2_ene_tot+sum(tot_int)-sum(inter_K)
run_time = time.time() - _t0

print("tot_ene",E_ass*27.21138)
print(f"run_time：{run_time:.3f}s")

