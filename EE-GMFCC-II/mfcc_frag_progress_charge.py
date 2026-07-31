import argparse
from MFCC import mfccprotein
from MFCC.pairamino import monotopair
import shutil
import os
import glob
import numpy as np
import pyframe
import time
from tqdm import tqdm
start_time = time.time()
##  parameter setting 
parser = argparse.ArgumentParser()
parser.add_argument("--pdb", type=str, help="pdbfile")
parser.add_argument("--parm", type=str, help="parm")
parser.add_argument("--crd", type=str, help="crd")
parser.add_argument("--dimerdir",type=str,help="dimer dirname",default="dimertmpfile")
parser.add_argument("--dist",type=float,help="dimer distance",default=5.0)
parser.add_argument("--no-progress", action="store_true", help="disable progress bar")
args = parser.parse_args()

# ------------------ Overall progress bar ------------------ #
# The total is updated dynamically because some workloads, e.g. the number of
# dimers, are only known after intermediate files have been generated.
pbar = tqdm(total=0, desc="Overall", unit="step", dynamic_ncols=True, disable=args.no_progress)

def add_progress_total(n):
    n = int(n)
    if n <= 0:
        return
    pbar.total = (pbar.total or 0) + n
    pbar.refresh()

def progress_update(n=1, msg=None):
    if msg:
        pbar.set_postfix_str(str(msg))
    pbar.update(n)

## MFCC deal with
pdb = mfccprotein.MFCCprotein(args.pdb)
add_progress_total(1)
pbar.set_description("MFCC")
pdb.mfcc()
progress_update(1, "mfcc done")
caps_list = glob.glob("tmpfile/caps/*xyz")
capped_list = glob.glob("tmpfile/capped/*xyz")
caps2_checkpro = mfccprotein.checkpro(caps_list, capped_list)
add_progress_total(2)
pbar.set_description("Check caps")
caps2_checkpro.checkcapslinktopro()
progress_update(1, "check caps link")
caps2_checkpro.fixpro()
progress_update(1, "fix pro")


cappednew_list = glob.glob("tmpfile/cappednew/*xyz")
cappednew_list = sorted(cappednew_list,key=lambda x:int(x.split("/")[2].split("_")[0]))

cappednew_list1 = sorted(cappednew_list,key=lambda x:int(x.split("/")[2].split("_")[0]))
cappednew_list1=cappednew_list1[1:-1]

add_progress_total(1)
pbar.set_description("Mono to pair")
monotopair(cappednew_list1)
progress_update(1, "pair amino generated")
cappedpair_list = glob.glob("tmpfile/pairamino/*xyz")
cappedpair_list = sorted(cappedpair_list,key=lambda x:int(x.split("/")[2].split("_")[0]))

#pdb_amber_charge=open(args.pdb).readlines()
#pdb_amber_charge=pdb_amber_charge[:-2]
pdb_amber_parm=open(args.parm).readlines()


def round_z(number, decimals):
    return f"{number:.{decimals}f}"

pdb_crd=open(args.crd).readlines()[2:]
pdb_crd=[round_z(float(j),3) for i in pdb_crd for j in i.strip().split()]
pdb_crd_new=[]

for pp in range(0,len(pdb_crd),3):
    pdb_crd_new.append(pdb_crd[pp:pp+3])

parm_charge_index=[]
for i in range(len(pdb_amber_parm)):
    if "%FLAG CHARGE" in pdb_amber_parm[i]:
        parm_charge_index.append(i)
    elif "%FLAG ATOMIC_NUMBER" in pdb_amber_parm[i]:
        parm_charge_index.append(i)
    elif "%FLAG MASS" in pdb_amber_parm[i]:
        parm_charge_index.append(i)

parm_charge = pdb_amber_parm[parm_charge_index[0]+2:parm_charge_index[1]]
parm_atomic_number = pdb_amber_parm[parm_charge_index[1]+2:parm_charge_index[2]]

charge_list=[]
for j in parm_charge:
    H=j.strip().split()
    for h in H:
        charge_list.append(float(h))

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
CHARGES = dict(((i,x) for i,x in enumerate(ELEMENTS)))

elements=[CHARGES[int(j)] for i in parm_atomic_number for j in i.strip().split()]

pdb_xyz_dict = {tuple(coord): ["1", str(i + 1), elem, coord, chg / np.sqrt(332.0637)]
                for i, (elem, coord, chg) in enumerate(zip(elements, pdb_crd_new, charge_list))}

# ------------------ write all-atom coordinate/charge/index file ------------------ #
def write_charge_txt(output_path="tmpfile/charge.txt"):
    """Write all atoms as: x y z charge atom_index.

    Coordinates are read from args.crd and rounded to 3 decimals above.
    Charges are read from args.parm and converted with the same factor used
    elsewhere in this script: charge / sqrt(332.0637).
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    if len(pdb_crd_new) != len(charge_list):
        print(
            f"Warning: coordinate count ({len(pdb_crd_new)}) and charge count "
            f"({len(charge_list)}) are different. charge.txt will contain "
            f"{min(len(pdb_crd_new), len(charge_list))} matched atoms."
        )

    with open(output_path, "w") as f_out:
        for atom_index, (coord, chg) in enumerate(zip(pdb_crd_new, charge_list), start=1):
            scaled_charge = chg / np.sqrt(332.0637)
            f_out.write(f"{' '.join(coord)} {scaled_charge} {atom_index}\n")

add_progress_total(1)
pbar.set_description("Write charge.txt")
write_charge_txt("tmpfile/charge.txt")
progress_update(1, "write tmpfile/charge.txt")

os.makedirs("tmpfile/charge_cappednew/", exist_ok=True)
os.makedirs("tmpfile/charge_pairamino/", exist_ok=True)

def process_file(file_list, output_dir, progress_desc=None):
    progress_desc = progress_desc or os.path.basename(os.path.normpath(output_dir))
    add_progress_total(len(file_list))
    pbar.set_description("Process files")
    for file_path in file_list:
        lines = [line.strip() for line in open(file_path).readlines()]
        atoms_lines = lines[2:]
        charge_repeat, charge_part = [], []

        for line in atoms_lines:
            parts = line.strip().split()
            coords = tuple([round_z(float(parts[1]), 3), round_z(float(parts[2]), 3), round_z(float(parts[3]), 3)])
            match = pdb_xyz_dict.get(coords)
            if match:
                entry = f"   {' '.join(match[3])} {match[4]} {match[1]}"
                charge_repeat.append(entry)
            else:
                for val in pdb_xyz_dict.values():
                    entry = f"   {' '.join(val[3])} {val[4]} {val[1]}"
                    charge_part.append(entry)

        # sort and filter
        charge_repeat_set = set(charge_repeat)
        unique_result = list(set(charge_part) - charge_repeat_set)
        # sort
        result_sorted = []
        for coord in pdb_xyz_dict.keys():
            coord_str = ' '.join(coord)
            for entry in unique_result:
                if coord_str in entry:
                    result_sorted.append(entry)
        # write file
        output_lines = lines + [""] + result_sorted + [""] + charge_repeat
        with open(os.path.join(output_dir, os.path.basename(file_path)), "w") as f_out:
            f_out.write("\n".join(output_lines))
        progress_update(1, progress_desc)

def check_seq_index(lst):
    index_new=[]
    for i in range(1, len(lst)):
        lst_f=lst[i].strip().split()
        lst_b=lst[i-1].strip().split()
        if int(lst_f[4]) - int(lst_b[4]) != 1:
            index_new.append([i-1,i])
            #return (i-1, i)
    return index_new

def process_file_mono(file_list, output_dir, progress_desc=None):
    progress_desc = progress_desc or os.path.basename(os.path.normpath(output_dir))
    add_progress_total(len(file_list))
    pbar.set_description("Process mono files")
    for file_path in file_list:
        lines = [line.strip() for line in open(file_path).readlines()]
        atoms_lines = lines[2:]
        charge_repeat, charge_part = [], []

        for line in atoms_lines:
            parts = line.strip().split()
            coords = tuple([round_z(float(parts[1]), 3), round_z(float(parts[2]), 3), round_z(float(parts[3]), 3)])
            match = pdb_xyz_dict.get(coords)
            if match:
                entry = f"   {' '.join(match[3])} {match[4]} {match[1]}"
                charge_repeat.append(entry)
            else:
                for val in pdb_xyz_dict.values():
                    entry = f"   {' '.join(val[3])} {val[4]} {val[1]}"
                    charge_part.append(entry)

        # sort and filter
        charge_repeat_set = set(charge_repeat)
        unique_result = list(set(charge_part) - charge_repeat_set)
        # sort
        result_sorted = []
        for coord in pdb_xyz_dict.keys():
            coord_str = ' '.join(coord)
            for entry in unique_result:
                if coord_str in entry:
                    result_sorted.append(entry)
        
        index=check_seq_index(result_sorted)
        print(index,file_path)
        if index!=None and len(index)==1:
            #print(index)
            index_del=[index[0][0]-1,index[0][0],index[0][1],index[0][1]+1]
            print(index)
            result_sorted1 = [item for idx, item in enumerate(result_sorted) if idx not in index_del]
            output_lines = lines + [""] + result_sorted1 + [""] + charge_repeat
            with open(os.path.join(output_dir, os.path.basename(file_path)), "w") as f_out:
                f_out.write("\n".join(output_lines))
        elif len(index)==0:
            if "ACE" in file_path:
                index_del=[0,1]
                print(index_del)
                result_sorted1 = [item for idx, item in enumerate(result_sorted) if idx not in index_del]
                output_lines = lines + [""] + result_sorted1 + [""] + charge_repeat
                with open(os.path.join(output_dir, os.path.basename(file_path)), "w") as f_out:
                    f_out.write("\n".join(output_lines))
            elif "NME" in file_path:
                index_del=[len(result_sorted)-2,len(result_sorted)-1]
                print(index_del)
                result_sorted1 = [item for idx, item in enumerate(result_sorted) if idx not in index_del]
                output_lines = lines + [""] + result_sorted1 + [""] + charge_repeat
                with open(os.path.join(output_dir, os.path.basename(file_path)), "w") as f_out:
                    f_out.write("\n".join(output_lines))
        progress_update(1, progress_desc)

def process_file_dimer(file_list, output_dir, progress_desc=None):
    progress_desc = progress_desc or os.path.basename(os.path.normpath(output_dir))
    add_progress_total(len(file_list))
    pbar.set_description("Process dimer files")
    for file_path in file_list:
        lines = [line.strip() for line in open(file_path).readlines()]
        atoms_lines = lines[2:]
        charge_repeat, charge_part = [], []

        for line in atoms_lines:
            parts = line.strip().split()
            coords = tuple([round_z(float(parts[1]), 3), round_z(float(parts[2]), 3), round_z(float(parts[3]), 3)])
            match = pdb_xyz_dict.get(coords)
            if match:
                entry = f"   {' '.join(match[3])} {match[4]} {match[1]}"
                charge_repeat.append(entry)
            else:
                for val in pdb_xyz_dict.values():
                    entry = f"   {' '.join(val[3])} {val[4]} {val[1]}"
                    charge_part.append(entry)

        # sort and filter
        charge_repeat_set = set(charge_repeat)
        unique_result = list(set(charge_part) - charge_repeat_set)
        # sort
        result_sorted = []
        for coord in pdb_xyz_dict.keys():
            coord_str = ' '.join(coord)
            for entry in unique_result:
                if coord_str in entry:
                    result_sorted.append(entry)
        # write file
        index=check_seq_index(result_sorted)
        if index!=None and len(index)==2:
            #print(index)
            index_del=[index[0][0]-1,index[0][0],index[0][1],index[0][1]+1,index[1][0]-1,index[1][0],index[1][1],index[0][1]+1]
            result_sorted1 = [item for idx, item in enumerate(result_sorted) if idx not in index_del]
            output_lines = lines + [""] + result_sorted1 + [""] + charge_repeat
            with open(os.path.join(output_dir, os.path.basename(file_path)), "w") as f_out:
                f_out.write("\n".join(output_lines))
        elif len(index)==1:
            if "ACE" in file_path:
                index_del=[index[0][0]-1,index[0][0],index[0][1],index[0][1]+1,0,1]
                #print(index_del)
                result_sorted1 = [item for idx, item in enumerate(result_sorted) if idx not in index_del]
                output_lines = lines + [""] + result_sorted1 + [""] + charge_repeat
                with open(os.path.join(output_dir, os.path.basename(file_path)), "w") as f_out:
                    f_out.write("\n".join(output_lines))
            elif "NME" in file_path:
                index_del=[index[0][0]-1,index[0][0],index[0][1],index[0][1]+1,len(result_sorted)-1,len(result_sorted)-2]
                #print(index_del)
                result_sorted1 = [item for idx, item in enumerate(result_sorted) if idx not in index_del]
                output_lines = lines + [""] + result_sorted1 + [""] + charge_repeat
                with open(os.path.join(output_dir, os.path.basename(file_path)), "w") as f_out:
                    f_out.write("\n".join(output_lines))
        progress_update(1, progress_desc)

process_file(cappedpair_list, "tmpfile/charge_pairamino/", "charge pair amino")
process_file(cappednew_list, "tmpfile/charge_cappednew", "charge cappednew")

shutil.rmtree("tmpfile/capped")
shutil.rmtree("tmpfile/cappednew")
shutil.rmtree("tmpfile/caps")
shutil.rmtree("tmpfile/capsnew")
shutil.rmtree("tmpfile/pairamino")

## dimer 
systems = pyframe.MolecularSystem(args.pdb)
dirname = args.dimerdir

os.makedirs(args.dimerdir + "/monomer", exist_ok=True)
os.makedirs(args.dimerdir + "/dimer", exist_ok=True)
os.makedirs(args.dimerdir + "/charge_monomer", exist_ok=True)
os.makedirs(args.dimerdir + "/charge_dimer", exist_ok=True)

#idnameL=[]
# generate monomer
fragments = list(systems.fragments.values())
add_progress_total(len(fragments))
pbar.set_description("Generate monomers")
for fragment in fragments:
    #idname=fragment.identifier
    fragment.create_mfcc_fragments(order=0)
    fragment.capped_fragment.write_xyz()
    progress_update(1, "generate monomer")

xyz_files = glob.glob("*.xyz")
add_progress_total(len(xyz_files))
pbar.set_description("Move monomers")
for file in xyz_files:
    shutil.move(file, os.path.join(args.dimerdir, "monomer", os.path.basename(file)))
    progress_update(1, "move monomer xyz")

#os.system("mkdir -p "+str(dirname)+"/monomer")
#os.system("mv *.xyz "+str(dirname)+"/monomer")
#os.system("mkdir -p "+str(dirname)+"/dimer")

def load_xyz_coordinates(file_list):
    coord_dict = {}
    add_progress_total(len(file_list))
    pbar.set_description("Load coordinates")
    for f in file_list:
        coord_dict[f] = np.loadtxt(f, skiprows=2, usecols=(1, 2, 3))
        progress_update(1, "load xyz coordinates")
    return coord_dict


monomer_files = sorted(glob.glob(os.path.join(args.dimerdir, "monomer", "*xyz")),
                       key=lambda x: int(os.path.basename(x).split("_")[0]))

monomer_coords = load_xyz_coordinates(monomer_files)

def compute_min_distance(coords1, coords2):
    distances = np.linalg.norm(coords1[:, np.newaxis] - coords2, axis=2)
    return np.min(distances)

dimerlist = []
num_dimer_candidates = max(0, (len(monomer_files) - 1) * (len(monomer_files) - 2) // 2)
add_progress_total(num_dimer_candidates)
pbar.set_description("Search dimers")
for i, fi in enumerate(monomer_files[:-2]):
    for fj in monomer_files[i + 2:]:
        min_dis = compute_min_distance(monomer_coords[fi], monomer_coords[fj])
        if min_dis < args.dist:
            dimerlist.append((fi, fj, min_dis))
        progress_update(1, "search dimer pairs")
# get dimer xyz
add_progress_total(len(dimerlist))
pbar.set_description("Write dimer xyz")
for fi, fj, min_dis in dimerlist:
    name1, name2 = os.path.splitext(os.path.basename(fi))[0], os.path.splitext(os.path.basename(fj))[0]
    if abs(int(name2.split("_")[0]) - int(name1.split("_")[0])) > 1:
        dimer_name = f"{name1}-{name2}"
        coords1 = open(fi).readlines()[2:]
        coords2 = open(fj).readlines()[2:]
        combined = [str(len(coords1) + len(coords2)), f"{dimer_name} {min_dis:.3f}"] + [line.strip() for line in coords1 + coords2]
        with open(os.path.join(args.dimerdir, "dimer", f"{dimer_name}.xyz"), "w") as fout:
            fout.write("\n".join(combined))
    progress_update(1, "write dimer xyz")

# get dimer charge and monomer charge ------------------ #
dimer_files = sorted(glob.glob(os.path.join(args.dimerdir, "dimer", "*xyz")),
                     key=lambda x: int(os.path.basename(x).split("_")[0]))
dimer_files=dimer_files
process_file(dimer_files, os.path.join(args.dimerdir, "charge_dimer"), "charge dimer")
process_file(monomer_files, os.path.join(args.dimerdir, "charge_monomer"), "charge monomer")

#print("Caps: ",len(cappednew_list),"Frag:",len(cappedpair_list), "Dimer:", len(dimer_files), "Monomer:", len(monomer_files),"time:",time.time()-start_time)

shutil.rmtree(os.path.join(args.dimerdir, "monomer"))
shutil.rmtree(os.path.join(args.dimerdir, "dimer"))

#result_dir = os.path.splitext(os.path.basename(args.pdb))[0]
result_dir = os.path.normpath(args.pdb).split(os.sep)[0]
os.makedirs(result_dir, exist_ok=True)
add_progress_total(2)
pbar.set_description("Finalize")
shutil.move(args.dimerdir, result_dir)
progress_update(1, "move dimerdir")
shutil.move("tmpfile", result_dir)
progress_update(1, "move tmpfile")
pbar.close()
print("Caps:", len(cappednew_list), "Frag:", len(cappedpair_list), "Dimer:", len(dimer_files), "Monomer:", len(monomer_files), "time:", time.time() - start_time)
