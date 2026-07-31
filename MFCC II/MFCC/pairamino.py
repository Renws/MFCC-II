import os
import glob

def monotopair(capped_list):
    for i in range(0,len(capped_list)-1):
        if "PRO" in capped_list[i+1]:
            T0=open(capped_list[i]).readlines()[2:]
            T1=open(capped_list[i+1]).readlines()[2:]
            frag1_index=[-1,-2,-3,-4,-5]
            del_index_infrag1=[len(T0)-5, len(T0)-4,len(T0)-3,len(T0)-2,len(T0)-1]
            T0=[n for i,n in enumerate(T0) if i not in del_index_infrag1]
            start_index_frag2=[0]
            end_index_frag2=[-12, -11, -10, -9, -8, -7]
            end_index_frag2=[i + len(T1) for i in end_index_frag2]
            del_index_infrag2=start_index_frag2+end_index_frag2
            #print(del_index_infrag2)
            T1=[n for i,n in enumerate(T1) if i not in del_index_infrag2]
            name0=capped_list[i].split("/")[2].split(".")[0]
            name1=capped_list[i+1].split("/")[2].split(".")[0]
            name=name0+"-"+name1
            Tot_two=T0+T1
            Tot_two=[i.strip() for i in Tot_two]
            xyztitle=[str(len(Tot_two)),name]
            Tot_amino=xyztitle+Tot_two
            os.system("mkdir -p tmpfile/pairamino/")
            f=open("tmpfile/pairamino/"+name+".xyz","w")
            f.write("\n".join(Tot_amino))

        else:
            T0=open(capped_list[i]).readlines()[2:]
            T1=open(capped_list[i+1]).readlines()[2:]
            frag1_index=[-1,-2,-3,-4]
            del_index_infrag1=[len(T0),len(T0)-4,len(T0)-3,len(T0)-2,len(T0)-1]
            T0=[n for i,n in enumerate(T0) if i not in del_index_infrag1]
            start_index_frag2=[0, 2]
            end_index_frag2=[-12, -11, -10, -9, -8, -7]
            end_index_frag2=[i + len(T1) for i in end_index_frag2]
            del_index_infrag2=start_index_frag2+end_index_frag2
            T1=[n for i,n in enumerate(T1) if i not in del_index_infrag2]
            name0=capped_list[i].split("/")[2].split(".")[0]
            name1=capped_list[i+1].split("/")[2].split(".")[0]
            name=name0+"-"+name1
            Tot_two=T0+T1
            Tot_two=[i.strip() for i in Tot_two]
            xyztitle=[str(len(Tot_two)),name]
            Tot_amino=xyztitle+Tot_two
            os.system("mkdir -p tmpfile/pairamino/")
            f=open("tmpfile/pairamino/"+name+".xyz","w")
            f.write("\n".join(Tot_amino))
