import pyframe
import os
import glob
import numpy as np

# import the package 

class MFCCprotein:
    def __init__(self, pdbfile):
        self.pdbfile = pdbfile

    def mfcc(self):
        #print(self.pdbfile)
        systems = pyframe.MolecularSystem(self.pdbfile)
        os.system("mkdir -p tmpfile/caps tmpfile/capped")
        for fragment in systems.fragments.values():
            fragment.create_mfcc_fragments(bond_threshold=1.15)
            capsname = list(fragment.concaps.keys())[-1]
            fragment.concaps.write_xyz("tmpfile/caps/"+capsname)
            cappedname = fragment.capped_fragment.identifier
            fragment.capped_fragment.write_xyz("tmpfile/capped/"+cappedname)
        T=glob.glob("tmpfile/caps/*xyz")
        for i in T:
            jandk = i.split("/")[2].split(".xyz")[0].split("-")
            j = int(jandk[0].split("_")[0])
            k = int(jandk[1].split("_")[0])
            if j > k:
                #print("mv "+i+" caps/"+jandk[1]+"-"+jandk[0]+".xyz")
                os.system("mv "+i+" tmpfile/caps/"+jandk[1]+"-"+jandk[0]+".xyz")

class checkpro:
    def __init__(self, capslist, cappedlist):
        self.capslist = capslist
        self.cappedlist = cappedlist
    """
       capslist is the caps geometry
       cappedlist is the capped geometry
    """
    
    def sortcaps(self):
        capslist = self.capslist
        oldline = []
        for i in capslist:
            j = int(i.split("/")[2].split("_")[0])
            oldline.append([i, j])
        newline = sorted(oldline, key=lambda x :x[1])
        sortline = [h[0] for h in newline]
        return sortline
    """
    give caps list with number follow, follow the pyframe rules

    """
    def sortcapped(self):
        cappedlist = self.cappedlist
        oldline = []
        for i in cappedlist:
            j = int(i.split("/")[2].split("_")[0])
            oldline.append([i, j])
        newline = sorted(oldline, key=lambda x :x[1])
        sortline = [h[0] for h in newline]
        return sortline
    """
    give capped list with number follow
    """

    def checkcapslinktopro(self):
        sortcapslist = checkpro.sortcaps(self)

        def checkxyz(xyz):
            f = open(xyz).readlines()
            num = f[0].strip()
            return num
            
        def replaceCH3toH(xyz):
            index = [6, 7, 10]
            G = [n for i, n in enumerate(xyz) if i not in index]
            Noord = np.array([float(i) for i in G[2].split()[1:4]])
            Coord = np.array([float(i) for i in G[4].split()[1:4]])
            Hoord = (Coord-Noord)*(1.010/(np.linalg.norm(Coord-Noord)))+Noord
            Hoord = Hoord.tolist()
            Hoord = [str(round(i, 6)) for i in Hoord]
            Hoord = "    ".join(Hoord)
            # G[4]=G[4].replace("C","H")
            G[4] = "H    "+Hoord+"\n"
            return G
        """
        PRO residule need to consider its hydrogen in pyframe
        
        """

        def splitcaps(xyz, num):
            f = open(xyz).readlines()
            caps = f[:num]
            if num == 17:
                caps[0] = "12\n"
            else:
                caps[0] = "12\n"
            return caps
        
        Tot_caps = []
        for index in range(len(sortcapslist)):
            namelist = sortcapslist[index].split("/")[2].split("-")[0].split("_")
            name = namelist[0]+"_"+namelist[1]+"_CAP_"+namelist[2]
            if "PRO" in sortcapslist[index].split("_")[-3]:
                number1 = int(checkxyz(sortcapslist[index-1]))
                number2 = int(checkxyz(sortcapslist[index]))
                if number1 == 27 and number2 == 27:
                    pro_before = splitcaps(sortcapslist[index], 17)
                    fixxyz = replaceCH3toH(pro_before)
                    Tot_caps.append([name, fixxyz])
            else:
                Tot_caps.append([name, splitcaps(sortcapslist[index], 14)])
    
        os.system("mkdir -p tmpfile/capsnew")
        for i in Tot_caps:
            #print(i[0])
            extxyz = open("tmpfile/capsnew/"+i[0]+".xyz", "w")
            extxyz.write("".join(i[1]))
            extxyz.close()
    """
    need to fix pro and add into new 
    """

    def fixpro(self):
        sortcappedlist = checkpro.sortcapped(self)
        def fixpdblinktopro(xyz):
            def replaceCH3toHcapped(xyz):
                index = [3, 4, 6]
                G = [n for i, n in enumerate(xyz) if i not in index]
                Noord = np.array([float(i) for i in G[0].split()[1:4]])
                Coord = np.array([float(i) for i in G[1].split()[1:4]])
                Hoord = (Coord-Noord)*(1.010/(np.linalg.norm(Coord-Noord))) + Noord
                Hoord = Hoord.tolist()
                Hoord = [str(round(i,6)) for i in Hoord]
                Hoord = "    ".join(Hoord)
                # G[4]=G[4].replace("C","H")
                G[1] = "H     "+Hoord+"\n"
                return G
            xyzfile = open(xyz).readlines()
            unfixcoord = xyzfile[:-9]
            fixcoord = xyzfile[-9:]
            #print(fixcoord)
            #print(xyz)
            new = replaceCH3toHcapped(fixcoord)
            data = unfixcoord+new
            data[0] = str(len(data[2:]))+"\n"
            return data
        L2 = sortcappedlist
        os.system("mkdir -p tmpfile/cappednew")
        for index in range(len(L2)):
            if "PRO" in L2[index]:
                os.system("cp -r "+L2[index]+" tmpfile/cappednew")
                name = L2[index-1].split("/")[2]
                data = fixpdblinktopro(L2[index-1])
                H = open("tmpfile/cappednew/" + name, "w")
                H.write("".join(data))
                H.close()
            else:
                os.system("cp -r "+L2[index]+" tmpfile/cappednew")



