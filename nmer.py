import sys
from itertools import combinations
#import psi4


cov_rad = {   'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71,
  'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 'Br': 1.14, 'I' : 1.33, 'He': 0.30,
  'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 'Na': 1.02,
  'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75,
  'Ti': 0.86, 'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
  'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
  'Se': 1.17, 'Br': 1.14, 'Kr': 1.03}

def write_geom(name, geom):
    with open(name, 'w') as fil:
        lg = len(geom)
        fil.write(str(lg)+'\n\n')
        for i in range(lg):
            atom = geom[i]
            lin = atom[0]
            for c in range(1,4): lin+= " "+str(atom[c])
            fil.write(lin+'\n')

def dist(a,b):
    a_coords = a[1:] 
    b_coords = b[1:]
    xdist = (a_coords[0]-b_coords[0])**2
    ydist = (a_coords[1]-b_coords[1])**2
    zdist = (a_coords[2]-b_coords[2])**2
    dist = (xdist+ydist+zdist)**0.5
    return dist


def bound(a,b):
    if cov_rad[a[0]]+cov_rad[b[0]] >=  dist(a,b): return True
    else: return False

def intersect(frag,frags):
    intersects = set() 
    for f in frags:
        if frag & f == set([]): continue
        intersects = intersects | frag | f
        frags.remove(f) 
    return intersects
    
def bfs(geom):
    minifrags = []
    #make sets of bonded atoms
    for i in range(len(geom)):
        for j in range(i):
            if bound(geom[i],geom[j]): minifrags.append(set([str(geom[i]),str(geom[j])]))
    frags = []
    while minifrags != []:
        one = intersect(minifrags[0],minifrags)
        last = []
        while one != set([]):
            one = intersect(one,minifrags)
            if one != set([]): last = one
        if last != []: frags.append(last)
    
    return frags

def frag2str(frag):
    #newfrag = ''
    #frag = list(frag)
    #for i in frag:
    #    el = i.split("'")[1]
    #    coords = i.split(']')[0]
    #    coords = ''.join(coords.split(',')[1:])
    #    newfrag += el + coords + '\n'
    #return newfrag
    newfrag = []
    frag = list(frag)
    for i in frag:
        el = i.split("'")[1]
        coords = i.split(']')[0]
        coords = ''.join(coords.split(',')[1:])
        coords = [float(x) for x in coords.split()]
        newfrag.append([el,coords[0],coords[1],coords[2]])
    return newfrag


def rmsd(A,B):
    ans = 0.
    for i in range(len(A)): ans += ((A[i][0] - B[i][0])**2 + (A[i][1] - B[i][1])**2 + (A[i][2] - B[i][2])**2)
    ans /= float(len(A))
    return ans**0.5

el2mass = {'H': 1.00782503207, 'HE': 4.00260325415, 'LI': 7.016004548, 'BE': 9.012182201, 'B': 11.009305406, 'C': 12, 'N': 14.00307400478, 'O': 15.99491461956, 'F': 18.998403224, 'NE': 19.99244017542, 'NA': 22.98976928087, 'MG': 23.985041699, 'AL': 26.981538627, 'SI': 27.97692653246, 'P': 30.973761629, 'S': 31.972070999, 'CL': 34.968852682, 'AR': 39.96238312251, 'K': 38.963706679, 'CA': 39.962590983, 'SC': 44.955911909, 'TI': 47.947946281, 'V': 50.943959507, 'CR': 51.940507472, 'MN': 54.938045141, 'FE': 55.934937475, 'CO': 58.933195048, 'NI': 57.935342907, 'CU': 62.929597474, 'ZN': 63.929142222, 'GA': 68.925573587, 'GE': 73.921177767, 'AS': 74.921596478, 'SE': 79.916521271, 'BR': 78.918337087, 'KR': 85.910610729, 'RB': 84.911789737, 'SR': 87.905612124, 'Y': 88.905848295, 'ZR': 89.904704416, 'NB': 92.906378058, 'MO': 97.905408169, 'TC': 98.906254747, 'RU': 101.904349312, 'RH': 102.905504292, 'PD': 105.903485715, 'AG': 106.90509682, 'CD': 113.90335854, 'IN': 114.903878484, 'SN': 119.902194676, 'SB': 120.903815686, 'TE': 129.906224399, 'I': 126.904472681, 'XE': 131.904153457, 'CS': 132.905451932, 'BA': 137.905247237, 'LA': 138.906353267, 'CE': 139.905438706, 'PR': 140.907652769, 'ND': 141.907723297, 'PM': 144.912749023, 'SM': 151.919732425, 'EU': 152.921230339, 'GD': 157.924103912, 'TB': 158.925346757, 'DY': 163.929174751, 'HO': 164.93032207, 'ER': 165.930293061, 'TM': 168.93421325, 'YB': 173.938862089, 'LU': 174.940771819, 'HF': 179.946549953, 'TA': 180.947995763, 'W': 183.950931188, 'RE': 186.955753109, 'OS': 191.96148069, 'IR': 192.96292643, 'PT': 194.964791134, 'AU': 196.966568662, 'HG': 201.970643011, 'TL': 204.974427541, 'PB': 207.976652071, 'BI': 208.980398734, 'PO': 208.982430435, 'AT': 210.987496271, 'RN': 222.017577738, 'FR': 222.01755173, 'RA': 228.031070292, 'AC': 227.027752127, 'TH': 232.038055325, 'PA': 231.03588399, 'U': 238.050788247, 'NP': 237.048173444, 'PU': 242.058742611, 'AM': 243.06138108, 'CM': 247.07035354, 'BK': 247.07030708, 'CF': 251.079586788, 'ES': 252.082978512, 'FM': 257.095104724, 'MD': 258.098431319, 'NO': 255.093241131, 'LR': 260.105504, 'RF': 263.112547, 'DB': 255.107398, 'SG': 259.1145, 'BH': 262.122892, 'HS': 263.128558, 'MT': 265.136151, 'DS': 281.162061, 'RG': 272.153615, 'UUB': 283.171792, 'UUT': 283.176451, 'UUQ': 285.183698, 'UUP': 287.191186, 'UUH': 292.199786, 'UUS': 291.206564, 'UUO': 293.21467}

def com(A,els):
    totmass = 0.
    center = [0.0, 0.0, 0.0]
    for i,el in enumerate(els): 
        for c in range(3): center[c] += el2mass[el] * A[i][c]  
        totmass += el2mass[el]
    for c in range(3): center[c] = center[c]/totmass
    return center

def kabsch(sys0,sys1):
        nsys0 = [] 
        els0 = []
        for i in sys0: els0.append(i[0].upper())
        for i in sys0: nsys0.append([x for x in i[1:]])
        nsys1 = [] 
        els1 = []
        for i in sys1: els1.append(i[0].upper())
        for i in sys1: nsys1.append([x for x in i[1:]])
        #Translate sys1 and sys0 to origin
        center = com(nsys0,els0)
        center1 = com(nsys1,els1)
        for i in range(len(nsys1)):
            for j in range(3): 
                nsys0[i][j] -= center[j]
                nsys1[i][j] -= center1[j]
        #compute cross-covariance matrix
        H = np.dot(nsys0,np.array(nsys1).T)
        #H = np.dot(np.array(nsys0).T,np.array(nsys1))
        #SVD to get optimal rotation
        V,s,W = np.linalg.svd(H)
        d = np.linalg.det(np.dot(W,V))
        if d < 0: V[:,-1] *= -1.
        R = np.dot(V,W) 
        #Rotate
        #nsys1 = np.dot(nsys1,R)
        nsys1 = np.dot(R,nsys1)
        diff = rmsd(nsys0,nsys1)
        return diff

import numpy as np

#for i in range(1):
#    mol0 = sys.argv[1]
#    mol1 = sys.argv[2]
#    sys0 = []
#    sys1 = []
#    with open(mol0,'r') as fil:
#        for lin in fil:
#            if len(lin.split()) != 4: continue
#            sys0.append(lin.split('\n')[0])
#    with open(mol1,'r') as fil:
#        for lin in fil:
#            if len(lin.split()) != 4: continue
#            sys1.append(lin.split('\n')[0])
## make all dimers and print nre
##temporarily taking in two xyz's
##maxn = 5
##for n in range(2,maxn+1):
#        #for i in range(len(nsys0)): print els0[i],' '.join(map(str,nsys0[i]))
#        #for i in range(len(nsys1)): print els0[i],' '.join(map(str,nsys1[i]))
#    rms = kabsch(sys0,sys1)       
#    print(rms)
def distance(atom1, atom2):
    x = atom1[1]-atom2[1]
    y = atom1[2]-atom2[2]
    z = atom1[3]-atom2[3]
    d = (x*x + y*y + z*z)**0.5
    d = 1.88973*d
    return d

temp_symbol = ["X", "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG",
"AL", "SI", "P", "S", "CL", "AR", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO",
"NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR", "NB",
"MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE", "CS",
"BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM",
"YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI",
"PO", "AT", "RN", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK",
"CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG",
"UUB", "UUT", "UUQ", "UUP", "UUH", "UUS", "UUO"]

z_num = {}
for i,el in enumerate(temp_symbol): z_num[el] = i

def nuclear_repulsion_energy(mol):
    e = 0
    for a1 in range(len(mol)):
        for a2 in range(a1):
            z1 = z_num[mol[a1][0]]
            z2 = z_num[mol[a2][0]]
            dist = distance(mol[a1],mol[a2])
            e += z1 * z2 / dist
    return e 

def choose(n,k):
    topprod = 1
    for i in range(n-k+1,n+1): topprod *= i
    bottomprod = 1
    for i in range(1,k+1): bottomprod *= i
    return topprod // bottomprod

def get_nmers(geo_fil_name, N):
    """ Returns a dictionary of lists that contain Molecule objects for
        all unique n-mers from n=1-N. Key values correspond to the value of 
        n.
    >>> H2O_20.get_nmers(1)
    {1: [<qcdb.molecule.Molecule object at 0x7f3a66929f10>, <qcdb.molecule.Molecule object at 0x7f3a66929ed0>]}
    """
    #Read Geometry
    with open(geo_fil_name,'r') as fil: 
        n = int(next(fil))
        next(fil)
        geom = []
        for i in range(n):
            atom = next(fil).split()
            for i in range(1,4): atom[i] = float(atom[i])
            geom.append(atom)
    #Generate Monomers
    frags = bfs(geom)
    geoms = [frag2str(i) for i in frags]

    from itertools import combinations
    
    
    nfrags = len(geoms)
    if N > nfrags: 
        raise Exception("BFS only found %d fragments. Cannot find %d-mer with %d fragments." % (nfrags,N,nfrags)) 

    inds = [x for x in range(nfrags)]
    
    #dictionary with lists for each class of n-mer
    unique = {}
    for i in range(N): unique[i+1] = []

    #for each type of n-mer
    for i in range(N):
        #add each unique n-mer
        nres = {}
        for combo in combinations(inds,i+1):
            mol = []
            for nm in list(combo): mol += geoms[nm]
            nre = nuclear_repulsion_energy(mol)
            #tweak to whatever tightness we'd likee
            nre = round(nre,2)                    
            try: nres[nre].append(mol)
            except: nres[nre] = [mol]
        for nre in nres.keys():
            lnre = len(nres[nre])
            #last element not checked for duplication 
            unique[i+1].append(nres[nre][-1])
            #nothing to compare too if only one 
            if lnre == 1: continue
            for m1 in range(lnre-1):
                duplicate = False
                for m2 in range(m1+1,lnre):
                    mol1 = nres[nre][m1] 
                    mol2 = nres[nre][m2]
                    rmsd = kabsch(mol1,mol2) 
                    #make keyword for this option
                    if rmsd < 0.01: 
                        duplicate = True
                        break
                if not duplicate: unique[i+1].append(mol1)
    return unique, len(geom) 

unis, tot_mons = get_nmers(sys.argv[1], 3) 

uk = sorted(unis.keys())


import os
try: os.mkdir("Unique_Nmers")
except: pass

print("Total number of monomers: %d" % (tot_mons))

os.chdir("Unique_Nmers")
for k in uk: 
    print("Total number of unique %d-mers out of a possible %10d: %5d" % (k,choose(tot_mons,k),len(unis[k])))
    try: os.mkdir(str(k))
    except: pass
    os.chdir(str(k))
    for num, geom in enumerate(unis[k]): write_geom(str(num+1)+'.xyz', geom)
    os.chdir("..")

