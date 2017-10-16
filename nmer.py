import sys
from itertools import combinations
#import psi4

geo_fil = sys.argv[1]
with open(geo_fil,'r') as fil: 
    n = int(next(fil))
    next(fil)
    geom = []
    for i in range(n):
        atom = next(fil).split()
        for i in range(1,4): atom[i] = float(atom[i])
        geom.append(atom)

cov_rad = {   'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71,
  'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 'Br': 1.14, 'I' : 1.33, 'He': 0.30,
  'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 'Na': 1.02,
  'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75,
  'Ti': 0.86, 'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
  'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
  'Se': 1.17, 'Br': 1.14, 'Kr': 1.03}

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
    newfrag = ''
    frag = list(frag)
    for i in frag:
        el = i.split("'")[1]
        coords = i.split(']')[0]
        coords = ''.join(coords.split(',')[1:])
        newfrag += el + coords + '\n'
    return newfrag

frags = bfs(geom)
geoms = [frag2str(i) for i in frags]

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

# make all dimers and print nre
for i in combinations(geoms,2):
    sys0 = i[0].split('\n')[:-1]
    sys1 = i[1].split('\n')[:-1]
    nsys0 = [] 
    els0 = []
    for i in sys0: els0.append(i.split()[0].upper())
    for i in sys0: nsys0.append([float(x) for x in i.split()[1:]])
    nsys1 = [] 
    els1 = []
    for i in sys1: els1.append(i.split()[0].upper())
    for i in sys1: nsys1.append([float(x) for x in i.split()[1:]])
    for i in range(len(els0)): 
        if els0[i] != els1[i]: 
            print("ALL ELEMENTS NEED ALIGN, TECH NOT READY YET")
            sys.exit()
    #Translate sys1 to com of sys0
    center = com(nsys0,els0)
   
