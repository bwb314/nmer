import sys
from itertools import combinations
import psi4

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

# make all dimers and print nre
for i in combinations(geoms,2):
    nmer = ''.join(i)
    pmol = psi4.geometry(nmer)
    nre = pmol.nuclear_repulsion_energy()
    print(nre)

