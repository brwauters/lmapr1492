#!/usr/bin/env python
# coding: utf-8

# In[57]:


from pymatgen.ext.matproj import MPRester
import numpy

with MPRester("K9id1PwVOvn68mCE") as m:
    results = m.query("mp-561619", ['cif'])[0]


# In[58]:


a=open("local.cif",'w')
b=results.items()
for d,w in b:
    a.write(str(w))
a.close()


# In[59]:


from pymatgen.io.cif import CifParser
parser = CifParser("local.cif")
structure = parser.get_structures(True)[0]


# In[60]:


from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from IPython.display import Image
from pymatgen.symmetry.structure import SymmetrizedStructure



#La cellule symmetrised est une cellule conventionnelle
n=SpacegroupAnalyzer(structure).get_conventional_standard_structure(international_monoclinic=True)


#Retourne le dataset des symétries sous forme d'un dict
j=SpacegroupAnalyzer(n).get_symmetry_dataset()

#Dans le cas où les vecteurs de translations sortiraient de la cellule conventionnelle on fait -1 pour être conforme avec la représentation donnée sur le site du prof
#on transforme ainsi les vecteurs(1/3,2/3,2/3)et (2/3,1/3,1/3) en (1/3,-1/3,-1/3) et (-1/3,1/3,1/3)
for i in range(0,36):
    for k in range(0,3):
        if (j['translations'][i][k]>0.4):
            j['translations'][i][k]=j['translations'][i][k]-1

for i in range(0,36):
    for k in range(0,3):
        if (j['translations'][i][k]<-0.4):
            j['translations'][i][k]=j['translations'][i][k]+1
            
            
#créé un objet SymmOp qui est issu d'une matrice de rotation donnée dans "j" et d'un vecteur de translation
l1=SymmOp.from_rotation_and_translation(j['rotations'][7],j['translations'][30],tol=0.01)#32

#retourne les équations caractérisant l'application
u1=l1.as_xyz_string()

#Effectue la transformation sur l'atome a la position donnée en argument ici :[0.33,0.67,0.67] et retorune sa position après la transformation
p1=l1.operate([0.33,0.67,0.67])#9

#Affiche l'image en png de la transformation via le site linker par le professeur notons que ces images sont disponibles en grand format dans le dossier img du notebook
print (" le premier élément de symétrie est : ",u1,"\n","les coord de l'atome sur lequel elle est appliquée  : ", [0.33,0.67,0.67],"\n","les coord de l'atome après application de l'élément de symmétrie : ",p1)
Image(filename="img/SYM1-32-9.PNG",width=500,height=500)




# In[61]:


from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from IPython.display import Image
from pymatgen.symmetry.structure import SymmetrizedStructure



#La cellule symmetrised est une cellule conventionnelle
n=SpacegroupAnalyzer(structure).get_conventional_standard_structure(international_monoclinic=True)


#Retourne le dataset des symétries sous forme d'un dict
j=SpacegroupAnalyzer(n).get_symmetry_dataset()

#Dans le cas où les vecteurs de translations sortiraient de la cellule conventionnelle on fait -1 pour être conforme avec la représentation donnée sur le site du prof
for i in range(0,36):
    for k in range(0,3):
        if (j['translations'][i][k]>0.4):
            j['translations'][i][k]=j['translations'][i][k]-1

for i in range(0,36):
    for k in range(0,3):
        if (j['translations'][i][k]<-0.4):
            j['translations'][i][k]=j['translations'][i][k]+1
            
            
#créé un objet SymmOp qui est issu d'une matrice de rotation donnée dans "j" et d'un vecteur de translation
l2=SymmOp.from_rotation_and_translation(j['rotations'][7],j['translations'][15],tol=0.01)#23

#retourne les équations caractérisant l'application
u2=l2.as_xyz_string()

#Effectue la transformation sur l'atome a la position donnée en argument ici :[0.00,1.00,0.50] et retourne sa position après la transformation
#Les coordonnées de l'atome sont recopiées depuis l'objet 'structure'
p2=l2.operate([1.00,1.00,1.00])#15

#Affiche l'image en png de la transformation via le site linker par le professeur notons que ces images sont disponibles en grand format dans le dossier img du notebook
print (" le second élément de symétrie est : ",u2,"\n","les coord de l'atome sur lequel elle est appliquée  : ", [1.00,1.00,1.00],"\n","les coord de l'atome apès application de l'élément de symmétrie : ",p1)
Image(filename="img/SYM2.PNG",width=500,height=500)




# In[62]:


from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from IPython.display import Image



#La cellule symmetrised est une cellule conventionnelle
n=SpacegroupAnalyzer(structure).get_conventional_standard_structure(international_monoclinic=True)

#Retourne le dataset des symétries sous forme d'un dict
j=(SpacegroupAnalyzer(n).get_symmetry_dataset())

#créé un objet SymmOp qui est issu d'une matrice de rotation donnée dans "j" et d'un vecteur de translation
l3=SymmOp.from_rotation_and_translation(j['rotations'][5],j['translations'][0],tol=0.01)#32

#retourne les équations caractérisant l'application
u3=l3.as_xyz_string()

#Effectue la transformation sur l'atome a la position donnée en argument ici :[1.00,0.00,0.00] et retorune sa position après la transformation
p3=l3.operate([1.00,0.00,0.00])#9

#Affiche l'image en png de la transformation via le site linker par le professeur notons que ces images sont disponibles en grand format dans le dossier img du notebook
print (" le troisième élément de symétrie est : ",u3,"\n","les coord de l'atome sur lequel elle est appliquée  : ", [1.00,0.00,0.00],"\n","les coord de l'atome apès application de l'élément de symmétrie : ",p3)
Image(filename="img/SYM3-.PNG",width=500,height=500)


