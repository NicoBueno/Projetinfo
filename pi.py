# Projet informatique 
# GROUPE n°8 : Mendez De Maria, Bernon, Bueno
# Sujet : ETUDE DU MOUVEMENT D'UNE SONDE DANS UN CHAMP GRAVITATIONNEL

import matplotlib.pyplot as plt
import scipy
from scipy import integrate

####### Variables

m=20               # masse de la sonde
Mp=5.97*(10**24)   # masse planete
Rt=6.6356*(10**6)  # rayon terre
v0=7451.9
V=[0.8*v0,0.9*v0,v0,1.1*v0,1.2*v0]
h=800000
#alpha=-G*Mp # donner valeur a G

####### fonctions 

def recurrence(ti,tf,n,y0,z0,a0,b0,i):
    p=(tf-ti)/n    # pas
    T,Y,Z,A,B=[t0],[y0],[z0],[a0],[b0]
    t,y,z,a,b=t0,y0,z0,a0,b0 
    while (t+p)<=tf:
        t+=p
        y,z,a,b=y+p*z,z+p*(y*V[i]-alpha)/(y**2),a+p*V[i]/y,b+p*(-2*z*a)/y
        T.append(t)
        Y.append(y)
        Z.append(z)
        A.append(a)
        B.append(b)
    return T,Y,Z,A,B
    
def trace(T,Y,Z):
    t,x,y,a,b=recurrence(ti,tf,n,y0,z0,a0,b0,i)  # à compléter
    plt.title("Angle et rayon en fonction du temps")
    plt.plot(x,t,'r',y,t,'b')     #on trace les solutions obtenues grâce à python
    
def vitesse(V):
    v1,v2,v3,v4,v5=V[0],V[1],V[2],V[3],V[4]
    return v1,v2,v3,v4,v5