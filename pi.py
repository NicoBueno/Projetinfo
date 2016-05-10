# -*- coding: utf-8 -*-
# Projet informatique 
# GROUPE nÂ°8 : Mendez De Maria, Bernon, Bueno
# Sujet : ETUDE DU MOUVEMENT D'UNE SONDE DANS UN CHAMP GRAVITATIONNEL

from matplotlib.pyplot import *
import scipy.integrate as spi

##### Variables
m=20               # masse de la sonde
Mp=5.97*(10**24)   # masse planete
Rt=6.6356*(10**6)  # rayon terre
v0=7451.9
V=[0.8*v0,0.9*v0,v0,1.1*v0,1.2*v0]
h=800000
G=9.81
alpha=G*Mp 


# fonctions 

def f1(z,t):
    return z 

def f2(y,i):
    return (y*(V[i]**2)-(alpha/(y**2)))
    
def f3(b,t):
     return b
     
def f4(z,a,y):
    return ((-2*z*a)/y)
   
def recurrence(ti,tf,n,z0,a0,i):
    y0=Rt+h
    b0=(V[i])/(Rt+h)
    p=(tf-ti)/n                        # pas
    T,Y,Z,A,B=[ti],[y0],[z0],[a0],[b0]
    t,y,z,a,b=ti,y0,z0,a0,b0 
    while (t+p)<=tf:
        t+=p
        y,z,a,b=y+p*f1(z),z+p*f2(y,i),a+p*f3(b),b+p*f4(z,a,y)
        T.append(t)
        Y.append(y)
        Z.append(z)
        A.append(a)
        B.append(b)
    return T,Y,Z,A,B
    
def trace(n,i): #on prend en argument le nombre de points d'acquisition et la vitesse que l'on veut 
    t,y,z,a,b=recurrence(0,1,n,0,0,i)  # à compléter
    title("Rayon en fonction du temps")
    plot(t,y,'r')     #on trace les solutions obtenues grâce à python
    show()
    plot(t,a,'b')    
    show()
    
#def vitesse(V):
#    v1,v2,v3,v4,v5=V[0],V[1],V[2],V[3],V[4]
#    return v1,v2,v3,v4,v5
    
def resolution(f,tf):
    t=np.linspace(0.,tf,10)
    y=spi.odeint(f,tf,t)
    Y=y[:,0]
    plot(t,Y,'g')

def compare():
    
