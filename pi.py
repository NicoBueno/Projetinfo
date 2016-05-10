# -*- coding: utf-8 -*-
# Projet informatique 
# GROUPE nÂ°8 : Mendez De Maria, Bernon, Bueno
# Sujet : ETUDE DU MOUVEMENT D'UNE SONDE DANS UN CHAMP GRAVITATIONNEL

import matplotlib.pyplot as plt
import scipy.integrate as spi
import numpy as np
from numpy import cos 

#------------------------------------------------------- Variables -------------------------------------------------

Mp=5.97*(10**24)   # masse planete
Rt=6.6356*(10**6)  # rayon terre
h=800000
G=9.81
y0=Rt+h
v0=7451.9
V=[0.8*v0,0.9*v0,v0,1.1*v0,1.2*v0]
alpha=-G*Mp 


#------------------------------------------------------ Fonctions -------------------------------------------------- 

def f1(z,t):
    return z 

def f2(y,b):
    return (y*(b**2)+(alpha/(y**2)))
    
def f3(b):
     return b
     
def f4(z,b,y):
    return ((2*z*b)/y)
   
def recurrence(ti,tf,n,z0,a0,i):
    b0=V[i]/y0
    p=float((tf-ti))/float(n)
    T,Y,Z,A,B=[ti],[y0],[z0],[a0],[b0]    #Y correspond au rayon, Z à la vitesse
    t,y,z,a,b=ti,y0,z0,a0,b0              #A correspond à l'angle, B à la vitesse angulaire
    while (t+p)<=tf:
        t+=p
        y,z,a,b=y+p*f1(z,t),z+p*f2(y,b),a+p*f3(b),b-p*f4(z,b,y)
        T.append(t)
        Y.append(y)
        Z.append(z)
        A.append(a)
        B.append(b)
    return T,Y,A                #on veut uniquement l'évolution de y et de a en fonction du temps

def solscipy(f,tf,n):           #on trace la solution à l'equation differentielle grace a scipy
    ts=np.linspace(0.,tf,n)
    y=spi.odeint(f,tf,ts)
    Y=y[:,0]
    plt.plot(ts,Y,'g')

def solexacte(n,i):               #l'argument permet de choisir la vitesse initiale
    Ye=[y0]
    T,Y,A=recurrence(0,1,n,0,0,i)
    y=y0
    for k in A:
        y+=(((y0*V[i])**2)/(G*Mp))/(1+((((((y0**2)*V[i])**2)/(G*Mp))/y0)-1)*cos(k))
        Ye.append(y)
    return Ye



    
def traceeuler(tf,n,i):      #on prend en argument le nombre de points d'acquisition et la vitesse que l'on veut 
    t,y,a=recurrence(0,tf,n,0,0,i)
    plt.title("Rayon en fonction du temps")
    plt.plot(t,y,'r')
    plt.show()
    plt.title("Angle en fonction du temps")
    plt.plot(t,a,'b')
    plt.show()


def compare1(i,n,tf,f):       #compare les différentes solutions en fonction de la vitesse initiale et du nombre de points d'acquisitions
    T,Y,A=recurrence(0,1,n,0,0,i)
    ts=np.linspace(0.,tf,n)
    Ys=spi.odeint(f,tf,ts)
    Ye=solexacte(n,i)
    T=np.array(T)
    Y=np.array(Y)
    plt.title("Rayons (Euler en rouge, Scipy en vert et Exacte en bleu) en fonction du temps")
    plt.xlabel("Temps (s)")
    plt.ylabel("Rayons (m)")
    plt.plot(T,Y,'r',T,Ys,'g',T,Ye,'b')
    plt.legend()
    plt.show()

#def compare2
    
