# -*- coding: utf-8 -*-
# Projet informatique 
# GROUPE nÂ°8 : Mendez De Maria, Bernon, Bueno
# Sujet : ETUDE DU MOUVEMENT D'UNE SONDE DANS UN CHAMP GRAVITATIONNEL

import matplotlib.pyplot as plt
import scipy.integrate as spi
import numpy as np
from math import cos 

#------------------------------------------------------- Variables -------------------------------------------------

Mp=5.97*(10**24)   # masse planete
Rt=6.6356*(10**6)  # rayon terre
h=800000
G=6.67*10**(-11)
y0=Rt+h
v0=7451.9
V=[0.8*v0,0.9*v0,v0,1.1*v0,1.2*v0]
alpha=-G*Mp 


#------------------------------------------------------ Fonctions -------------------------------------------------- 

def fy(z,t):
    return z 

def fz(y,b):
    return (y*(b**2)+(alpha/(y**2)))
    
def fa(b):
     return b
     
def fb(z,b,y):
    return ((2*z*b)/y)

   
def recurrence(ti,tf,n,i):
    z0=0
    a0=0
    b0=V[i]/y0
    p=float((tf-ti))/float(n)
    T,Y,Z,A,B=[float(ti)],[y0],[z0],[a0],[b0]    #Y correspond au rayon, Z à la vitesse
    t,y,z,a,b=ti,y0,z0,a0,b0              #A correspond à l'angle, B à la vitesse angulaire
    while (t+p)<=tf:
        t+=p
        y,z,a,b=y+p*fy(z,t),z+p*fz(y,b),a+p*fa(b),b-p*fb(z,b,y)
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

def solexacte(tf,n,i):               #l'argument permet de choisir la vitesse initiale
    Ye=[]
    d=(y0*(V[i]**2))/(G*Mp)
    e=(d/y0)-1
    T,Y,A=recurrence(0,tf,n,i)
    y=y0
    k=0
    for k in A:
        y+=d/(1+e*cos(k))
        Ye.append(y)
    plt.plot(T,Ye,'b')
    plt.show()
    #return Ye


    
def traceeuler(tf,n,i):      #on prend en argument le nombre de points d'acquisition et la vitesse que l'on veut 
    t,y,a=recurrence(0,tf,n,i)
    plt.title("Rayon en fonction du temps")
    plt.plot(t,y,'r')
    plt.show()
    plt.title("Angle en fonction du temps")
    plt.plot(t,a,'b')
    plt.show()
    
    
    
def compare1(i,n,tf,f):       #compare les différentes solutions en fonction de la vitesse initiale et du nombre de points d'acquisitions
    T,Y,A=recurrence(0,1,n,i)
    print(T)
    print(Y)
    ts=np.linspace(0.,tf,n)
    Ys=[i[0] for i in spi.odeint(f,tf,ts)]
    te=T
    te.append()
    Ye=solexacte(n,i)
    print(Ye)
    print(Ys)
    plt.title("Rayons (Euler en rouge, Scipy en vert et Exacte en bleu) en fonction du temps")
    plt.xlabel("Temps (s)")
    plt.ylabel("Rayons (m)")
    plt.plot(T,Y,'r',ts,Ys,'g',T,Ye,'b')
    plt.legend()
    plt.show()

def compare2(tf,n):         #trace les rayons et angles en fonction de v0
    T,Y0,A0=recurrence(0,tf,n,0)
    T,Y1,A1=recurrence(0,tf,n,1)
    T,Y2,A2=recurrence(0,tf,n,2)
    T,Y3,A3=recurrence(0,tf,n,3)
    T,Y4,A4=recurrence(0,tf,n,4)
    R=[Rt for i in range(len(T))]
    plt.plot(T,Y0,'r',T,Y1,'g',T,Y2,'b',T,Y3,'y',T,Y4,'black',T,R,'purple')
    plt.show()
    plt.plot(T,A0,'r',T,A1,'g',T,A2,'b',T,A3,'y',T,A4,'black')
    plt.show()