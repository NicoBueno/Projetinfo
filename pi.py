# -*- coding: utf-8 -*-
# Projet informatique 
# GROUPE nÂ°8 : Mendez De Maria, Bernon, Bueno
# Sujet : ETUDE DU MOUVEMENT D'UNE SONDE DANS UN CHAMP GRAVITATIONNEL

import matplotlib.pyplot as plt
import scipy.integrate as spi
import numpy as np
from math import cos 

#-------------------------------- Variables --------------------------------

Mp=5.97*(10**24)   # masse planete
Rt=6.6356*(10**6)  # rayon terre
h=800000
G=6.67*10**(-11)
y0=Rt+h
v0=7451.9
V=[0.8*v0,0.9*v0,v0,1.1*v0,1.2*v0]
alpha=-G*Mp 


#--------------------------- Fonctions de récurrence ------------------------- 

def fy(z,t):
    return z 

def fz(y,b):
    return (y*(b**2)+(alpha/(y**2)))
    
def fa(b):
     return b
     
def fb(z,b,y):
    return ((2*z*b)/y)
    
def f(var,t):
    [y,z,a,b]=var
    sol=[z,y*(b**2)+(alpha/(y**2)),b,-(2*z*b)/y]
    return sol

#--------------------- Fonctions donnant les solutions ------------------------- 
   
def recurrence(tf,n,i):    # ti et tf correspondent respectivement aux temps initial et final, ils sont exprimés en secondes
    z0=0
    a0=0
    b0=V[i]/y0
    p=float(tf)/float(n)
    T,Y,Z,A,B=[0.],[y0],[z0],[a0],[b0]    #Y correspond au rayon, Z à la vitesse
    t,y,z,a,b=0.,y0,z0,a0,b0              #A correspond à l'angle, B à la vitesse angulaire
    while (t+p)<=tf:
        t+=p
        y,z,a,b=y+p*fy(z,t),z+p*fz(y,b),a+p*fa(b),b-p*fb(z,b,y)
        T.append(t)
        Y.append(y)
        Z.append(z)
        A.append(a)
        B.append(b)
    return T,Y,A     #on veut uniquement l'évolution de y et de a en fonction du temps

def solscipy(tf,n,i):  #on trace la solution à l'equation differentielle grace a scipy
    b0=V[i]/y0
    CI=[y0,0,0,b0]
    ts=np.linspace(0.,tf,n)
    y=spi.odeint(f,CI,ts)
    Y=y[:,0]
    A=y[:,2]
    return ts,Y,A   


def solexacte(tf,n,i):  #l'argument permet de choisir la vitesse initiale
    Ye=[y0]
    b0=V[i]/y0
    C=((y0**2)*b0)
    d=(C**2)/(G*Mp)
    e=(d/y0)-1
    T,Y,A=recurrence(tf,n,i)
    y=y0
    for k in range(len(A)-1):
        y=d/(1+e*cos(A[k]))
        Ye.append(y)
    return T,Ye

#------------- Fonctions pour comparer graphiquement les ----------------------
#-----------------rayons et angles en fonction des CI ------------------------- 

def traceeuler(tf,n,i):   #trace rayon et angle en fonction d'une CI donnée
    t,y,a=recurrence(tf,n,i)
    plt.title("Rayon en fonction du temps")
    plt.plot(t,y,'r')
    plt.show()
    plt.title("Angle en fonction du temps")
    plt.plot(t,a,'b')
    plt.show()
    
def tracescipy(tf,n,i):
    t,y,a=solscipy(tf,n,i)
    plt.title("Rayon en fonction du temps")
    plt.plot(t,y,'r')
    plt.show()
    plt.title("Angle en fonction du temps")
    plt.plot(t,a,'b')
    plt.show()
    
def traceexacte(tf,n,i):
    t,y,a=solexacte(tf,n,i)
    plt.title("Rayon en fonction du temps")
    plt.plot(t,y,'r')
    plt.show()
    plt.title("Angle en fonction du temps")
    plt.plot(t,a,'b')
    plt.show()
    
#-------- Fonctions pour comparer graphiquement l'évolution des rayon et ----------------------
#------------- angle en fonction des diférentes CI avec une méthode------------------------- 
    
def compare1(tf,n):         #trace les rayons et angles obtenus avec Euler
    T,Y0,A0=recurrence(0,tf,n,0)
    T,Y1,A1=recurrence(0,tf,n,1)
    T,Y2,A2=recurrence(0,tf,n,2)
    T,Y3,A3=recurrence(0,tf,n,3)
    T,Y4,A4=recurrence(0,tf,n,4)
    R=[Rt for i in range(len(T))]
    plt.plot(T,Y0,'r',label="V(t=0)=5491 m/s")
    plt.plot(T,Y1,'g',label="V(t=0)=5491 m/s")
    plt.plot(T,Y2,'b',label="V(t=0)=5491 m/s")
    plt.plot(T,Y3,'y',label="V(t=0)=5491 m/s")
    plt.plot(T,Y4,'black',label="V(t=0)=5491 m/s")
    plt.plot(T,R,'purple',label="Surface de la Terre")
    plt.ylim(0.4*10**(7),1.2*10**(7))
    plt.xlabel("Temps(s)")
    plt.ylabel("Rayons(m)")
    plt.legend(loc=2)
    plt.show()
    plt.plot(T,A0,'r')
    plt.plot(T,A1,'g')
    plt.plot(T,A2,'b')
    plt.plot(T,A3,'y')
    plt.plot(T,A4,'k')
    plt.show()
    
def compare2(tf,n):         #trace les rayons et angles obtenus avec Odeint
    T,Y0,A0=solscipy(tf,n,0)
    T,Y1,A1=solscipy(tf,n,1)
    T,Y2,A2=solscipy(tf,n,2)
    T,Y3,A3=solscipy(tf,n,3)
    T,Y4,A4=solscipy(tf,n,4)
    R=[Rt for i in range(len(T))]
    plt.plot(T,Y0,'r',label="V(t=0)=5491 m/s")
    plt.plot(T,Y1,'g',label="V(t=0)=5491 m/s")
    plt.plot(T,Y2,'b',label="V(t=0)=5491 m/s")
    plt.plot(T,Y3,'y',label="V(t=0)=5491 m/s")
    plt.plot(T,Y4,'black',label="V(t=0)=5491 m/s")
    plt.plot(T,R,'purple',label="Surface de la Terre")
    plt.ylim(0.45*10**(7),1.2*10**(7))
    plt.xlabel("Temps(s)")
    plt.ylabel("Rayons(m)")
    plt.legend(loc=2)
    plt.show()
    plt.plot(T,A0,'r')
    plt.plot(T,A1,'g')
    plt.plot(T,A2,'b')
    plt.plot(T,A3,'y')
    plt.plot(T,A4,'k')
    plt.show()
    
def compare3(tf,n,i):   #trace les rayons et angles obtenus avec la solution exacte
    T,Y0=solexacte(tf,n,0)
    T,Y1=solexacte(tf,n,1)
    T,Y2=solexacte(tf,n,2)
    T,Y3=solexacte(tf,n,3)
    T,Y4=solexacte(tf,n,4)
    R=[Rt for i in range(len(T))]
    plt.plot(T,Y0,'r',label="V(t=0)=5491 m/s")
    plt.plot(T,Y1,'g',label="V(t=0)=5491 m/s")
    plt.plot(T,Y2,'b',label="V(t=0)=5491 m/s")
    plt.plot(T,Y3,'y',label="V(t=0)=5491 m/s")
    plt.plot(T,Y4,'black',label="V(t=0)=5491 m/s")
    plt.plot(T,R,'purple',label="Surface de la Terre")
    plt.ylim(0.4*10**(7),1.2*10**(7))
    plt.xlabel("Temps(s)")
    plt.ylabel("Rayons(m)")
    plt.legend(loc=2)
    plt.show()
    
#-------- Fonctions pour comparer graphiquement les solutions ------------------
#------------- obtenues grâce aux différentes méthodes ------------------------- 

def compare4(tf,n,i):         #compare les solutions obtenues avec les 3 méthodes 
    T,Y0,A0=recurrence(tf,n,i)
    T1,Y1,A1=solscipy(tf,n,i)
    T2,Y2=solexacte(tf,n,i)
    R=[Rt for i in range(len(T))]
    plt.plot(T,Y0,'r',label="V(0)=5491 m/s")
    plt.plot(T1,Y1,'g',label="V(0)=5491 m/s")
    plt.plot(T2,Y2,'b',label="V(0)=5491 m/s")
    plt.plot(T,R,'purple',label="Surface de la Terre")
    plt.xlabel("Temps(s)")
    plt.ylabel("Rayons(m)")
    plt.show()
    plt.plot(T,A0,'r',label="V(0)=5491 m/s")
    plt.plot(T1,A1,'g',label="V(0)=5491 m/s")
    plt.xlabel("Temps(s)")
    plt.ylabel("Angle(°)")
    plt.show()
    
def compare5(i,n,tf):  #compare les différentes solutions en fonction de la vitesse initiale et du nombre de points d'acquisitions
    T,Y,A=recurrence(tf,n,i)
    ts=np.linspace(0.,tf,n)
    Ys=[k[0] for k in spi.odeint(f,tf,ts)]
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