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
G=6.67*10**(-11)   # reste à voir pour avoir la valeur par scipy
r0=Rt+h
v0=7451.9
V=[0.8*v0,0.9*v0,v0,1.1*v0,1.2*v0]
al=-G*Mp
# ECRITURE DES VARIABLES
# r pour r
# R pour rpoint
# th pour theta
# TH pour thetapoint
# al pour alpha
# un L devant le nom de la variable quand il s agit d une liste
# un 0 a la fin du nom de la variable pour la valeur initiale
# t pour le temp, et T pour liste de temp
#------------------------------------------------------ Fonctions -------------------------------------------------- 
def fr(R):                # je pense que cette fonction est inutile, on devrait mettre fr directment dans recurrence
    return R 
def fR(r,TH):
    return (r*(TH**2)+(al/(r**2)))
def fth(TH):       # je pense que cette fonction est inutile, on devrait mettre TH directment dans recurrence
     return TH   
def fTH(R,TH,r):    # dans odeint, on ne peut entrer que deux argument, dont le temps en deuxieme ...
    return ((2*R*TH)/r)
    
#---------------------------------------------------------Euler-------------------------------------------------------- 
def recurrence(ti,tf,n,i):                     # ti => t init  tf => t final   secondes
    R0=0
    th0=0
    TH0=V[i]/r0
    p=float((tf-ti))/float(n)
    T,Lr,LR,Lth,LTH=[float(ti)],[r0],[R0],[th0],[TH0]   #Y:rayon, Z:vitesse, A:angle, B:vitesse angulaire
    t,r,R,th,TH=ti,r0,R0,th0,TH0                    
    for i in range (n):                                  # et non pas while (t+p)<=tf:
        t+=p
        r,R,th,TH=r+p*fr(R),R+p*fR(r,TH),th+p*fth(TH),TH-p*fTH(R,TH,r)   # je pense que l erreur est ici mais je vois pas où. Je pense que c est au niveau des arguments que ça foire
        T.append(t)
        Lr.append(r)
        LR.append(R)
        Lth.append(th)
        LTH.append(TH)
    return T,Lr,Lth                #on veut uniquement l'évolution de r et de theta en fonction du temps


#--------------------------------------------------------------------------------------------------------
# en fait odeint a besoin dune fonction, de la condition initale du terme dont tu veut faire la liste, et dune liste de temps
# qui correspond aux points où le calcul sera effectue .Mais si , comme dans la fonction fR ou dans la fonction fTH
# plusieurs variables changent en meme temps, je ne pense pas que ça fonctionne, il faut donc, je pense,
#pour utiliser odeint, creer les listes des focntions fth et fTH en même temps grace a une boucle for
# je met donc l ancienne def en commentaire et jessaie une nouvelle def
#def solscipy(f,tf,n,condinit):           #on trace la solution à l'equation differentielle grace a scipy
 #   ts=np.linspace(0.,tf,n)              # de 0 a tf en n points
    #y=spi.odeint(f,condinit,ts)          # fonction, condition init (et non pas tf !!!), liste de temps
    #Y=y[:,0]                           # sert à rien
    #plt.plot(ts,y,"b")
    #plt.show()
def solscipy(tf,n,i):
    R0=0
    th0=0          # probleme, on utilise pas cette cond init alors qu on devrait je pense
    TH0=V[i]/r0
    ts=np.linspace(0.,tf,n) 
    yr,yTH,yR,yth=[],[],[],[]
    for i in ts:                   # ça fonctionne pas encore tres bien mais on devrait garder l idee
        yr.append((spi.odeint(fr,R0,[i]))[-1])
        yTH.append((spi.odeint(fTH,TH0,[i]))[-1])
        yR.append((spi.odeint(fR,[r0,yTH[-1]],[i]))[-1])
        yth.append((spi.odeint(fTH,[R0,TH0,yr[-1]],[i]))[-1])
    plt.plot(ts,yR,"b")
    plt.show()
    plt.plot(ts,yth,"b")
    plt.show()
    plt.plot(ts,yr,"b")
    plt.show()
    plt.plot(ts,yTH,"b")
    plt.show()
    
def solexacte(tf,n,i):               #l'argument permet de choisir la vitesse initiale
    Ye=[r0]
    b0=V[i]/r0
    C=((r0**2)*b0)
    d=(C**2)/(G*Mp)
    e=(d/r0)-1
    T,Lr,Lth=recurrence(0,tf,n,i)
    r=r0
    for k in range(len(Lth)-1):
        y=d/(1+e*cos(Lth[k]))
        Ye.append(y)
    return T,Ye


    
def traceeuler(tf,n,i):      #on prend en argument le nombre de points d'acquisition et la vitesse que l'on veut 
    t,r,th=recurrence(0,tf,n,i)
    plt.title("Rayon en fonction du temps")
    plt.plot(t,r,'r')
    plt.show()
    plt.title("Angle en fonction du temps")
    plt.plot(t,th,'theta')
    plt.show()
    
    
    
def compare1(i,n,tf,f,condinit):       #compare les différentes solutions en fonction de la vitesse initiale et du nombre de points d'acquisitions
    T,Lr,Lth=recurrence(0,1,n,i)
    ts=np.linspace(0.,tf,n)
    Ys=[k[0] for k in spi.odeint(f,condinit,ts)]    # ce n est pas tf qu il faut rentrer en deuxieme argument, mais la condition initiale de la variable a etudier
    # te=T                sert a rien
    #te.append()          sert a rien
    Ye=solexacte(tf,n,i)      #  solexacte prend 3 arguments donc jai rajoute tf 
    plt.title("Rayons (Euler en rouge, Scipy en vert et Exacte en bleu) en fonction du temps")
    plt.xlabel("Temps (s)")
    plt.ylabel("Rayons (m)")
    plt.plot(T,Lr,'r',ts,Ys,'g',T,Ye,'b')
    plt.legend()
    plt.show()

def compare2(tf,n):         #trace les rayons et angles en fonction de v0 (avec Euler)
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
    
def compare3(tf,n,i):         #trace les rayons et angles en fonction de v0 (avec la solution exacte)
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
    
def compare4(tf,n,i):         #compare les solutions oobtenues avec les 3 méthodes 
    T,Y0,A0=recurrence(0,tf,n,i)
    T1,Y1=solscipy(fy,tf,n)
    T2,Y2=solexacte(tf,n,i)
    R=[Rt for i in range(len(T))]
    plt.plot(T,Y0,'r',label="V(t=0)=5491 m/s")
    plt.plot(T1,Y1,'g',label="V(t=0)=5491 m/s")
    plt.plot(T2,Y2,'b',label="V(t=0)=5491 m/s")
    plt.plot(T,R,'purple',label="Surface de la Terre")
    plt.xlabel("Temps(s)")
    plt.ylabel("Rayons(m)")
    #plt.legend(loc=2)
    plt.show()
