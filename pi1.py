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
P0,P1,P2,P3,P4=0.001,0.1,1,10,1000

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
    
#--------------------------- Fonction qui trace sur un même graphique --------------------
#--------------------------les différentes listes qu'elle prend en argument ------------------------- 
    
def trace1(L1,L2,L3,L4,L5,R,n):
    if L1[2]==L2[2]:
        plt.plot(L1[0],L1[1],'r',label="P="+str(0.01*n))
        plt.plot(L2[0],L2[1],'g',label="P="+str(0.1*n))
        plt.plot(L3[0],L3[1],'b',label="P="+str(n))
        plt.plot(L4[0],L4[1],'y',label="P="+str(10*n))
        plt.plot(L5[0],L5[1],'black',label="P="+str(100*n))
        plt.title("V(0)="+str(V[L1[2]])+"m/s")
    else:
        plt.plot(L1[0],L1[1],'r',label="V(0)="+str(V[L1[2]])+"m/s")
        plt.plot(L2[0],L2[1],'g',label="V(0)="+str(V[L2[2]])+"m/s")
        plt.plot(L3[0],L3[1],'b',label="V(0)="+str(V[L3[2]])+"m/s")
        plt.plot(L4[0],L4[1],'y',label="V(0)="+str(V[L4[2]])+"m/s")
        plt.plot(L5[0],L5[1],'black',label="V(0)="+str(V[L5[2]])+"m/s")
    plt.plot(R[0],R[1],'purple',label="Surface de la Terre")
    plt.ylim(0.4*10**(7),1.2*10**(7))
    plt.xlabel("Temps(s)")
    plt.ylabel("Rayons(m)")
    plt.legend(loc=2)
    plt.show()

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
    return T,Ye,A

#------------- Fonctions pour comparer graphiquement les ----------------------
#-----------------rayons et angles en fonction des CI ------------------------- 

def traceeuler(tf,n,i):   #trace rayon et angle en fonction d'une CI donnée
    t,y,a=recurrence(tf,n,i)
    plt.title("Rayon en fonction du temps")
    plt.plot(t,y,'r')
    plt.show()
    plt.title("Angle en fonction du temps")
    plt.plot(t,a,'b')
    plt.savefig("GrapheEuler.png")
    plt.show()
    
def tracescipy(tf,n,i):
    t,y,a=solscipy(tf,n,i)
    plt.title("Rayon en fonction du temps")
    plt.plot(t,y,'r')
    plt.show()
    plt.title("Angle en fonction du temps")
    plt.plot(t,a,'b')
    plt.savefig("GrapheOdeint")
    plt.show()
    
def traceexacte(tf,n,i):
    t,y=solexacte(tf,n,i)
    plt.title("Rayon en fonction du temps")
    plt.plot(t,y,'r')
    plt.show()
    plt.title("Angle en fonction du temps")
    plt.plot(t,'b')
    plt.savefig("GrapheExacte")
    plt.show()
    
#-------- Fonctions pour comparer graphiquement l'évolution des rayons et ----------------------
#------------- angles en fonction des diférentes CI avec une méthode------------------------- 
    
def compare1(tf,n):         #trace les rayons et angles obtenus avec Euler
    T,Y0,A0=recurrence(0,tf,n,0)
    T,Y1,A1=recurrence(0,tf,n,1)
    T,Y2,A2=recurrence(0,tf,n,2)
    T,Y3,A3=recurrence(0,tf,n,3)
    T,Y4,A4=recurrence(0,tf,n,4)
    R=[Rt for i in range(len(T))]
    trace1([T,Y0,0],[T,Y1,1],[T,Y2,2],[T,Y3,3],[T,Y4,4],[T,R,0],n)
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
    trace1([T,Y0,0],[T,Y1,1],[T,Y2,2],[T,Y3,3],[T,Y4,4],[T,R,0],n)
    plt.plot(T,A0,'r')
    plt.plot(T,A1,'g')
    plt.plot(T,A2,'b')
    plt.plot(T,A3,'y')
    plt.plot(T,A4,'k')
    plt.show()
    
def compare3(tf,n):   #trace les rayons et angles obtenus avec la solution exacte
    T,Y0,A0=solexacte(tf,n,0)
    T,Y1,A1=solexacte(tf,n,1)
    T,Y2,A3=solexacte(tf,n,2)
    T,Y3,A4=solexacte(tf,n,3)
    T,Y4,A5=solexacte(tf,n,4)
    R=[Rt]*(len(T))
    trace1([T,Y0,0],[T,Y1,1],[T,Y2,2],[T,Y3,3],[T,Y4,4],[T,R,0],n)

#--------------- Fonctions pour comparer graphiquement l'évolution des rayon et ----------------------
#---------- angle en fonction des différents pas d'intégration avec une même méthode ------------------------- 

def compare4(tf,n,i):         #trace les rayons et angles obtenus avec Euler
    T0,Y0,A0=recurrence(tf,P0*n,i)
    T1,Y1,A1=recurrence(tf,P1*n,i)
    T2,Y2,A2=recurrence(tf,P2*n,i)
    T3,Y3,A3=recurrence(tf,P3*n,i)
    T4,Y4,A4=recurrence(tf,P4*n,i)
    R=[Rt]*(len(T0))
    trace1([T0,Y0,i],[T1,Y1,i],[T2,Y2,i],[T3,Y3,i],[T4,Y4,i],[T0,R,i],n)
    plt.figure(2)    
    plt.plot(T0,A0,'r')
    plt.plot(T1,A1,'g')
    plt.plot(T2,A2,'b')
    plt.plot(T3,A3,'y')
    plt.plot(T4,A4,'k')
    plt.show()
    
def compare5(tf,n,i):         #trace les rayons et angles obtenus avec Odeint
    T0,Y0,A0=solscipy(tf,P0*n,i)
    T1,Y1,A1=solscipy(tf,P1*n,i)
    T2,Y2,A2=solscipy(tf,P2*n,i)
    T3,Y3,A3=solscipy(tf,n*P3,i)
    T4,Y4,A4=solscipy(tf,n*P4,i)
    R=[Rt]*(len(T0))
    trace1([T0,Y0,i],[T1,Y1,i],[T2,Y2,i],[T3,Y3,i],[T4,Y4,i],[T0,R,i],n)
    plt.figure(2)
    plt.plot(T0,A0,'r')
    plt.plot(T1,A1,'g')
    plt.plot(T2,A2,'b')
    plt.plot(T3,A3,'y')
    plt.plot(T4,A4,'k')
    plt.show()
    
def compare6(tf,n,i):   #trace les rayons et angles obtenus avec la solution exacte
    T0,Y0,A0=solexacte(tf,P0*n,i)
    T1,Y1,A1=solexacte(tf,P1*n,i)
    T2,Y2,A2=solexacte(tf,P2*n,i)
    T3,Y3,A3=solexacte(tf,P3*n,i)
    T4,Y4,A4=solexacte(tf,P4*n,i)
    R=[Rt]*(len(T1))
    trace1([T0,Y0,i],[T1,Y1,i],[T2,Y2,i],[T3,Y3,i],[T4,Y4,i],[T0,R,i],n)

#-------- Fonction pour comparer graphiquement les solutions ------------------
#------------- obtenues grâce aux différentes méthodes ------------------------- 

def compare7(tf,n,i): #compare les solutions obtenues avec les 3 méthodes pour une CI donnée
    T,Y0,A0=recurrence(tf,n,i)
    T1,Y1,A1=solscipy(tf,n,i)
    T2,Y2=solexacte(tf,n,i)
    R=[Rt for i in range(len(T))]
    plt.plot(T,Y0,'r',label="Euler")
    plt.plot(T2,Y2,'b',label="Solution exacte")
    plt.plot(T,R,'purple',label="Surface de la Terre")
    plt.xlabel("Temps(s)")
    plt.ylabel("Rayons(m)")
    plt.title("V(0)="+str(V[i])+"m/s")
    plt.legend()
    plt.show()
    plt.figure(2)
    plt.plot(T1,Y1,'g',label="Odeint")
    plt.plot(T2,Y2,'b',label="Solution exacte")
    plt.plot(T1,R,'purple',label="Surface de la Terre")
    plt.xlabel("Temps(s)")
    plt.ylabel("Rayons(m)")
    plt.title("V(0)="+str(V[i])+"m/s")
    plt.legend()
    plt.show()
    
#-------- Fonction pour comparer graphiquement les solutions ------------------
#------------- obtenues grâce aux différentes méthodes ------------------------- 

# j'ai effectué mes tests avec (1000,10000,solexacte,2) comme arguments. ecart4 fonctionne meme si je pense pas que ça renvoie les bons resultats. à voir      
           
def ecart(tf,n,fonction,i):
    Tref,Yref,Aref=fonction(tf,P0*n,i)
    T1,Y1,A1=fonction(tf,P1*n,i)
    T2,Y2,A2=fonction(tf,P2*n,i)
    T3,Y3,A3=fonction(tf,P3*n,i)
    T4,Y4,A4=fonction(tf,P4*n,i)
    LY=[Y1,Y2,Y3,Y4]
    LA=[A1,A2,A3,A4]
    LT=[T1,T2,T3,T4]
    Tableau=[]
    print(Tref)
    for j in range (4):  # de 0 a 3 pour parcourir LT
        for k in Tref:
            #for l in range (len(LT[j])):  solution2
            #   LT[j][l]=int(LT[j][l])     solution2
            #Ycorrespondant=LY[j][LT[j].index(int(k))]  solution2
            #Acorrespondant=LA[j][LT[j].index(int(k))]  solution2
            Ycorrespondant=LY[j][LT[j].index(k)]  # solution1 k, c'est un temps donné. LY[j]n c'est soit Y1, soit Y2 etc ... et le [LT.index[k]] correspond à la position de k dans la liste LT correspondante à celle de Y. donc on récupere la valeur de la bone liste Y qui a la meme position que la valeur k de la liste de temps qui correspond à la liste Y
            Acorrespondant=LA[j][LT[j].index(k)]  # solution1 PROBLEME LT[j] ne contient pas k exactement car k, c'est des valeurs parfaites, genre 0.0, 100.0 etc alors que LT[j] contient des valeurs "imparfaites" genre 100,000001 etc 
            Yreference=Yref[Tref.index(k)]
            Areference=Aref[Tref.index(k)]
            ecartA=Areference-Acorrespondant
            ecartY=Yreference-Ycorrespondant
            Tableau.append([k,j,ecartY]) # k : temp , j : vaut 0,1,2 ou 3 en fonction du pas, ecart est l'ecart entre ref et correspondant
# vu que celles cis fonctionnent pas; une autre solution possible mais un peut prise de tete : solution3
def ecartsol3(tf,n,fonction,i):    # ne fonctionne pas car beaucoup trop long
    Tref,Yref,Aref=fonction(tf,P0*n,i)
    T1,Y1,A1=fonction(tf,P1*n,i)
    T2,Y2,A2=fonction(tf,P2*n,i)
    T3,Y3,A3=fonction(tf,P3*n,i)
    T4,Y4,A4=fonction(tf,P4*n,i)
    LY=[Y1,Y2,Y3,Y4]
    LA=[A1,A2,A3,A4]
    LT=[T1,T2,T3,T4]   
    Tableau=[]
    imprecision=(Tref[1]-Tref[0])/2
    for j in range(4):
        for k in Tref: 
            for l in LT[j]:  # ces trois for c'est vraiment long, mais je vois pas comment faire autrement
                if l-imprecision<=k<l+imprecision: # on met le inferieur ou egal que dun cote pour pqs se retrouver avec deux if possibles
                    kcorrespondantimprecis=l
                    Ycorrespondant=LY[j][LT[j].index(l)]  
                    Acorrespondant=LA[j][LT[j].index(l)] 
                    Yreference=Yref[Tref.index(k)]
                    Areference=Aref[Tref.index(k)]
                    ecartA=Areference-Acorrespondant
                    ecartY=Yreference-Ycorrespondant
                    Tableau.append([k,j,ecartY]) 
    return (Tableau)
def ecartsol4(tf,n,fonction,i):    # tentative de reduction du temps de ecartsol3   
    Tref,Yref,Aref=fonction(tf,P0*n,i)
    T1,Y1,A1=fonction(tf,P1*n,i)
    T2,Y2,A2=fonction(tf,P2*n,i)
    T3,Y3,A3=fonction(tf,P3*n,i)
    T4,Y4,A4=fonction(tf,P4*n,i)
    P=[P1,P2,P3,P4]
    LY=[Y1,Y2,Y3,Y4]
    LA=[A1,A2,A3,A4]
    LT=[T1,T2,T3,T4]   
    Tableau=[]
    imprecision=(Tref[1]-Tref[0])/2
    # on sait que Tref contient peu de valeurs, de 0 à 1000 de 100 en 1000 donc au lieu de chercher tous les LT[j] on va chercher ceux dont la position est suseptible de s'apprecher des valeurs qu'on veut, c'est a dire celles de T ref
    for j in range(4):
        for k in Tref:
            pas=(LT[j][-1]-LT[j][0])/(P[j]*n)
            positionsusceptible=k/pas       # k, c'est la valeur qu'on veut. il faut faire positionsusceptible fois le pas pour atteindre k. comme dans la liste on se deplace de pas en pas. on aura k pour une position environ egale à positionsuseptible
            if k !=Tref[0] and k!=Tref[-1]:           
                posdepart=positionsusceptible-10
                posfin=positionsusceptible+10
            if k ==Tref[0] :           
                posdepart=0
                posfin=positionsusceptible+10
            if k==Tref[-1]:           
                posdepart=positionsusceptible-10
                posfin=len(Tref)-1
            for l in range (int(posdepart),int(posfin)):  # ces trois for c'est vraiment long, mais je vois pas comment faire autrement
                if  LT[j][l]-imprecision<=k< LT[j][l]+imprecision: # on met le inferieur ou egal que dun cote pour pqs se retrouver avec deux if possibles
                    print("ok")
                    kcorrespondantimprecis=LT[j][l]
                    Ycorrespondant=LY[j][LT[j].index(LT[j][l])]  
                    Acorrespondant=LA[j][LT[j].index(LT[j][l])] 
                    Yreference=Yref[Tref.index(k)]
                    Areference=Aref[Tref.index(k)]
                    ecartA=Areference-Acorrespondant
                    ecartY=Yreference-Ycorrespondant
                    Tableau.append([k,j,ecartY]) # je sait pas si c'est les ecartsY ou les ecartsA qu'il faut faire. pour les ecarts Y je trouve des valeurs d'ecarts qui valent environ -20000. ça parait un peu lourd. en tout cas la def arrive à se finir en une minute ou deux donc c'est bien
    return (Tableau)
