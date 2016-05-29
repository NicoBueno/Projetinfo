# -*- coding: utf-8 -*-
# Projet informatique 
# GROUPE nÂ°8 : Mendez De Maria, Bernon, Bueno
# Sujet : ETUDE DU MOUVEMENT D'UNE SONDE DANS UN CHAMP GRAVITATIONNEL

import matplotlib.pyplot as plt
import scipy.integrate as spi
import numpy as np
from math import cos 
from tkinter import *

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
    trace1([T0,A0,i],[T1,A1,i],[T2,A2,i],[T3,A3,i],[T4,A4,i],[T0,R,i],n)

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

def ecartsol(tf,n,fonction,i):      
    Tref,Yref,Aref=fonction(tf,P0*n,i) # ceci represente notre courbe de reference, avec un P0*n le plus petit
    T1,Y1,A1=fonction(tf,P1*n,i)
    T2,Y2,A2=fonction(tf,P2*n,i)
    T3,Y3,A3=fonction(tf,P3*n,i)
    T4,Y4,A4=fonction(tf,P4*n,i)
    if fonction==solscipy:
        Tiref,Yiref,Airef,Ti1,Yi1,Ai1,Ti2,Yi2,Ai2,Ti3,Yi3,Ai3,Ti4,Yi4,Ai4=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
        J=[Tiref,Yiref,Airef,Ti1,Yi1,Ai1,Ti2,Yi2,Ai2,Ti3,Yi3,Ai3,Ti4,Yi4,Ai4]
        H=[Tref,Yref,Aref,T1,Y1,A1,T2,Y2,A2,T3,Y3,A3,T4,Y4,A4]
        for h in range(len(H)):
            for g in range(len(H[h])):  
                J[h].append(H[h].tolist()[g])
        Tref,Yref,Aref,T1,Y1,A1,T2,Y2,A2,T3,Y3,A3,T4,Y4,A4=Tiref,Yiref,Airef,Ti1,Yi1,Ai1,Ti2,Yi2,Ai2,Ti3,Yi3,Ai3,Ti4,Yi4,Ai4
    P=[P1,P2,P3,P4] # pour les pas
    LY=[Y1,Y2,Y3,Y4] #liste de y
    LA=[A1,A2,A3,A4]# liste de a
    LT=[T1,T2,T3,T4] # liste de temps 
    Tableau=[]# on initialise un tableau vide
    imprecision=(Tref[1]-Tref[0])/8 # on a une liste de temps de reference Tref, et il va falloir trouver les mêmes valeurs de Tref dans les autres listes T1, T2, T3 et T4 sauf que dans ces listes y'aura pas EXACTEMENT les mêmes valeurs que pour Tref, donc on va considérer ces valeurs à un intervalle de tolérance près. ce pas, on le defini plus petit qu'un pas de temps de reference, sinon on pourrait avoir les mêmes valeurs pour deux Tref indentiques.plusieurs valeurs des listes T1,T2,T3,T4 seront dans cet intervalle, et on sélectionnera une seule  à la fin. 
    # on sait que Tref contient peu de valeurs, de 0 à 1000 de 100 en 1000 donc au lieu de chercher tous les LT[j] on va chercher ceux dont la position est suseptible de s'apprecher des valeurs qu'on veut, c'est a dire de celles de Tref
    for j in range(4): # car on teste pour 4 valeurs de pas, donc pour les listes T1,T2,T3,et T4 respectivement en positions de 0 à 3 dans la liste LT
        pas=(LT[j][1]-LT[j][0]) # on a le pas dont on cherche a quantifier l'impact sur limprecision
        for k in Tref: # pour chaque valeur de Treference
            positionsusceptible=k/pas # k, c'est la valeur de Tref qu'on veut retrouver dans les autres LT. il faut faire positionsusceptible fois le pas pour atteindre k. comme dans la liste on se deplace de pas en pas. on aura k pour une position dans la liste environ egale à positionsuseptible
            if k!=Tref[0] and k!=Tref[-1]:           # si k n'est ni 0 ni 1000, on etend notre recherche aux positions de LT entre positionsusceptible-10 et position susceptible+10
                posdepart=positionsusceptible-1 # ainsi, dans la première for bouclera 20 fois au lieu des len(LT) fois. ce qui est un enorme gain de temps
                posfin=positionsusceptible+1
            if k ==Tref[0] :           # si k =0, alors on va chercher la valeur k dans les valeurs de LT comprises aux positions allant de 0 à 10 environ
                posdepart=0
                posfin=positionsusceptible+1
            if k==Tref[-1]:           # si k=1000, on va chercher la valeur k dans les valeurs de LT comprises au positions allant de 990 à 1000 environ
                posdepart=positionsusceptible-1
                posfin=len(LT[j])
            for l in range (int(posdepart),int(posfin)): #○ on va donc parcourir nos listes LT sur les positions où on pense trouver la bonne valeur de k
                if  LT[j][l]-imprecision<=k< LT[j][l]+imprecision: # Sauf que, comme on a dit, on aura jamais exatement k, d'où l'interet de trouver une valeur à une imprecision près. on met le inferieur ou egal que dun cote pour pqs se retrouver avec deux if possibles
                    kcorrespondantimprecis=LT[j][l] # on a donc trouvé un k dans la LT qui est "egal" au k de Tref (de façon imprecise)
                    Ycorrespondant=LY[j][LT[j].index(kcorrespondantimprecis)] # le y correspondant est donc celui de la liste LY[j] qui est à la même position que la valeur de k dans sa liste respective LT[j]
                    Acorrespondant=LA[j][LT[j].index(kcorrespondantimprecis)] # pareil pour A
                    Yreference=Yref[Tref.index(k)]# le Y de reference correspondant est celui qui est à la même position dans sa liste Yref que le k dans sa liste Tref
                    Areference=Aref[Tref.index(k)]# pareil pour A
                    ecartA=Areference-Acorrespondant# on voit quel est l'ecart
                    ecartY=Yreference-Ycorrespondant
                    Tableau.append([k,j,ecartA]) # k, c'est la valeur de Tref à laquelle on a mesuré l'ecart, j ça nous dit si on a pris la liste de temps T1(j=0), T2(j=1) T3 ou T4, on peut donc en deduire le pas, et ecart A, c'est les ecart. je sait pas si c'est les ecartsY ou les ecartsA qu'il faut faire. pour les ecarts Y je trouve des valeurs d'ecarts qui valent environ -20000. ça parait un peu lourd. en tout cas la def arrive à se finir en une minute ou deux donc c'est bien
    vraiTableau=[Tableau[1]] # cette partie de code sert à recuperer une seule liste par valeur de Tref. parce qu'avec l'histoire de l'imprecision, plusieurs trucs ont etes calculées pour la même valeur de Tref. donc on s'epargne les doublons en ne recuperant qu'une seule fois chaque valeur de Tempsref
    for i in Tableau:
        if i[0]!=vraiTableau[-1][0]:
            vraiTableau.append(i)# ce tableau renvoie, pour chaque valeur de Tref, et pour chaque  j (qu'on relie au pas) , la valeur de l'imprecision (pour les Angles)
    # pour les angles j'ai l'impression que ça marche bien, pour les Y les ecartes me paraissent grand
    
    fenetre=Tk()    # TEST AFFICHAGE GRAPHIQUE AVEC tKINTER
    fenetre.geometry("1400x200")
    c=Canvas(fenetre,width=900,height=900)
    LT2=[]
    Lposcol=[]
    compteur=1
    for a in Tref:
        if a in LT2: 
            pos=LT2.index(a)+1
        else:
            pos=compteur
            compteur+=1
        LT2.append(a)
        Temps=Label(c,text=a)
        Temps.grid(row=0,column=pos)
        Lposcol.append([a,pos])
    Lj2=[]
    Lposrow=[]
    compteur=1
    for a in range(4):
        if a in Lj2: 
            pos=Lj2.index(a)+1
        else:
            pos=compteur
            compteur+=1
        Lj2.append(a)
        Pas=Label(c,text=a)
        Pas.grid(row=pos,column=0)
        Lposrow.append([a,pos])
    for a in vraiTableau:
        posrow=0
        poscol=0
        for d in Lposcol:
            if a[0]==d[0]:
                poscol=d[1]
        for b in Lposrow:
            if a[1]==b[0]:
                posrow=b[1]
        Ecart=Label(c,text=a[2])
        Ecart.grid(row=posrow,column=poscol)
    c.pack()
    fenetre.mainloop()
    fenetre.destroy()
        
 
