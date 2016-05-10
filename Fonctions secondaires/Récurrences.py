# -*- coding: utf-8 -*-
Mp=5.97*(10**24)   # masse planete
Rt=6.6356*(10**6)  # rayon terre
h=800000
G=9.81
y0=Rt+h
v0=7451.9
V=[0.8*v0,0.9*v0,v0,1.1*v0,1.2*v0]
alpha=-G*Mp 



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