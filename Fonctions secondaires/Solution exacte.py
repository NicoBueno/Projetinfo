Mp=5.97*(10**24)   # masse planete
Rt=6.6356*(10**6)  # rayon terre
h=800000
G=6.67*10**(-11)
y0=Rt+h
v0=7451.9
V=[0.8*v0,0.9*v0,v0,1.1*v0,1.2*v0]
alpha=-G*Mp 

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