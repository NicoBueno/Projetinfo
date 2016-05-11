Mp=5.97*(10**24)   # masse planete
Rt=6.6356*(10**6)  # rayon terre
h=800000
G=6.67*10**(-11)
y0=Rt+h
v0=7451.9
V=[0.8*v0,0.9*v0,v0,1.1*v0,1.2*v0]
alpha=-G*Mp

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