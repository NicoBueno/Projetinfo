Mp=5.97*(10**24)   # masse planete
Rt=6.6356*(10**6)  # rayon terre
h=800000
G=6.67*10**(-11)
y0=Rt+h
v0=7451.9
V=[0.8*v0,0.9*v0,v0,1.1*v0,1.2*v0]
alpha=-G*Mp

def solscipy(f,tf,n):           #on trace la solution à l'equation differentielle grace a scipy
    ts=np.linspace(0.,tf,n)
    y=spi.odeint(f,tf,ts)
    Y=y[:,0]
    plt.plot(ts,Y,'g')