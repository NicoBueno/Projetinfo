def trace(T,Y,Z):
    t,x,y,a,b=recurrence(ti,tf,n,y0,z0,a0,b0,i)  # à compléter
    plt.title("Angle et rayon en fonction du temps")
    plt.plot(x,t,'r',y,t,'b')     #on trace les solutions obtenues grâce à python