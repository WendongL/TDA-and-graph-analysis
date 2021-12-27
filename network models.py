import numpy as np
import scipy as sp
import networkx as nx
import matplotlib.pyplot as plt
import gudhi
import numpy as np

def RSG(N,r,s): #sur S^2
    lat= (np.random.rand(N)-0.5)*np.pi
    lon= np.random.rand(N)*np.pi*2
    G=nx.Graph()
    for i in range(N):
        for j in range(i):
            d=2*r*np.arcsin(np.sqrt(np.sin((lat[j]-lat[i])/2)**2 + np.cos(lat[i])* np.cos(lat[j]) * np.sin((lon[j]- lon[i])/2)**2))
            if d<=s:
                G.add_edge(i,j)
    return G

G=RSG(100,1,0.5)
nx.draw(G, with_labels=True, font_weight='bold')
plt.show()