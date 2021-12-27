import matplotlib
import matplotlib.pyplot as plt
import gudhi
import numpy as np
import scipy as sp
import networkx as nx
import itertools

def kcore_decroit_voisin(G, max_dim):
    st=gudhi.SimplexTree()
    layer = nx.algorithms.core.core_number(G)
    max_k = max(layer.values())
    for i in range(0,max_k+1):
        i_core= nx.k_core(G,k=i)
        for v in i_core.nodes():
            voisins=list(i_core.neighbors(v))
            if len(voisins)<=max_dim+1:
                if st.find(voisins):
                    st.remove_maximal_simplex(voisins)
                st.insert(voisins, filtration=i)
            else:
                for simplex in itertools.combinations(voisins, max_dim+1):
                    if st.find(voisins):
                        st.remove_maximal_simplex(voisins)
                    st.insert(simplex, filtration=i)
    return st

def kcore_croit_voisin(G, max_dim):
    st=gudhi.SimplexTree()
    layer = nx.algorithms.core.core_number(G)
    max_k = max(layer.values())
    for i in range(0,max_k+1):
        i_core= nx.k_core(G,k=i)
        for v in i_core.nodes():
            voisins=list(i_core.neighbors(v))
            if len(voisins)<=max_dim+1:
                if st.find(voisins):
                    st.remove_maximal_simplex(voisins)
                st.insert(voisins, filtration=i)
            else:
                for simplex in itertools.combinations(voisins, max_dim+1):
                    if st.find(voisins):
                        st.remove_maximal_simplex(voisins)
                    st.insert(simplex, filtration=max_k-i)
    return st

G = nx.Graph()
G.add_edge(1, 2)
G.add_edge(1, 3)
G.add_edge(1, 4)
G.add_edge(2, 3)
G.add_edge(2, 4)
G.add_edge(4, 3)
G.add_edge(3, 5)
G.add_edge(4, 5)
G.add_edge(6, 5)

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

# G=RSG(50,1,0.5)
G=nx.barabasi_albert_graph(20,2)

nx.draw(G, with_labels=True)
st=kcore_decroit_voisin(G,2)
diag=st.persistence(homology_coeff_field=2,min_persistence=0,persistence_dim_max=True)
print("persistent pair:",diag)
print("betti_numbers()=",st.betti_numbers())

gudhi.plot_persistence_barcode(diag)
gudhi.plot_persistence_diagram(diag)

plt.show()