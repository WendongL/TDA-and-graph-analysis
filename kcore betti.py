import numpy as np
import scipy as sp
import networkx as nx
import matplotlib.pyplot as plt
import gudhi


G=nx.karate_club_graph()
G = nx.barbell_graph(30, 10)
G = nx.watts_strogatz_graph(30, 13, 0.1)

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

G.add_edge(1, 7)
G.add_edge(7, 8)
G.add_edge(7, 9)
G.add_edge(7, 10)
G.add_edge(8, 10)
G.add_edge(9, 10)
G.add_edge(8, 9)
G.add_edge(2, 10)
G.add_edge(2, 7)
G.add_edge(1, 10)
G.add_edge(3, 9)

'''
from networkx.generators.degree_seq import expected_degree_graph

# make a random graph of 500 nodes with expected degrees of 50
n = 500  # n nodes
p = 0.1
w = [p * n for i in range(n)]  # w = p*n for all nodes
G = expected_degree_graph(w)  # configuration model
G.remove_edges_from(nx.selfloop_edges(G))

z = [5, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1]
G = nx.configuration_model(z)  # configuration model

plt.figure()
plt.clf()
nx.draw(G, with_labels=True, font_weight='bold')
plt.show()'''
#I wanted to use expected_degree_graph, which takes the degrees and generates a graph. It is useful but it creates a multigraph, which makes impossible to create a corresponding SimplexTree.


layer = nx.algorithms.core.core_number(G)
max_k = max(layer.values())
N=list(set(layer.values())) #les nombres k existants de k-core
X=[] #nombre de sommets dans chaque core
Y1=[] #nombre de cycles dans chaque core
Y2=[] #nombre de voids dans chaque core
Y3=[] #nombre de 4-dimensional voids dans chaque core

for i in N:
    i_core= nx.k_core(G,k=i)
    X.append(len(list(i_core.nodes)))
    st=gudhi.SimplexTree()
    for edge in list(i_core.edges):
        st.insert(list(edge),filtration=0)

    st.persistence(homology_coeff_field=2,min_persistence=0, persistence_dim_max=True) #c'est nécessaire d'écrire persistence_dim_max=True, sinon il calcule seulement H_0
    if len(list(st.betti_numbers()))<2:
        Y1.append(0)
    else:
        Y1.append(st.betti_numbers()[1])

    st=gudhi.SimplexTree() #Il semble que st.persistence ne peut pas être mis à jours. J'ai donc refait le st.
    for edge in list(i_core.edges):
        st.insert(list(edge),filtration=0)
    st.expansion(2)
    st.persistence(homology_coeff_field=2,min_persistence=0, persistence_dim_max=True) #c'est nécessaire d'écrire persistence_dim_max=True, sinon il calcule seulement H_0
    if len(list(st.betti_numbers()))<3:
        Y2.append(0)
    else:
        Y2.append(st.betti_numbers()[2])


    st=gudhi.SimplexTree() #Il semble que st.persistence ne peut pas être mis à jours. J'ai donc refait le st.
    for edge in list(i_core.edges):
        st.insert(list(edge),filtration=0)
    st.expansion(3)
    diagram=st.persistence(homology_coeff_field=2,min_persistence=0, persistence_dim_max=True) #c'est nécessaire d'écrire persistence_dim_max=True, sinon il calcule seulement H_0
    if len(list(st.betti_numbers()))<4:
        Y3.append(0)
    else:
        Y3.append(st.betti_numbers()[3])


X=np.array(X)
Y1=np.array(Y1)
print(X-Y1)
Y2=np.array(Y2)
print(X-Y2)
Y3=np.array(Y3)
print(X-Y3)


plt.figure()
plt.clf()
plt.plot(N,X,label='nombre de sommets') #straight line which is distinguished from Y1. They sometimes overlap.
plt.plot(N,Y1,'o',label='H1')
plt.plot(N,Y2,'o',label='H2')
plt.plot(N,Y3,'o',label='H3')
plt.legend(loc='best')


plt.figure()
plt.clf()
nx.draw(G, with_labels=True, font_weight='bold')

plt.show()