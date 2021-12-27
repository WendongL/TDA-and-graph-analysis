import numpy as np
import scipy as sp
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import gudhi


def kcore(G, expansion_dim):
    st=gudhi.SimplexTree()
    layer = nx.algorithms.core.core_number(G)
    max_k = max(layer.values())
    for i in range(0,max_k+1):
        i_core= nx.k_core(G,k=i)
        for edge in list(i_core.edges):
            if st.find(list(edge)):
                st.remove_maximal_simplex(list(edge)) #"insert" can not update the value of filtration of a simplex already existing. The user's manual says only the high filtration value could be replaced with low value, but when I tried, neither high to low nor low to high works. So remove it and readd it with new values!
            st.insert(list(edge),filtration=i)
    st.expansion(expansion_dim)#this is not a wise way to add faces of higher dimensions, because the filtration value of faces are made equal to the maximal value of it stars. The result is that, the cycles are immediately filled with volumes when they are born. Well, maybe it is because the k-core filtration is not a good filtration...
    return st

def kcore_k_decroit(G, expansion_dim):
    st=gudhi.SimplexTree()
    layer = nx.algorithms.core.core_number(G)
    max_k = max(layer.values())
    for i in range(0,max_k+1):
        i_core= nx.k_core(G,k=i)
        for edge in list(i_core.edges):
            if st.find(list(edge)):
                st.remove_maximal_simplex(list(edge))
            st.insert(list(edge),filtration=max_k-i)
    st.expansion(expansion_dim)
    return st


#G = nx.barbell_graph(30, 10)





# G = nx.Graph()
# G.add_edge(1, 2)
# G.add_edge(1, 3)
# G.add_edge(1, 4)
# G.add_edge(2, 3)
# G.add_edge(2, 4)
# G.add_edge(4, 3)
# G.add_edge(3, 5)
# G.add_edge(4, 5)
# G.add_edge(6, 5)
#
# G.add_edge(1, 7)
# G.add_edge(7, 8)
# G.add_edge(7, 9)
# G.add_edge(7, 10)
# G.add_edge(8, 10)
# G.add_edge(9, 10)
# G.add_edge(8, 9)


# G.add_edge(2, 10)
# G.add_edge(2, 7)
# G.add_edge(1, 10)
# G.add_edge(3, 9)



# G = nx.Graph()
# G.add_edge(1, 2)
# G.add_edge(1, 3)
# G.add_edge(1, 4)
# G.add_edge(2, 3)
# G.add_edge(2, 4)
# G.add_edge(4, 3)
# G.add_edge(4, 5)
# G.add_edge(3, 5)
# G.add_edge(5, 6)

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

# G=RSG(100,1,0.5)

#G = nx.Graph()
#file=open(r'''/home/wendong/code/test2.hg''', 'r')
'''
i=0
for line in file.readlines():
    if i>0:
        (head, tail) = line.split()
        G.add_edge(int(head), int(tail))
    i+=1
G.remove_edges_from(nx.selfloop_edges(G))
'''
# G = nx.watts_strogatz_graph(100, 4, 0.5)
# G=nx.barabasi_albert_graph(100,2)
G=nx.karate_club_graph()
G2=nx.subgraph(G,[0,1,2,3,4,5,6,7,10,11,12,13,16,17,19,21])
G3=nx.subgraph(G,[8,9,14,15,18,20,22,23,24,25,26,27,28,29,30,31,32,33])

color_map = []
for node in G:
    if node in G2.nodes():
        color_map.append('lime')
    else:
        color_map.append('cyan')
nx.draw(G, node_color=color_map, with_labels=True)

plt.savefig('karate.pgf')


# diag=st.persistence(homology_coeff_field=2,min_persistence=0,persistence_dim_max=False)
# gudhi.plot_persistence_barcode(diag)
# # plt.savefig('RSG(100,1,0.5) croit bar dim1.pgf')
# gudhi.plot_persistence_diagram(diag)
# # plt.savefig('RSG(100,1,0.5) croit diag dim1.pgf')
st=kcore(G,2)
diag=st.persistence(homology_coeff_field=2,min_persistence=0,persistence_dim_max=True)
print("persistent pair:",diag)
print("betti_numbers()=",st.betti_numbers())

gudhi.plot_persistence_barcode(diag)
plt.savefig('karate_club_graph croit bar dim2.pgf')
gudhi.plot_persistence_diagram(diag)
plt.savefig('karate_club_graph croit diag dim2.pgf')

st=kcore(G2,2)
diag=st.persistence(homology_coeff_field=2,min_persistence=0,persistence_dim_max=True)
gudhi.plot_persistence_barcode(diag)
plt.savefig('karate_club_graph_1 croit bar dim2.pgf')
gudhi.plot_persistence_diagram(diag)
plt.savefig('karate_club_graph_1 croit diag dim2.pgf')

st=kcore(G3,2)
diag=st.persistence(homology_coeff_field=2,min_persistence=0,persistence_dim_max=True)
gudhi.plot_persistence_barcode(diag)
plt.savefig('karate_club_graph_2 croit bar dim2.pgf')
gudhi.plot_persistence_diagram(diag)
plt.savefig('karate_club_graph_2 croit diag dim2.pgf')


'''

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
        st.insert(list(edge),filtration=i)

    st.persistence(homology_coeff_field=2,min_persistence=0, persistence_dim_max=True) #c'est nécessaire d'écrire persistence_dim_max=True, sinon il calcule seulement H_0
    if len(list(st.betti_numbers()))<2:
        Y1.append(0)
    else:
        Y1.append(st.betti_numbers()[1])

    st=gudhi.SimplexTree()
    for edge in list(i_core.edges):
        st.insert(list(edge),filtration=i)
    st.expansion(2)
    st.persistence(homology_coeff_field=2,min_persistence=0, persistence_dim_max=True) #c'est nécessaire d'écrire persistence_dim_max=True, sinon il calcule seulement H_0
    if len(list(st.betti_numbers()))<3:
        Y2.append(0)
    else:
        Y2.append(st.betti_numbers()[2])

    st=gudhi.SimplexTree()
    for edge in list(i_core.edges):
        st.insert(list(edge),filtration=i)
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
plt.savefig('RSG(100,1,0.5) betti.pgf')
# plt.show()
'''
