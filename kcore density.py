#kcore et densité

import numpy as np
import scipy as sp
import networkx as nx
import matplotlib.pyplot as plt
import gudhi
import numpy as np
import itertools


#file=open(r'''E:\CODES\reseau\CA-AstroPh.txt''', 'r')



# file=open(r'''E:\CODES\reseau\email-Eu\email.txt''', 'r')
# i=0
# for line in file.readlines():
#     i+=1
#     if i>4:
#         (head, tail) = line.split()
#         G.add_edge(int(head), int(tail))
# G.remove_edges_from(nx.selfloop_edges(G))

G=nx.karate_club_graph()


layer = nx.algorithms.core.core_number(G)
max_k = max(layer.values())

D = []
for i in range(0,max_k+1):
    i_core = nx.k_core(G,k=i)
    num_edge= len(list(i_core.edges()))
    num_node= len(list(i_core.nodes()))
    d = num_edge/num_node
    D.append(d)

plt.figure(10)
plt.clf()
plt.plot(np.arange(0,max_k+1),D,label='k-core')

st=gudhi.SimplexTree()
layer = nx.algorithms.core.core_number(G)
max_k = max(layer.values())

C=[] #array contenant les sommets de chaque étape de k-core
for i in range(0,max_k+1):
    C.append([])
for i in range(0,max_k+1):
    i_core= nx.k_core(G,k=i)
    for edge in list(i_core.edges):
        if st.find(list(edge)):
            st.remove_maximal_simplex(list(edge)) #"insert" can not update the value of filtration of a simplex already existing. The user's manual says only the high filtration value could be replaced with low value, but when I tried, neither high to low nor low to high works. So remove it and readd it with new values!
        st.insert(list(edge),filtration=i)
    # st.insert(list(i_core.nodes()), filtration = i)
    print(i)
    if i ==0:
        for vert in G.nodes:
            st.insert([vert], filtration = i)
            C[0].append(vert)
    else:
        C[i]=list(i_core.nodes())
        print(i_core.nodes())


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


D=[] #nodes in new filtration
for i in range(0,max_k):
    D.append([])
D.append(C[max_k]) #la dernière case
D2=np.zeros(max_k+1) #densité
for i in reversed(range(0,max_k+1)):
    # st.expansion(i+1)
    # simplices=[]
    # fil = st.get_filtration()
    # for x in fil:
    #     simplices.append(x[0])
    # shell= set(nx.k_core(G,k=i-1).edges())-set(nx.k_core(G,k=i).edges())
    # shell_nodes=[]
    # for edge in shell:
    #     shell_nodes.append(edge[0])
    #     shell_nodes.append(edge[1])

    # border_nodes=list(set(shell_nodes)-set(nx.k_core(G,k=i).nodes()))
    # for sigma in simplices:
    #     if len(sigma)>=i-1 and set(sigma).issubset(shell_nodes) and len(intersection(sigma, border_nodes))>=2 and i>2:
    #         for v in sigma:
    #             D[i].append(v)

    #List of 2-simplices
    triangles = list(set(itertools.combinations(C[i], 3)))
    for triangle in triangles:
        # if st.find(triangle):
        #     st.remove_maximal_simplex(triangle)
        st.insert(triangle, filtration=i)

    #List of 3-simplices
    tetras = list(set(itertools.combinations(C[i], 4)))
    for tetra in tetras:
        # if st.find(tetra):
        #     st.remove_maximal_simplex(tetra)
        st.insert(tetra, filtration =i)

for i in reversed(range(0,max_k+1)):
    fil=st.get_filtration()
    for x in fil:
        if len(x[0])>= i-1 and x[1]<=i-1:
            for v in x[0]:
                if not v in D[i]:
                    D[i].append(v)
    if i> 0 : D[i-1]=list(set(D[i]).union(set(C[i-1])))
    S=nx.subgraph(G,D[i])
    num_edge= len(S.edges)
    num_node= len(list(D[i]))
    # plt.figure(i)
    # plt.clf()
    # nx.draw(S, with_labels=True, font_weight='bold')
    D2[i] = num_edge/num_node
    print(set(C[i]).issubset(set(D[i])))
    if i < max_k: print(set(D[i+1]).issubset(set(C[i])))

plt.figure(10)
plt.plot(np.arange(0,max_k+1),D2,label='k-core density friendly')





plt.title('Densité suivant k')
plt.legend(loc='best')

plt.show()



