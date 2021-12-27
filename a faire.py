import numpy as np
import scipy as sp
import networkx as nx
import matplotlib.pyplot as plt
import gudhi

def kcore(G, expansion_dim):
    st=gudhi.SimplexTree()
    for edge in list(G.edges):
        st.insert(list(edge))
    return st


G = nx.Graph()
G.add_edge(1, 2)
G.add_edge(1, 3)
G.add_edge(1, 4)
G.add_edge(2, 3)
G.add_edge(2, 4)
G.add_edge(4, 3)
G.add_edge(4, 5)
G.add_edge(3, 5)
G.add_edge(5, 6)

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

plt.figure()
nx.draw(G, with_labels=True, font_weight='bold')


st=kcore(G,1)
diag=st.persistence(homology_coeff_field=2,min_persistence=0,persistence_dim_max=True)
print("persistent pair:",diag)
print("betti_numbers()=",st.betti_numbers())

gudhi.plot_persistence_barcode(diag)
gudhi.plot_persistence_diagram(diag)
plt.show()
