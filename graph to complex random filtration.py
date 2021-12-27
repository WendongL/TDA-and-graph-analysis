import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import gudhi

plt.figure()
plt.clf()

st = gudhi.SimplexTree()
G = nx.watts_strogatz_graph(30, 13, 0.1)

nx.draw(G, with_labels=True, font_weight='bold')
layer = nx.algorithms.core.core_number(G)
max_k = max(layer.values())

for edge in list(G.edges):
    st.insert(list(edge),filtration=np.random.rand(1))
st.expansion(3) #add faces of dimension as you want, for example 3

diag=st.persistence(homology_coeff_field=2,persistence_dim_max=True)
print(diag)

print("betti_numbers()=")
print(st.betti_numbers())
gudhi.plot_persistence_diagram(diag, band=0)
plt.show()
