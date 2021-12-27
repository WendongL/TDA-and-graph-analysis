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
import numpy as np

def deg_filtre(G, expansion_dim):
    st=gudhi.SimplexTree()
    layer = nx.algorithms.core.core_number(G)
    max_k = max(layer.values())
    for i in range(0,max_k+1):
        i_core= nx.k_core(G,k=i)
        for edge in list(i_core.edges):
            if st.find(list(edge)):
                st.remove_maximal_simplex(list(edge)) #"insert" can not update the value of filtration of a simplex already existing. The user's manual says only the high filtration value could be replaced with low value, but when I tried, neither high to low nor low to high works. So remove it and readd it with new values!
            st.insert(list(edge),filtration=1/(i+1))
    st.expansion(expansion_dim)#this is not a wise way to add faces of higher dimensions, because the filtration value of faces are made equal to the maximal value of it stars. The result is that, the cycles are immediately filled with volumes when they are born. Well, maybe it is because the k-core filtration is not a good filtration...
    return st

def laplacian(G, expansion_dim, max_radius):
    #Laplacian embedding then Rips complex filtration
    #expansion_dim: the max dimension of simplices that you want to get
    #max_radius: max radius for Rips complex
    L=nx.linalg.laplacianmatrix.laplacian_matrix(G)
    L=sp.sparse.csr_matrix.toarray(L)#sparse matrices failed in small dimensions. When your G is big network, use the code in comment.
    #L= L.asfptype()
    #d,U=sp.sparse.linalg.eigs(L)
    d,U=np.linalg.eig(L)
    D=np.diag(d)
    X=U.dot(np.sqrt(D))
    X=np.nan_to_num(X)
    rips = gudhi.RipsComplex(points=X, max_edge_length=max_radius)
    st = rips.create_simplex_tree(max_dimension=expansion_dim)
    return st

def exp_ker(G, alpha, expansion_dim, max_radius):
    #Exponential diffusion kernel
    #expansion_dim: the max dimension of simplices that you want to get
    #max_radius: max radius for Rips complex
    A=nx.linalg.adjacency_matrix(G)
    A=sp.sparse.csr_matrix.toarray(A) #sparse matrix cannot have exp function

    K=np.exp(alpha*A)
    d,U=np.linalg.eig(K)
    D=np.diag(d)
    X=U.dot(np.sqrt(D))
    X=np.nan_to_num(X)
    print(X)
    rips = gudhi.RipsComplex(points=X, max_edge_length=max_radius)
    st = rips.create_simplex_tree(max_dimension=expansion_dim)
    return st

def lap_exp(G, alpha, expansion_dim, max_radius):
    #Laplacian exponential diffusion kernel
    #expansion_dim: the max dimension of simplices that you want to get
    #max_radius: max radius for Rips complex
    L=nx.linalg.laplacianmatrix.laplacian_matrix(G)
    L=sp.sparse.csr_matrix.toarray(L)#sparse matrices failed in small dimensions. When your G is big network, use the code in comment.
    #L= L.asfptype()
    #d,U=sp.sparse.linalg.eigs(L)
    K=np.exp(-alpha*L)
    d,U=np.linalg.eig(K)
    D=np.diag(d)
    X=U.dot(np.sqrt(D))
    X=np.nan_to_num(X)
    rips = gudhi.RipsComplex(points=np.transpose(X), max_edge_length=max_radius)
    st = rips.create_simplex_tree(max_dimension=expansion_dim)
    return st

def neumann(G, alpha, expansion_dim, max_radius):
    #von Neumann diffusion kernel
    #expansion_dim: the max dimension of simplices that you want to get
    #max_radius: max radius for Rips complex
    A=nx.linalg.adjacency_matrix(G)
    A=sp.sparse.csr_matrix.toarray(A) #sparse matrix cannot have exp function

    K=np.linalg.inv(np.eye(A.shape[0])-alpha*A)
    d,U=np.linalg.eig(K)
    D=np.diag(d)
    U=np.transpose(U)
    X=U.dot(np.sqrt(D))
    X=np.nan_to_num(X)
    print(X.shape)
    rips = gudhi.RipsComplex(points=X, max_edge_length=max_radius)
    st = rips.create_simplex_tree(max_dimension=expansion_dim)
    return st

def reg_lap(G, alpha, expansion_dim, max_radius):
    #Regularized Laplacian kernel
    #expansion_dim: the max dimension of simplices that you want to get
    #max_radius: max radius for Rips complex
    L=nx.linalg.laplacianmatrix.laplacian_matrix(G)
    # L=sp.sparse.csr_matrix.toarray(L)#sparse matrices failed in small dimensions. When your G is big network, use the code in comment.
    L= L.asfptype()
    #d,U=sp.sparse.linalg.eigs(L)
    K=np.linalg.inv(np.eye(L.shape[0]) + alpha*L)
    d,U=np.linalg.eig(K)
    D=np.diag(d)
    X=U.dot(np.sqrt(D))
    X=np.nan_to_num(X)
    print(X.shape)
    rips = gudhi.RipsComplex(points=X, max_edge_length=max_radius)
    st = rips.create_simplex_tree(max_dimension=expansion_dim)
    return st

def com_time(getattr, expansion_dim, max_radius):
    #Commute time kernel
    #expansion_dim: the max dimension of simplices that you want to get
    #max_radius: max radius for Rips complex
    L=nx.linalg.laplacianmatrix.laplacian_matrix(G)
    L=sp.sparse.csr_matrix.toarray(L)#sparse matrices failed in small dimensions. When your G is big network, use the code in comment.
    #L= L.asfptype()
    #d,U=sp.sparse.linalg.eigs(L)
    K=np.linalg.pinv(L)
    d,U=np.linalg.eig(K)
    D=np.diag(d)
    X=U.dot(np.sqrt(D))
    X=np.nan_to_num(X)
    rips = gudhi.RipsComplex(points=X, max_edge_length=max_radius)
    st = rips.create_simplex_tree(max_dimension=expansion_dim)
    return st

# G = nx.Graph()
# file=open(r'''E:\CODES\reseau\email-Eu\email.txt''', 'r')
# for line in file.readlines():
#     (head, tail) = line.split()
#     G.add_edge(int(head), int(tail))
# G.remove_edges_from(nx.selfloop_edges(G))

# G=nx.karate_club_graph()
# G=nx.lollipop_graph(5,10)
# G=nx.barbell_graph(5, 10, create_using=None)

G1 = nx.complete_graph(range(0, 5))
G2 = nx.complete_graph(range(5, 10))
G3 = nx.complete_graph(range(10, 15))
G4 = nx.complete_graph(range(15, 20))


H=nx.disjoint_union_all([G1,G2,G3])
H.add_edge(0, 20)
H.add_edge(20, 30)
H.add_edge(5, 21)
H.add_edge(21, 30)
H.add_edge(10, 22)
H.add_edge(22, 30)
# H.add_edge(15, 23)
# H.add_edge(23, 30)

plt.figure()
plt.clf()
nx.draw(H, with_labels=True, font_weight='bold')
plt.savefig('graphe_testé.pgf')
# st=exp_ker(G,3,2,5)
# diag=st.persistence(homology_coeff_field=2,min_persistence=0,persistence_dim_max=True)
# print("persistent pair:",diag)
# print("betti_numbers()=",st.betti_numbers())
#
# gudhi.plot_persistence_barcode(diag)
# gudhi.plot_persistence_diagram(diag)
# plt.show()