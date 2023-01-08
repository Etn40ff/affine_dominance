from collections import defaultdict


def mutate_inequalities(B, seq, ieqs_list):
    #this function assumes B is 2n x n
    #this pulls the inequalities back along the sequence, starting from the end

    if seq == []:
        return ieqs_list

    n = B.ncols()
    k = seq[0]

    def E(sgn):
        E = identity_matrix(2*n)
        E.set_column(k, list(map( lambda b: max(sgn*b,0), B.column(k))))
        E[k,k] = -1
        return E

    B_p = copy(B)
    B_p.mutate(k)

    ieqs_list_p = mutate_inequalities(B_p, seq[1:], ieqs_list)

    ieqs_list = []
    for ieqs_p in ieqs_list_p:
        for sgn in [-1, 1]:
            ieqs = [ (tuple(E(sgn).transpose()*vector(normal_p)), const) for (normal_p, const) in ieqs_p ] + [ ((0,)*(k) + (sgn,) + (0,)*(2*n-k-1),0) ]
            #check the sign of sgn above
            P = Polyhedron(ieqs=[ [-rhs]+list(lhs) for (lhs,rhs) in ieqs ] )
            if P.dimension() == 2*n:
                ieqs_list.append( [ (v[1:],-v[0]) for v in P.inequalities_list() ] )
    return ieqs_list


def get_inequalities(B, seq, la):
    #this function assumes B is 2n x n and the bottom is invertible and la is a vector of length 2n

    n = B.ncols()

    M = block_matrix([[B,la.column()]])
    for k in seq:
        M.mutate(k)

    la_p = vector(M[:,-1])

    return mutate_inequalities(B, seq, [ [ ((0,)*n+tuple(normal), vector(normal)*la_p[n:]) for normal in M[n:,:n].inverse().rows() ] ])


def parabolic_sequences(B, seq, k):
    A = ClusterAlgebra(B)
    A.current_seed().mutate(seq, mutating_F=False)
    seeds = A.seeds(mutating_F=False, allowed_directions=[i for i in range(A.rank()) if i != k], from_current_seed=True)
    return [ seed.path_from_initial_seed() for seed in seeds ]


def get_polyhedron(B, seq, la, projection=None):
    #this function assumes B is 2n x n and the bottom is invertible and la is a vector of length 2n
    ieqs_list = get_inequalities(B, seq, la)
    dim = B.nrows()
    if projection:
        dim = B.ncols()
        ieqs_list = project_inequalities(B, ieqs_list, projection=projection)
    polyhedra = [ Polyhedron(ieqs=[ [-rhs]+list(lhs) for (lhs,rhs) in ieqs ] ) for ieqs in ieqs_list ]
    Q = Polyhedron(ambient_dim=dim,base_ring=QQ)
    for P in polyhedra:
        Q = Q.convex_hull(P)
    return Q


def project_inequalities(B, ieqs_list, projection='c'):
    #this function assumes B is 2n x n and the bottom is invertible
    n = B.ncols()
    top = B[:n,:n]
    bottom = B[n:,:n]
    if projection == 'c':
        elimination_matrix = block_matrix([[-identity_matrix(n), top*bottom.inverse()]])
        start = 0
        end = n
    if projection == 'g':
        #requires top to be invertible
        elimination_matrix = block_matrix([[bottom*top.inverse(), -identity_matrix(n)]])
        start = n
        end = 2*n
    new_ieqs_list = []
    for ieqs in ieqs_list:
        new_ieqs = []
        for ieq in ieqs:
            #take first entries in one and last entries in the other
            new_ieqs.append((tuple(vector(ieq[0])-vector(ieq[0][start:end])*elimination_matrix)[2*n-end:2*n-start],ieq[1]))
        new_ieqs_list.append(new_ieqs)

    return new_ieqs_list


def global_function(B, seq, la):
    n = B.ncols()
    if B.is_square():
        B = B.stack(identity_matrix(n))

    if len(la) == n:
        la = tuple(la)+(0,)*n
    if len(la) != B.nrows():
        raise ValueError("la has the wrong length")
    la = vector(la)
