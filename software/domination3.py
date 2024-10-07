def p_lambda_seq(B, la, seq):
    m = B.nrows()
    rk = B.rank()
    n = B.ncols()
    B = block_matrix([[B,matrix(la).transpose()]])
    for k in seq:
        B.mutate(k)
    B, la  = B[:,:-1], B[:,-1] 
    Ps = [Polyhedron(rays=B.columns(),base_ring=QQ).translation(la)]
    for k in reversed(seq):
        Ep = E(B, k, 1)
        Em = E(B, k, -1)
        Hp = Polyhedron(ieqs=[(0,)*(k+1)+(1,)+(0,)*(m-k-1)])
        Hm = Polyhedron(ieqs=[(0,)*(k+1)+(-1,)+(0,)*(m-k-1)])
        new_Ps = []
        for P in Ps:
            Pp = P.intersection(Hp)
            if Pp.dimension() == rk:
                new_Ps.append(Ep*Pp)
            Pm = P.intersection(Hm)
            if Pm.dimension() == rk:
                new_Ps.append(Em*Pm)
        Ps = new_Ps
        B.mutate(k)
#    return Ps
    return PolyhedralComplex(Ps).union_as_polyhedron()

def p_lambda_seq_as_list_any_rank(B, la, seq):
    m = B.nrows()
    rk = B.rank()
    n = B.ncols()
    B = block_matrix([[B,matrix(la).transpose()]])
    for k in seq:
        B.mutate(k)
    B, la  = B[:,:-1], B[:,-1] 
    Ps = [Polyhedron(rays=B.columns(),base_ring=QQ).translation(la)]
    for k in reversed(seq):
        Ep = E(B, k, 1)
        Em = E(B, k, -1)
        Hp = Polyhedron(ieqs=[(0,)*(k+1)+(1,)+(0,)*(m-k-1)])
        Hm = Polyhedron(ieqs=[(0,)*(k+1)+(-1,)+(0,)*(m-k-1)])
        new_Ps = []
        for P in Ps:
            Pp = P.intersection(Hp)
            if Pp.dimension() == rk:
                new_Ps.append(Ep*Pp)
            Pm = P.intersection(Hm)
            if Pm.dimension() == rk:
                new_Ps.append(Em*Pm)
        Ps = new_Ps
        B.mutate(k)
    P = Ps[0]
    for Q in Ps:
        P = P.convex_hull(Q)
    return P

def p_lambda_seq_convexhull(B, la, seq):
    m = B.nrows()
    rk = B.rank()
    n = B.ncols()
    B = block_matrix([[B,matrix(la).transpose()]])
    for k in seq:
        B.mutate(k)
    B, la  = B[:,:-1], B[:,-1] 
    Ps = [Polyhedron(rays=B.columns(),base_ring=QQ).translation(la)]
    for k in reversed(seq):
        Ep = E(B, k, 1)
        Em = E(B, k, -1)
        Hp = Polyhedron(ieqs=[(0,)*(k+1)+(1,)+(0,)*(m-k-1)])
        Hm = Polyhedron(ieqs=[(0,)*(k+1)+(-1,)+(0,)*(m-k-1)])
        new_Ps = []
        for P in Ps:
            Pp = P.intersection(Hp)
            if Pp.dimension() == rk:
                new_Ps.append(Ep*Pp)
            Pm = P.intersection(Hm)
            if Pm.dimension() == rk:
                new_Ps.append(Em*Pm)
        Ps = new_Ps
        B.mutate(k)
    P = Ps[0]
    for Q in Ps:
        P = P.convex_hull(Q)
    return P


def p_lambda_seq_as_list(B, la, seq):
    m = B.nrows()
    n = B.ncols()
    B = block_matrix([[B,matrix(la).transpose()]])
    for k in seq:
        B.mutate(k)
    B, la  = B[:,:-1], B[:,-1] 
    Ps = [Polyhedron(rays=B.columns()).translation(la)]
    for k in reversed(seq):
        Ep = E(B, k, 1)
        Em = E(B, k, -1)
        Hp = Polyhedron(ieqs=[(0,)*(k+1)+(1,)+(0,)*(m-k-1)])
        Hm = Polyhedron(ieqs=[(0,)*(k+1)+(-1,)+(0,)*(m-k-1)])
        new_Ps = []
        for P in Ps:
            Pp = P.intersection(Hp)
            if Pp.dimension() == n:
                new_Ps.append(Ep*Pp)
            Pm = P.intersection(Hm)
            if Pm.dimension() == n:
                new_Ps.append(Em*Pm)
        Ps = new_Ps
        B.mutate(k)
    return Ps

def intersect_polytope_lists(A,B):
    return [ P for Q in A for R in B for P in [Q.intersection(R)] if not P.is_empty() ]

def E(B, k, eps):
    n = B.nrows()
    M = identity_matrix(n)
    for i in range(n):
        M[i,k] = max(eps*B[i,k], 0)
    M[k,k] = -1
    return M
