# inequalities are stored in the same format as Polyhedron expect them:
# [-1,7,3,4] represents the inequality 7x_1+3x_2+4x_3>= 1

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
            ieqs = [ (const,) + tuple(E(sgn).transpose()*vector(normal_p)) for (const, *normal_p) in ieqs_p ] + [ ((0,)*(k+1) + (sgn,) + (0,)*(2*n-k-1)) ]
            #check the sign of sgn above
            P = Polyhedron(ieqs=ieqs)
            if P.dimension() == 2*n:
                ieqs_list.append( P.inequalities_list() )
    return ieqs_list


def get_inequalities(B, seq, la):
    #this function assumes B is 2n x n and the bottom is invertible and la is a vector of length 2n

    n = B.ncols()

    M = block_matrix([[B,la.column()]])
    for k in seq:
        M.mutate(k)

    la_p = vector(M[:,-1])

    return mutate_inequalities(B, seq, [ [ (-vector(normal)*la_p[n:],) + (0,)*n+tuple(normal)  for normal in M[n:,:n].inverse().rows() ] ])


def parabolic_sequences(B, seq, k):
    A = ClusterAlgebra(B)
    A.reset_current_seed()
    A.current_seed().mutate(seq, mutating_F=False)
    seeds = A.seeds(mutating_F=False, allowed_directions=[i for i in range(A.rank()) if i != k], from_current_seed=True)
    return [ seed.path_from_initial_seed() for seed in seeds ]


def get_polyhedron(B, seq, la, projection=None):
    #this function assumes B is 2n x n and the bottom is invertible and la is a vector of length 2n

    ieqs_list = get_inequalities(B, seq, la)
    dim = B.nrows()
    if projection:
        dim = B.ncols()
        ieqs_list = project_inequalities(B, la, ieqs_list, projection=projection)
    polyhedra = [ Polyhedron(ieqs=ieqs) for ieqs in ieqs_list ]
    Q = Polyhedron(ambient_dim=dim,base_ring=QQ)
    for P in polyhedra:
        Q = Q.convex_hull(P)
    return Q


def project_inequalities(B, la, ieqs_list, projection=None):
    # TODO: this function still uses the old notation for inequalities
    # TODO: the g projection needs a shift
    #this function assumes B is 2n x n and the bottom is invertible
    if not projection:
        return ieqs_list

    n = B.ncols()
    top = B[:n,:n]
    bottom = B[n:,:n]
    if projection == 'c':
        elimination_matrix = block_matrix([[-identity_matrix(n), top*bottom.inverse()]])
        start = 0
        end = n
        shift = lambda normal: vector(normal).dot_product(la)
    if projection == 'g':
        #requires top to be invertible
        bti = bottom*top.inverse()
        elimination_matrix = block_matrix([[bti, -identity_matrix(n)]])
        start = n
        end = 2*n
        shift = lambda normal: -vector(normal[n:])*bti*la[:n]
    new_ieqs_list = []
    for ieqs in ieqs_list:
        #take first entries in one and last entries in the other
        new_ieqs = [ (const+shift(normal),)+tuple(vector(normal)+vector(normal[start:end])*elimination_matrix)[2*n-end:2*n-start] for (const, *normal) in ieqs ]
        new_ieqs_list.append(new_ieqs)

    return new_ieqs_list


def parabolic_intersection(B, seq, la, k, projection=None):
    n = B.ncols()
    if B.is_square():
        B = B.stack(identity_matrix(n))

    if len(la) == n:
        la = tuple(la)+(0,)*n
    if len(la) != B.nrows():
        raise ValueError("la has the wrong length")
    la = vector(la)

    seqs = parabolic_sequences(B, seq, k)
    seq = seqs.pop(0)

    P = get_polyhedron(B, seq, la, projection=projection)
    for seq in seqs:
        P = P.intersection(get_polyhedron(B, seq, la, projection=projection))

    return P


def E(B,k,sgn):
     n = B.nrows()
     E = identity_matrix(n)
     E.set_column(k, list(map( lambda b: max(sgn*b,0), B.column(k))))
     E[k,k] = -1
     return E


#hard coded for now
#b_list = [[b,-1],[c,-1]]
#var('a','b','c','ap','bp','cp')
#b_list = [[ap,bp,-1],[cp,a,-1],[b,c,-1]]
n=4
var(['b'+str(i)+str(j) for i in range(1,n) for j in range(1,n+1)])
b_list = [ [ eval('b'+str(i)+str(j)) for i in range(1,n) ] + [-1] for j in range(1,n+1) ]
def gen_chebyshev(k,m,eps,n=len(b_list)):
    if eps == 0:
        eps = n
    if m < 1 or m > n:
        raise ValueError("m must be between 1 and " + str(n))
    if eps <1 or eps > n:
        raise ValueError("eps must be between 1 and " + str(n))
    if k in range(-n+1,0):
        return 0
    if k == 0:
        return 1
    if k == m:
        return b_list[eps-1][k-1]
    return sum(gen_chebyshev(l,l,eps)*gen_chebyshev(k-l,1,(eps+l)%n) for l in range(m,n+1))

leading_terms = {}
chebs = { (k,m,eps): u for k in range(0,9) for m in range(1,n+1) for eps in range(1,n+1) for u in [expand(gen_chebyshev(k,m,eps))] if u != 0 }
for key1 in chebs:
    for key2 in chebs:
        if key1 > key2:
            continue
        prod = expand(chebs[key1]*chebs[key2])
        if prod in ZZ:
            continue
        leading_term = prod.polynomial(QQ).monomials()[0]
        if leading_term not in leading_terms:
            leading_terms[leading_term] = [(key1,key2)]
        else:
            leading_terms[leading_term].append((key1,key2))

def expansion(poly, decomp):
    if poly == 0:
        return [decomp]
    out = []
    mons = poly.polynomial(QQ).monomials()
    leading_term = mons[0]
    if leading_term not in leading_terms:
        return -1
    for (key1,key2) in leading_terms[leading_term]:
        if key1[1] == n and key2[1] == n:
            continue
        new_poly = expand(poly - gen_chebyshev(*key1)*gen_chebyshev(*key2))
        if new_poly == 0 or (new_poly.polynomial(QQ).monomials()[0] != leading_term and all(mon in mons for mon in new_poly.polynomial(QQ).monomials())):
            exp = expansion(new_poly,decomp+[(key1,key2)])
            if exp != -1:
                out += exp
    return out
