class DominanceOrder(SageObject):
    r"""
    A class to compute the dominace order as defined in `arXiv:1902.09507`
    """

    def __init__(self, B):
        #B is n x n
        self.n = B.ncols()
        self.B = copy(B).stack(identity_matrix(self.n))
        self.A = ClusterAlgebra(self.B, principal_coefficients=True)
        self.Aop = ClusterAlgebra((-B).stack(identity_matrix(self.n)), principal_coefficients=True)

    def paths_up_to_length(self, k):
        paths = [ s.path_from_initial_seed() for s in self.A.seeds(mutating_F=False, depth=k) ]
        prefixes = [ p[:-1] for p in paths ]
        return [ p for p in paths if p not in prefixes ]

    def cut_along_sequence(self, g, B, seq):
        r"""
        Return the intersection of the cones dominated by all the translates of g along the sequence of mutations seq
        """
        if len(g) == self.n:
            g = tuple(g) + tuple([0 for _ in range(self.n)])
        g = vector(g)
        current_cone = Polyhedron(rays=B.columns(),base_ring=QQ).translation(g)
        if seq == []:
            return current_cone
        k = seq.pop()
        Hp = Polyhedron(ieqs=[(0,)*(k+1)+(1,)+(0,)*(2*self.n-k-1)])
        Ep = matrix(2*self.n, lambda i,j: (1 if i == j else 0) if j != k else (max(B[i,k],0) if i != k else -1) )
        Hm = Polyhedron(ieqs=[(0,)*(k+1)+(-1,)+(0,)*(2*self.n-k-1)])
        Em = matrix(2*self.n, lambda i,j: (1 if i == j else 0) if j != k else (max(-B[i,k],0) if i != k else -1) )
        new_g = (Ep if g in Hp else Em) * g
        new_B = copy(B)
        new_B.mutate(k)
        polytope = self.cut_along_sequence(new_g, new_B, seq)
        polytope_p = Em*(polytope.intersection(Hp))
        polytope_m = Ep*(polytope.intersection(Hm))
        polytope = polytope_p.convex_hull(polytope_m)
        return polytope
        return polytope.intersection(current_cone)

    def cut_along_sequence(self, g, B, seq, return_steps=1000000):
        r"""
        Return the intersection of the cones dominated by all the translates of g along the sequence of mutations seq
        """
        if len(g) == self.n:
            g = tuple(g) + tuple([0 for _ in range(self.n)])
        g = vector(g)
        current_cone = Polyhedron(rays=B.columns(),base_ring=QQ).translation(g)
        if seq == []:
            return return_steps, current_cone
        k = seq.pop()
        Hp = Polyhedron(ieqs=[(0,)*(k+1)+(1,)+(0,)*(2*self.n-k-1)])
        Ep = matrix(2*self.n, lambda i,j: (1 if i == j else 0) if j != k else (max(B[i,k],0) if i != k else -1) )
        Hm = Polyhedron(ieqs=[(0,)*(k+1)+(-1,)+(0,)*(2*self.n-k-1)])
        Em = matrix(2*self.n, lambda i,j: (1 if i == j else 0) if j != k else (max(-B[i,k],0) if i != k else -1) )
        new_g = (Ep if g in Hp else Em) * g
        new_B = copy(B)
        new_B.mutate(k)
        remaining_steps, polytope = self.cut_along_sequence(new_g, new_B, seq, return_steps)
        if remaining_steps > 0:
            polytope_p = Em*(polytope.intersection(Hp))
            polytope_m = Ep*(polytope.intersection(Hm))
            polytope = polytope_p.convex_hull(polytope_m)
            remaining_steps -= 1
        return remaining_steps, polytope
        return polytope.intersection(current_cone)

    @cached_method(key=lambda self, g, depth: (tuple(g), depth))
    def dominated_polytope(self, g, depth):
        paths = self.paths_up_to_length(depth)
        p = paths.pop()
        polytope = self.cut_along_sequence(g,self.B,p)
        while paths:
            p = paths.pop()
            polytope = polytope.intersection(self.cut_along_sequence(g,self.B,p))
        return polytope

    def dominated_g_vectors(self, g, depth):
        polytope = self.dominated_polytope(g, depth)
        pts = polytope.integral_points()
        lattice = Polyhedron(rays=self.B.columns(),base_ring=QQ)
        g = vector(g)
        return [ p for p in pts if p-g in lattice ]

    def show_domination(self, g, depth):
        g = vector(g)
        pts = self.dominated_g_vectors(g, depth)
        polytope = self.dominated_polytope(g, depth)
        return polytope.plot(zorder=-1) + list_plot(pts, color='red', size=50) 

    def generic_F_polynomial(self, g):
        dominated = self.dominated_g_vectors(g, 10)
        dominated.remove(vector(g))
        u = self.A._U.gens()
        return self.A.theta_basis_F_polynomial(g)+sum( var('r%d'%i) * self.A.theta_basis_F_polynomial(h) * u[0]**((g[1]-h[1])/self.B[0,1]) * u[1]**((g[0]-h[0])/self.B[1,0]) for (i,h) in enumerate(dominated)) 

    def reorder(self, f):
        u = self.A._U.gens()
        if self.B[0,1]<0:
            u = list(reversed(u))
        coefficients = f.collect(u[0]).coefficients(u[0])
        collect = list(map(lambda x: [x[0].collect(u[1]),x[1]], coefficients))
        return sum( x[0]*u[0]**x[1] for x in collect)  

    def plot_F_polynomial(self, f):
        u = self.A._U.gens()
        vrs = f.variables()[:-2]
        R = QQ[vrs][u]
        f = SR(f).polynomial(ring=R)
        G = Graphics()
        for p,cf in f.dict().items():
            G += text(cf, p)
        G.axes(show=False)
        return G

    def coxeter_matrix(self):
        return -self.A.euler_matrix()*self.Aop.euler_matrix().inverse()

    def eigenvectors(self):
        return self.coxeter_matrix().eigenvectors_right()

    def opposite_g_vectors(self, g):
        g = vector(g)
        return [ (1+self.B*self.A.euler_matrix().inverse())*g, (1-self.B*self.Aop.euler_matrix().inverse())*g ]

    def find_opposit_corner(self, g):
        g = vector(g)
        C = 2-self.B.apply_map(abs)
        return C.adjugate().transpose()*g
            



def cut_along_sequence(g, B, seq):
    r"""
    Return the intersection of the cones dominated by all the translates of g along the sequence of mutations seq
    """
    current_cone = Polyhedron(rays=B.columns(),base_ring=QQ).translation(g)
    if seq == []:
        return current_cone
    k = seq.pop()
    n = B.ncols()
    Hp = Polyhedron(ieqs=[(0,)*(k+1)+(1,)+(0,)*(n-k-1)])
    Mp = matrix(n, lambda i,j: (1 if i == j else 0) if j != k else (max(B[i,k],0) if i != k else -1) )
    Hm = Polyhedron(ieqs=[(0,)*(k+1)+(-1,)+(0,)*(n-k-1)])
    Mm = matrix(n, lambda i,j: (1 if i == j else 0) if j != k else (max(-B[i,k],0) if i != k else -1) )
    new_g = (Mp if g in Hp else Mm) * vector(g)
    new_B = copy(B)
    new_B.mutate(k)
    polytope = cut_along_sequence(new_g, new_B, seq)
    polytope_p = Mm*(polytope.intersection(Hp))
    polytope_m = Mp*(polytope.intersection(Hm))
    polytope = polytope_p.convex_hull(polytope_m)
    return polytope.intersection(current_cone)

def region_vertices(b,c,g):
    #this is only for rank 2
    vertices = []
    vertices.append(g)
    vertices.append(vector((0,(b*c+sqrt(b*c*(b*c-4)))*g[0]/(2*b)+g[1])))
    vertices.append(vector(((b/sqrt(b*c*(b*c-4)))*((b*c+sqrt(b*c*(b*c-4)))*g[0]/b+(b*c+sqrt(b*c*(b*c-4)))*g[1]/2),(c/sqrt(b*c*(b*c-4)))*((b*c+sqrt(b*c*(b*c-4)))*g[1]/c+(b*c+sqrt(b*c*(b*c-4)))*g[0]/2))))
    vertices.append(vector(((b*c+sqrt(b*c*(b*c-4)))*g[1]/(2*c)+g[0],0)))
    return vertices


def foo(g, B, seq):
    r"""
    Return the intersection of the cones dominated by all the translates of g along the sequence of mutations seq
    """
    current_cone = Polyhedron(rays=B.columns(),base_ring=QQ).translation(g)
    if seq == []:
        return current_cone
    k = seq.pop()
    n = B.ncols()
    Hp = Polyhedron(ieqs=[(0,)*(k+1)+(1,)+(0,)*(n-k-1)])
    Mp = matrix(n, lambda i,j: (1 if i == j else 0) if j != k else (max(B[i,k],0) if i != k else -1) )
    Hm = Polyhedron(ieqs=[(0,)*(k+1)+(-1,)+(0,)*(n-k-1)])
    Mm = matrix(n, lambda i,j: (1 if i == j else 0) if j != k else (max(-B[i,k],0) if i != k else -1) )
    new_g = (Mp if g in Hp else Mm) * vector(g)
    new_B = copy(B)
    new_B.mutate(k)
    polytope = foo(new_g, new_B, seq)
    polytope_p = Mm*(polytope.intersection(Hp))
    polytope_m = Mp*(polytope.intersection(Hm))
    polytope = polytope_p.convex_hull(polytope_m)
    return polytope #.intersection(current_cone)

def get_inequalities(B, seq, la):
    r"""
    Return the dominance region based at `la` corresponding to the sequence of mutations seq
    """
    n = B.ncols()
    if B.is_square():
        B = B.stack(identity_matrix(n))

    if len(la) == n:
        la = vector(tuple(la)+(0,)*n)
    else:
        la = vector(la)

    if seq == []:
        O = (0,)*n
        normals = B[n:].inverse().rows()
        return [ [ (vector(O+tuple(normal)), vector(normal)*la[n:]) for normal in normals ] ]

    k = seq.pop(0)

    def E(sgn):
        E = identity_matrix(2*n)
        E.set_column(k, list(map( lambda b: max(sgn*b,0), B.column(k))))
        E[k,k] = -1
        return E

    new_la = E(sign(la[k])) * la

    new_B = copy(B)
    new_B.mutate(k)

    new_regions = get_inequalities(new_B, seq, new_la)

    regions = []
    for new_region in new_regions:
        #print(new_region)
        for sgn in [-1, 1]:
            #P = Polyhedron(ieqs=[ (-c+norm*new_la,)+tuple(norm) for (norm,c) in new_region]+[ (-sgn*new_la[k],)+(0,)*(k) + (-sgn,) + (0,)*(2*n-k-1)], eqns=[ (new_la[m]-(new_B[:n]*new_B[n:].inverse())[m]*new_la[n:],) + (0,)*m + (-1,) + (0,)*(n-m-1) + tuple((new_B[:n]*new_B[n:].inverse())[m]) for m in range(n) ])
            #print(sgn,P.inequalities(),P.equations())
            #if P.is_empty():
            #    continue
            region = []
            for (new_normal, new_const) in new_region+[ (vector((0,)*(k) + (-sgn,) + (0,)*(2*n-k-1)), 0)]:
                normal = E(sgn).transpose()*new_normal
                const = new_const
                region.append( (normal, const) )
            regions.append(region)
    return regions


def get_inequalities_first_attempt(B, seq, la):
    r"""
    Return the dominance region based at `la` corresponding to the sequence of mutations seq
    """
    n = B.ncols()
    if B.is_square():
        B = B.stack(identity_matrix(n))

    if seq == []:
        O = (0,)*n
        normals = B[n:].inverse().rows()
        return [ [ (vector(O+tuple(normal)), 0) for normal in normals ] ]

    if len(la) == n:
        la = vector(tuple(la)+(0,)*n)
    else:
        la = vector(la)

    k = seq.pop(0)

    def E(sgn):
        E = identity_matrix(2*n)
        E.set_column(k, list(map( lambda b: max(sgn*b,0), B.column(k))))
        E[k,k] = -1
        return E

    new_la = E(sign(la[k])) * la

    new_B = copy(B)
    new_B.mutate(k)

    new_regions = get_inequalities_first_attempt(new_B, seq, new_la)

    regions = []
    for new_region in new_regions:
        print(new_region)
        for sgn in [-1, 1]:
            P = Polyhedron(ieqs=[ (-c+norm*new_la,)+tuple(norm) for (norm,c) in new_region]+[ (-sgn*new_la[k],)+(0,)*(k) + (-sgn,) + (0,)*(2*n-k-1)], eqns=[ (new_la[m]-(new_B[:n]*new_B[n:].inverse())[m]*new_la[n:],) + (0,)*m + (-1,) + (0,)*(n-m-1) + tuple((new_B[:n]*new_B[n:].inverse())[m]) for m in range(n) ])
            print(sgn,P.inequalities(),P.equations())
            if P.is_empty():
                continue
            region = []
            for (new_normal, new_const) in new_region+[ (vector((0,)*(k) + (-sgn,) + (0,)*(2*n-k-1)), sgn*new_la[k])]:
                normal = E(sgn).transpose()*new_normal
                const = new_const - normal*(la-E(sgn)*new_la)
                region.append( (normal, const) )
            regions.append(region)
    return regions


def mutation_map(B, k, la):
    B = B.transpose().stack(la)
    B.mutate(k)
    return vector(B.row(-1))
