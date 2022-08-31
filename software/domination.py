from collections import defaultdict
from numpy import linspace

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

def get_inequalities_broken(B, seq, la):
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

    new_regions = get_inequalities_broken(new_B, seq, new_la)

    regions = []
    for new_region in new_regions:
        #print(new_region)
        for sgn in [-1, 1]:
            #P = Polyhedron(ieqs=[ (-c,)+tuple(norm) for (norm,c) in new_region]+[ (-sgn*new_la[k],)+(0,)*(k) + (-sgn,) + (0,)*(2*n-k-1)], eqns=[ (new_la[m]-(new_B[:n]*new_B[n:].inverse())[m]*new_la[n:],) + (0,)*m + (-1,) + (0,)*(n-m-1) + tuple((new_B[:n]*new_B[n:].inverse())[m]) for m in range(n) ])
            #print(sgn,P.inequalities(),P.equations())
            #if P.is_empty():
            #    print("It is empty")
            #   continue
            region = []
            for (new_normal, new_const) in new_region+[ (vector((0,)*(k) + (-sgn,) + (0,)*(2*n-k-1)), 0)]:
                normal = E(sgn).transpose()*new_normal
                const = new_const
                region.append( (normal, const) )
            regions.append(region)
    return regions

def get_region_broken(B, seq, la):
    n = B.ncols()
    if B.is_square():
        B = B.stack(identity_matrix(n))
    if len(la) == n:
        la = vector(tuple(la)+(0,)*n)
    else:
        la = vector(la)

    eqns = [ (la[m]-(B[:n]*B[n:].inverse())[m]*la[n:],) + (0,)*m + (-1,) + (0,)*(n-m-1) + tuple((B[:n]*B[n:].inverse())[m]) for m in range(n) ]
    polys = [ Polyhedron(ieqs=[ (-c,)+tuple(v) for (v,c) in ieqs],eqns=eqns) for ieqs in get_inequalities_broken(B,seq,la) ]
    P = polys[0]
    for Q in polys:
        P = P.convex_hull(Q)
    return P

def get_inequalities_with_return_steps(B, seq, la, return_steps=infinity):
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

    B = block_matrix([[B,la.column()]])
    while len(seq) > return_steps:
        k = seq.pop(0)
        B.mutate(k)

    la = vector(B[:,-1])
    B = B[:,:-1]

    return get_inequalities(B, seq, la)

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
        return [ ( tuple( (O+tuple(normal), vector(normal)*la[n:]) for normal in normals ), [vector(O)] ) ]

    k = seq.pop(0)

    def E(sgn):
        E = identity_matrix(2*n)
        E.set_column(k, list(map( lambda b: max(sgn*b,0), B.column(k))))
        E[k,k] = -1
        return E

    la_p = E(sign(la[k])) * la

    B_p = copy(B)
    B_p.mutate(k)

    inequalities_p = get_inequalities(B_p, seq, la_p)

    inequalities = []
    for (ieqs_p, region_p) in inequalities_p:
        for sgn in [-1, 1]:
            ieqs = tuple( (tuple(E(sgn).transpose()*vector(normal_p)), const) for (normal_p, const) in ieqs_p )
            region = [ E(sgn).transpose()[:n,:n]*normal_p for normal_p in region_p + [ vector((0,)*(k) + (-sgn,) + (0,)*(n-k-1))] ]
            C = Polyhedron(ieqs = [ (0,)+tuple(v) for v in region ])
            if C.dim() != n:
                continue
            inequalities.append( (ieqs, [ vector(v[1:]) for v in C.inequalities_list()]) )
    return inequalities

def get_regions(B, seq, la, return_steps=infinity):
    regions = defaultdict(list)
    for k, v in get_inequalities_with_return_steps(B, seq, la, return_steps):
        regions[Set(k)].append(v)
        #for kk in k:
        #    regions[kk].append(v)
    return regions
    print(len(regions))
    return table(list(regions.items()))


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



##### Try to understand pairings
def explore_function(B, length=5):
    A = ClusterAlgebra(B)
    AA = ClusterAlgebra(-B.transpose())
    n = B.ncols()
    good = [ 0, 0, 0, 0 ]
    negative = 0
    positive = 0
    not_coherent = 0
    bad = 0
    weird_prods = defaultdict(list)
    for s,t in cartesian_product([range(length+1)]*2):
        for i in cartesian_product([range(n)]*s):
            if all( i[k] != i[k+1] for k in range(len(i)-1)):
                S = A.initial_seed()
                S.mutate(i, mutating_F=False)
                for j in cartesian_product([range(n)]*t):
                    if all( j[k] != j[k+1] for k in range(len(j)-1)):
                        if i == j:
                            continue
                        SS = AA.initial_seed()
                        SS.mutate(j, mutating_F=False)
                        M = SS.g_matrix().transpose()*S.c_matrix()
                        if not all( all( x >= 0 for x in v ) or all( x <= 0 for x in v ) for v in M.columns() ):
                            not_coherent += 1
                            weird_prods['nsc'].append([i,j,M])
                            continue
                        As = []
                        As.append(ClusterAlgebra(-SS.b_matrix().transpose()))
                        As.append(ClusterAlgebra(SS.b_matrix().transpose()))
                        As.append(ClusterAlgebra(S.b_matrix()))
                        As.append(ClusterAlgebra(-S.b_matrix()))
                        Ss = [ _.initial_seed() for _ in As ]
                        is_good = False
                        k = 0
                        for s in Ss:
                            s.mutate(list(reversed(list(j)))+list(i), mutating_F=false)
                            if s.c_matrix() == M:
                                good[k] += 1
                                is_good = True
                                break
                            k += 1
                        if not is_good:
                            if all( x <=0 for x in M.list() ):
                                negative += 1
                                weird_prods['n'].append([i,j,M])
                            elif all( x >=0 for x in M.list() ):
                                positive += 1
                            else:
                                bad += 1
                                weird_prods['b'].append([i,j,M])
    print("g:", good, "p:", positive, "n:", negative, "nsc:", not_coherent, "b:", bad)
    return weird_prods


def _basis(north, right):
    north = _normalize(north)
    right = _normalize(right)
    right = _normalize(right - north.dot_product(right)*north)
    up = north.cross_product(right)
    return (north, right, up)

def _normalize(v):
    return vector(v,RR)/N(vector(v).norm())

def stereo_regions(regions, north=(-1,-1,-1), right=(1,0,0), distance=-1, list_normals=False):
    north, right, up = _basis(north, right)
    C = matrix((right, up, north)).inverse().transpose()

    def proj(pt):
        # apparently the code is robust even if we try to project the point at infinity... it's a kind of magic?
        #if _normalize(pt) == north:
        pt = _normalize(C*pt)
        t = (distance-pt[2])/(1-pt[2])
        return (1-t)*pt[:2]

    G = Graphics()

    colors = rainbow(len(regions))
    for region,color in zip(regions.values(),colors):
        for cone in region:
            rays = list(map( lambda x: x.vector(), Polyhedron(ieqs=list(map(lambda x : [0]+list(x), cone))).rays()))
            vertices = zip(rays, rays[1:]+[rays[0]])
            triangle = flatten( [ [ proj(r2*s+r1*(1-s)) for s in linspace(0,1,50*ceil(norm(proj(r1)-proj(r2)))) ] for (r1,r2) in vertices], max_level=1)
            if all( sum( a*b for (a,b) in zip(ieq, north)) >=0 for ieq in cone ):
                # HERE adjust the radius
                G += circle((0,0), max(norm(proj(r))**(3/2) for r in rays), color=color, zorder=0, fill=True)
                G += polygon( triangle, color='white', zorder=0)
            else:
                G += polygon( triangle, color=color)

    G.set_aspect_ratio(1)
    G.SHOW_OPTIONS['axes']=False
    return G

def stereo_fan(fan, north=-vector((1,1,1)), right=vector((1,0,0)), distance=-1, list_normals=False, color_initial=True, color_final=True):
    if fan.dim() != 3:
        raise ValueError("Can only stereographically project 3-dimensional fans.")

    north, right, up = _basis(north, right)
    C = matrix((right, up, north)).inverse().transpose()

    def proj(pt):
        # apparently the code is robust even if we try to project the point at infinity... it's a kind of magic?
        #if _normalize(pt) == north:
        pt = _normalize(C*pt)
        t = (distance-pt[2])/(1-pt[2])
        return (1-t)*pt[:2]

    G = Graphics()

    # plot walls
    if list_normals:
        seen_normals = []
        G.set_legend_options(title='c-vectors:')
        for w in fan.cones(2):
            r1, r2 = map(_normalize,w.rays())
            normal = tuple(w.orthogonal_sublattice().gens()[0])
            if normal in seen_normals:
                G += line([ proj(r1*s+r2*(1-s)) for s in linspace(0,1,30*ceil(norm(proj(r1)-proj(r2)))) ], rgbcolor=tuple(_normalize(normal)))
            else:
                seen_normals.append(normal)
                G += line([ proj(r1*s+r2*(1-s)) for s in linspace(0,1,30*ceil(norm(proj(r1)-proj(r2)))) ], rgbcolor=tuple(_normalize(normal)), legend_label=normal)
    else:
        for w in fan.cones(2):
            r1, r2 = map(_normalize,w.rays())
            G += line([ proj(r1*s+r2*(1-s)) for s in linspace(0,1,30*ceil(norm(proj(r1)-proj(r2)))) ], color='black')

    # plot rays
    for r in fan.rays():
        G += point(proj(r),color='black', zorder=len(G)+1)
        #G += text(tuple(r), proj(r), fontsize='x-small',horizontal_alignment='right')

    # plot initial cluster
    if color_initial:
        units = identity_matrix(3).columns()
        for r, color in zip(units,['red','green','blue']):
            G += point(proj(r),color=color, zorder=len(G)+1, size=20)
        vertices = zip(units, units[1:]+[units[0]])
        triangle = flatten( [ [ proj(r2*s+r1*(1-s)) for s in linspace(0,1,30*ceil(norm(proj(r1)-proj(r2)))) ] for (r1,r2) in vertices], max_level=1)
        if all(i>0 for i in list(north)):
            G += circle((0,0), max(sqrt(G.xmin()**2+G.ymin()**2), sqrt(G.xmax()**2+G.ymax()**2)), color='lightgreen', zorder=0, fill=True)
            G += polygon( triangle, color='white', zorder=0)
        else:
            G += polygon( triangle, color='lightgreen')

    # plot final cluster
    if color_final:
        units = (-identity_matrix(3)).columns()
        for r, color in zip(units,['red','green','blue']):
            G += point(proj(r),color=color, zorder=len(G)+1, size=20)
        vertices = zip(units, units[1:]+[units[0]])
        triangle = flatten( [ [ proj(r2*s+r1*(1-s)) for s in linspace(0,1,30*ceil(norm(proj(r1)-proj(r2)))) ] for (r1,r2) in vertices], max_level=1)
        if all(i<0 for i in list(north)):
            G += circle((0,0), max(sqrt(G.xmin()**2+G.ymin()**2), sqrt(G.xmax()**2+G.ymax()**2)), color='coral', zorder=0, fill=True)
            G += polygon( triangle, color='white', zorder=0)
        else:
            G += polygon( triangle, color='coral')


    G.set_aspect_ratio(1)
    G.SHOW_OPTIONS['axes']=False
    return G

def stereo(B, seq, la, return_steps=+Infinity, depth=6):
    if return_steps >= len(seq):
        mutated_B = B
    else:
        mutated_B = copy(B)
        for k in seq[:-return_steps]:
            mutated_B.mutate(k)
    A = ClusterAlgebra(mutated_B)
    F = A.cluster_fan(depth=depth)
    regions = get_regions(B, seq, la, return_steps)
    return stereo_regions(regions) + stereo_fan(F,color_initial=false, color_final=false)
