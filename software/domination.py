from collections import defaultdict
from numpy import linspace
from sage.repl.rich_output import get_display_manager
dm = get_display_manager()
setattr(dm.preferences,'graphics',u'vector')

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
    return polytope #.intersection(current_cone)


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
    return polytope.intersection(current_cone)


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
    seq = copy(seq)
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


def get_regions(B, seq, la, return_steps=infinity,return_table=False):
    regions = defaultdict(list)
    for k, v in get_inequalities_with_return_steps(B, seq, la, return_steps):
        regions[Set(k)].append(v)
        #for kk in k:
        #    regions[kk].append(v)
    #print(len(regions))
    if return_table:
        return table(list(regions.items()))
    return regions


def get_combined_regions(B, seq, la, return_steps=infinity):
    if B.ncols() != 3:
        raise ValueError("Only works for rank 3.")

    regions = get_regions(B, seq, la, return_steps=return_steps)

    for key in regions:
        boundaries = []
        for cone in regions[key]:
            rays = list(map( lambda x: x.vector(), Polyhedron(ieqs=list(map(lambda x : [0]+list(x), cone))).rays()))
            vertices = list(zip(rays, rays[1:]+[rays[0]]))
            for v in vertices:
                v = set( map( tuple, v) )
                if not v in boundaries:
                    boundaries.append(v)
                else:
                    boundaries.remove(v)
        rays = set(flatten(list(map(list,boundaries)),max_level=1))
        #regions[key] = [ tuple(n(vector(x[1:]).normalized(),digits=4)) for x in Polyhedron(rays=rays).inequalities_list()]
        regions[key] = [ vector(x[1:]) for x in Polyhedron(rays=rays).inequalities_list()]

    return list(regions.values())
    return sorted(set(flatten(list(regions.values()),max_level=1)))


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


def stereo_regions(regions, north=(-1,-1,-1), right=(1,0,0), distance=-1, list_normals=False, check_single_region=True, random_colors=True, precision=50, fixed_color=""):
    north, right, up = _basis(north, right)
    C = matrix((right, up, north)).inverse().transpose()

    def proj(pt):
        # apparently the code is robust even if we try to project the point at infinity... it's a kind of magic?
        #if _normalize(pt) == north:
        pt = _normalize(C*pt)
        t = (distance-pt[2])/(1-pt[2])
        return (1-t)*pt[:2]

    G = Graphics()

    if random_colors:
        colors = list(sage.plot.colors.colors.values())
        shuffle(colors)
    else:
        colors = rainbow(len(regions))
    for region,color in zip(regions.values(),colors):
        boundaries = []
        for cone in region:
            rays = list(map( lambda x: x.vector(), Polyhedron(ieqs=list(map(lambda x : [0]+list(x), cone))).rays()))
            vertices = list(zip(rays, rays[1:]+[rays[0]]))
            for v in vertices:
                v = set( map( tuple, v) )
                if not v in boundaries:
                    boundaries.append(v)
                else:
                    boundaries.remove(v)
            triangle = flatten( [ [ proj(r2*s+r1*(1-s)) for s in linspace(0,1,precision*ceil(norm(proj(r1)-proj(r2)))) ] for (r1,r2) in vertices], max_level=1)

            if not fixed_color:
                if all( sum( a*b for (a,b) in zip(ieq, north)) >=0 for ieq in cone ):
                    # HERE adjust the radius
                    G += circle((0,0), max(norm(proj(r))**(3/2) for r in rays), color=color, zorder=0, fill=True)
                    G += polygon( triangle, color='white', zorder=0)
                else:
                    G += polygon( triangle, color=color)

        for r1,r2 in boundaries:
            thickness = 0.3
            if fixed_color:
                thickness = 1
            G += line([ proj(vector(r1)*s+vector(r2)*(1-s)) for s in linspace(0,1,precision*ceil(norm(proj(vector(r1))-proj(vector(r2))))) ], thickness=thickness, color="black")

        #if check_single_region:
        #    bounderies_bkp = copy(boundaries)
        #    side = boundaries.pop()
        #    last_vertex = side.pop()
        #    current_vertex = side.pop()
        #    while current_vertex != last_vertex:
        #        sides = [ s for s in boundaries if current_vertex in s ]
        #        side = sides[0]
        #        if len(sides)>1:
        #            break
        #        boundaries.remove(side)
        #        current_vertex = [ v for v in side if v != current_vertex ][0]
        #    if boundaries:
        #        print("Warning: the following is a problematic region")
        #        print("region: ",region)
        #        print("boundaries: ",boundaries_bkp)

    G.set_aspect_ratio(1)
    G.SHOW_OPTIONS['axes']=False
    return G


def stereo_fan(fan, north=-vector((1,1,1)), right=vector((1,0,0)), distance=-1, list_normals=False, color_initial=True, color_final=True, precision=50, plot_vectors=[]):
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
                G += line([ proj(r1*s+r2*(1-s)) for s in linspace(0,1,precision*ceil(norm(proj(r1)-proj(r2)))) ], rgbcolor=tuple(_normalize(normal)))
            else:
                seen_normals.append(normal)
                G += line([ proj(r1*s+r2*(1-s)) for s in linspace(0,1,precision*ceil(norm(proj(r1)-proj(r2)))) ], rgbcolor=tuple(_normalize(normal)), legend_label=normal)
    else:
        for w in fan.cones(2):
            r1, r2 = map(_normalize,w.rays())
            G += line([ proj(r1*s+r2*(1-s)) for s in linspace(0,1,precision*ceil(norm(proj(r1)-proj(r2)))) ], color='black', thickness=0.001)

    # plot rays
    for r in fan.rays():
        G += point(proj(r),color='black', zorder=len(G)+1, size=0.5)
        #G += text(tuple(r), proj(r), fontsize='x-small',horizontal_alignment='right')

    # plot initial cluster
    if color_initial:
        units = identity_matrix(3).columns()
        for r, color in zip(units,['red','green','blue']):
            G += point(proj(r),color=color, zorder=len(G)+1, size=20)
        #vertices = zip(units, units[1:]+[units[0]])
        #triangle = flatten( [ [ proj(r2*s+r1*(1-s)) for s in linspace(0,1,precision*ceil(norm(proj(r1)-proj(r2)))) ] for (r1,r2) in vertices], max_level=1)
        #if all(i>0 for i in list(north)):
        #    G += circle((0,0), max(sqrt(G.xmin()**2+G.ymin()**2), sqrt(G.xmax()**2+G.ymax()**2)), color='lightgreen', zorder=0, fill=True)
        #    G += polygon( triangle, color='white', zorder=0)
        #else:
        #    G += polygon( triangle, color='white')

    # plot final cluster
    if color_final:
        units = (-identity_matrix(3)).columns()
        for r, color in zip(units,['red','green','blue']):
            G += point(proj(r),color=color, zorder=len(G)+1, size=20)
        #vertices = zip(units, units[1:]+[units[0]])
        #triangle = flatten( [ [ proj(r2*s+r1*(1-s)) for s in linspace(0,1,precision*ceil(norm(proj(r1)-proj(r2)))) ] for (r1,r2) in vertices], max_level=1)
        #if all(i<0 for i in list(north)):
        #    G += circle((0,0), max(sqrt(G.xmin()**2+G.ymin()**2), sqrt(G.xmax()**2+G.ymax()**2)), color='red', zorder=0, fill=True)
        #    G += polygon( triangle, color='white', zorder=0)
        #else:
        #    G += polygon( triangle, color='red')

    for plot_vector in plot_vectors:
        G += point(proj(vector(plot_vector)),color='black', zorder=len(G)+1, size=20, legend_label=str(plot_vector))
        for eps, color in zip(list(cartesian_product([[1,-1]]*B.ncols())),rainbow(8)):
            G += point(proj(tau(B,eps)*vector(plot_vector)),color='black', marker='d', markeredgecolor=color, zorder=len(G)+1, size=20, legend_label="ep="+str(eps))

    G.set_aspect_ratio(1)
    G.SHOW_OPTIONS['axes']=False
    return G


def stereo(B, seq, la, return_steps=+Infinity, depth=6, north=(-1,-1,-1), distance=-1, check_single_region=True, random_colors=False, precision=50, plot_vectors=[]):
    tau_linearity_regions = get_regions(B,source_sequence(B),(0,)*B.ncols())

    if return_steps >= len(seq):
        mutated_B = B
    else:
        mutated_B = copy(B)
        for k in seq[:-return_steps]:
            mutated_B.mutate(k)
    A = ClusterAlgebra(mutated_B)

    #plot limiting hyperplane
    ev = tau(B,(-1,-1,-1),inv=True).eigenvectors_left()[-1][1][0]
    right=vector((1,0,0))
    north, right, up = _basis(north, right)
    C = matrix((right, up, north)).inverse().transpose()

    def proj(pt):
        # apparently the code is robust even if we try to project the point at infinity... it's a kind of magic?
        #if _normalize(pt) == north:
        pt = _normalize(C*pt)
        t = (distance-pt[2])/(1-pt[2])
        return (1-t)*pt[:2]

    vecs = [vector((ev[1],-ev[0],0)),vector((ev[2],0,-ev[0])),vector((0,ev[2],-ev[1]))]
    vertices = list(zip(vecs, vecs[1:]+[vecs[0]]))
    triangle = flatten( [ [ proj(r2*s+r1*(1-s)) for s in linspace(0,1,precision*ceil(norm(proj(r1)-proj(r2)))) ] for (r1,r2) in vertices], max_level=1)
    G = polygon( triangle, fill=False, color="green")

    #seed = A.initial_seed()
    #seed.mutate(seq[-return_steps:],mutating_F=False)
    #north = sum(vector(g) for g in seed.g_vectors())
    F = A.cluster_fan(depth=depth)
    regions = get_regions(B, seq, la, return_steps)
    return stereo_regions(regions, distance=distance, check_single_region=check_single_region, random_colors=random_colors,precision=precision, north=north) + stereo_regions(tau_linearity_regions, distance=distance, check_single_region=check_single_region, random_colors=random_colors, precision=precision, north=north, fixed_color="blue") + stereo_fan(F,color_initial=True, color_final=True, distance=distance,precision=precision, north=north, plot_vectors=plot_vectors)+G


def regions_with_equalities(B, seq, la, return_table=False, projection='c'):
    if projection == 'g':
        if B.rank() != B.ncols():
            raise ValueError("B is not full rank")
    elif projection != 'c' and projection:
        raise ValueError("What projection are you trying to do??")
    n = B.ncols()
    la = vector(la)
    regions = get_regions(B, seq, la)
    new_regions = {}
    for (key,value) in regions.items():
        if projection == 'c':
            new_key = Set( (tuple(vector(lhs[:n])*B+vector(lhs[n:])),rhs-vector(lhs[:n]).dot_product(la)) for (lhs,rhs) in key)
            # lhs[:n][-I B]mu = -lhs[:n]*lambda
        elif projection == 'g':
            new_key = Set( (tuple(vector(lhs[:n])+vector(lhs[n:])*B.inverse()),rhs+vector(lhs[n:])*B.inverse()*la) for (lhs,rhs) in key)
            # lhs[n:][B.inverse() -I]mu = lhs[n:]*B.inverse()*lambda
        else:
            new_key = key
        new_value = []
        for cone in value:
            if projection == 'c':
                new_value.append([(vector(ieq)*B,-vector(ieq).dot_product(la)) for ieq in cone])
            elif projection == 'g':
                new_value.append([(ieq,0) for ieq in cone])
            else:
                new_value.append([(tuple(ieq)+(0,)*n,0) for ieq in cone])
        new_regions[new_key] = new_value
    if return_table:
        return table(list(new_regions.items()))
    return new_regions


def get_polyhedra(B, seq, la, projection='c'):
    regions =  regions_with_equalities(B, seq, la, projection=projection)
    polyhedra = flatten( [ [ Polyhedron(ieqs=[ [-rhs]+list(lhs) for (lhs,rhs) in key ] +[[-rhs] + list(lhs) for (lhs,rhs) in cone]) for cone in value ] for (key,value) in regions.items()], max_level=1 )
    return polyhedra


def get_polyhedron(B, seq, la, projection='c'):
    polyhedra = get_polyhedra(B, seq, la, projection=projection)
    dim = B.ncols()
    if not projection:
        dim = 2*dim
    Q = Polyhedron(ambient_dim=dim,base_ring=QQ)
    for P in polyhedra:
        Q = Q.convex_hull(P)
    return Q


def get_intersected_polyhedron(B, seq, la, projection='c'):
    P = get_polyhedron(B, seq, la, projection=projection)
    for k in range(len(seq)):
        old_P = copy(P)
        P = P.intersection(get_polyhedron(B, seq[:k], la, projection=projection))
    return P


def mutate_polyhedron(P, B, seq):
    r"""
    Actually broken since we're in the wrong lattice
    """
    polytope = P
    if seq == []:
        return polytope
    k = seq.pop(0)
    n = B.ncols()
    Hp = Polyhedron(ieqs=[(0,)*(k+1)+(1,)+(0,)*(n-k-1)])
    Mp = matrix(n, lambda i,j: (1 if i == j else 0) if j != k else (max(B[i,k],0) if i != k else -1) )
    Hm = Polyhedron(ieqs=[(0,)*(k+1)+(-1,)+(0,)*(n-k-1)])
    Mm = matrix(n, lambda i,j: (1 if i == j else 0) if j != k else (max(-B[i,k],0) if i != k else -1) )
    polytope_p = Mm*(polytope.intersection(Hp))
    polytope_m = Mp*(polytope.intersection(Hm))
    polytope = polytope_p.convex_hull(polytope_m)
    new_B = copy(B)
    new_B.mutate(k)
    polytope = mutate_polyhedron(polytope, new_B, seq)
    return polytope


def filtered_inequalities(ieqs, tolerance=0.0000001):
    ieqs = [ n(ieq.vector()/ieq.vector()[1:].norm()) for ieq in ieqs ]
    filtered_ieqs = {}
    for ieq in ieqs:
        found = False
        for normal in filtered_ieqs:
            if sum(abs(ieq[i+1]-normal[i]) for i in range(B.ncols())) < tolerance:
                filtered_ieqs[normal].append(ieq[0])
                found = True
                break
        if not found:
            filtered_ieqs[tuple(ieq[1:])] = [ieq[0]]
    return [ (min(filtered_ieqs[key]),)+tuple(key) for key in sorted(filtered_ieqs.keys()) ]


def find_inequalities(B, seq_type, num_steps, la, projection='c', tolerance=0.00000000000001):
    n = B.ncols()
    if seq_type == 'source':
        inv_seq = False
        taus_B = B.transpose()
    elif seq_type == 'sink':
        inv_seq = True
        taus_B = B
    else:
        raise ValueError("Unknown sequence type requested.")
    seq = source_sequence(B,inv=inv_seq)
    poly = get_intersected_polyhedron(B, seq*floor(num_steps/n)+seq[:num_steps%n], la, projection=projection)
    not_found_ieqs = filtered_inequalities(poly.inequalities(), tolerance=tolerance)

    ev = tau(taus_B,(-1,)*B.ncols(),inv=True).eigenvectors_right()[-1][1][0].normalized()
    if seq_type == 'source':
        ev = -ev

    working_ieqs = {}
    for ieq in not_found_ieqs:
        if sum(abs(ieq[i+1]-ev[i]) for i in range(B.ncols())) < tolerance:
            working_ieqs[ieq] = []
    for ieq in working_ieqs:
        not_found_ieqs.remove(ieq)

    found_ieqs = copy(working_ieqs)
    while working_ieqs:
        key = list(working_ieqs.keys())[0]
        eps_list = working_ieqs.pop(key)
        for eps in cartesian_product([[1,-1]]*n):
            checking = True
            if eps_list != []:
                if eps_list[-1] == eps:
                    checking = False
            if checking:
                new_normal = (tau(taus_B, eps, inv=True)*vector(key[1:])).normalized()
                found = []
                for ieq in not_found_ieqs:
                    if sum(abs(ieq[i+1]-new_normal[i]) for i in range(B.ncols())) < tolerance:
                        working_ieqs[ieq] = eps_list+[eps]
                        found_ieqs[ieq] = eps_list+[eps]
                        found.append(ieq)
                for ieq in found:
                    not_found_ieqs.remove(ieq)

    print("Not found: ",not_found_ieqs)
    print("Found: ",found_ieqs)


#Exchange matrix methods

def source_sequence(B,inv=False):
    seq = []
    Bp = copy(B)
    n = B.ncols()
    col_list = list(range(n))
    while col_list != []:
        for k in col_list:
            if all( Bp[l][k] >= 0 for l in range(n) ):
                seq.append(k)
                col_list.remove(k)
                Bp.mutate(k)
                break
    if inv:
        seq.reverse()
    return seq


def E(B,k,sgn):
     n = B.ncols()
     E = identity_matrix(n)
     E.set_column(k, list(map( lambda b: max(sgn*b,0), B.column(k))))
     E[k,k] = -1
     return E


def tau(B, eps, inv=False):
     n = B.ncols()
     Bp = copy(B)
     tau = identity_matrix(n)
     ss = source_sequence(B, inv)
     for k in range(n):
         tau *= E(Bp,ss[k],eps[k])
         Bp.mutate(ss[k])
     return tau

def chebyshev(k,eps,b,c):
    if k == 0:
        return 0
    if k == 1:
        return 1
    if eps == 1:
        return b*chebyshev(k-1,-eps,b,c) - chebyshev(k-2,eps,b,c)
    else:
        return c*chebyshev(k-1,-eps,b,c) - chebyshev(k-2,eps,b,c)


def phi(B, seq):
    B = copy(B)
    I = identity_matrix(B.nrows())
    
    if seq == []:
        I.set_immutable()
        return ({I:[[]]},B)

    k = seq.pop()
    current_regions, current_B = phi(B, seq)
    
    Ep = matrix(B.nrows(), lambda i,j: 0 if j!=k else current_B[i,j] if current_B[i,j] > 0 else 0) + 1
    Ep[k,k] = -1

    Em = matrix(B.nrows(), lambda i,j: 0 if j!=k else -current_B[i,j] if current_B[i,j] < 0 else 0) + 1
    Em[k,k] = -1
    
    new_regions = {}
    for transformation in current_regions:
        new_transformation = Ep*transformation
        new_transformation.set_immutable()
        for ieqs in current_regions[transformation]:
            P = Polyhedron(ieqs=ieqs+[[0]+list(transformation[k])])
            if P.dimension() == B.nrows():
                new_regions[new_transformation] = new_regions.get(new_transformation,[])+[P.inequalities_list()]
        new_transformation = Em*transformation
        new_transformation.set_immutable()
        for ieqs in current_regions[transformation]:
            P = Polyhedron(ieqs=ieqs+[[0]+list(-transformation[k])])
            if P.dimension() == B.nrows():
                new_regions[new_transformation] = new_regions.get(new_transformation,[])+[P.inequalities_list()]

    current_B.mutate(k)

    return (new_regions, current_B)
