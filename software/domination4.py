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

# Added by Nathan:

def p_lambda_int(B, la, seqs):
    P = p_lambda_faster(B, la, seqs[0])
    for s in seqs:
        P = P.intersection(p_lambda_faster(B, la, s))
        if P.dimension()==0:
            break
    return P


def p_lambda_faster(B, la, seq):  # faster?  I think so
    m = B.nrows()
    rk = B.rank()
    n = B.ncols()
    B = block_matrix([[B,matrix(la).transpose()]])
    for k in seq:
        B.mutate(k)
    B, la  = B[:,:-1], B[:,-1] 
    P = Polyhedron(rays=B.columns(),base_ring=QQ).translation(la)
    #print("P:   ",P,P.vertices(),P.rays(),P.lines(),"\n")
    for k in reversed(seq):
        Ep = E(B, k, 1)
        Em = E(B, k, -1)
        Hp = Polyhedron(ieqs=[(0,)*(k+1)+(1,)+(0,)*(m-k-1)])
        Hm = Polyhedron(ieqs=[(0,)*(k+1)+(-1,)+(0,)*(m-k-1)])
        #print("Hp:   ",Hp,Hp.vertices(),Hp.rays(),Hp.lines(),"\n")
        #print("Hm:   ",Hm,Hm.vertices(),Hm.rays(),Hm.lines(),"\n\n")
        Pp = P.intersection(Hp)
        Pm = P.intersection(Hm)
        if Pp.dimension() < rk:
            P=Em*Pm
        elif Pm.dimension() < rk:
            P=Ep*Pp
        else:
            #print("Pp:   ",Pp,Pp.vertices(),Pp.rays(),Pp.lines(),"\n")
            #print("Pm:   ",Pm,Pm.vertices(),Pm.rays(),Pm.lines(),"\n\n")
            P=(Ep*Pp).convex_hull(Em*Pm)
        B.mutate(k)
    return P


# Given a polyhedron P and a vertex v, finds the smallest cone at v that contains P.
# The point is that you need to keep the lines and the rays and add a ray in the 
# direction v'-v for each other vertex v'.
# This will give strange results unless v is a vertex of P or is on a line or ray of P.
# To make it better, you would need to do a lot of checking about whether P is in the relative 
# interior of a face and I don't want to waste that computation time.  We'll only use this when 
# v is the lambda that defines P as a P_lambda,k
# The vertex should be given as a tuple, e.g. (1,2,3)
def cone_at(P,v):
    newrays=list(r.vector() for r in P.rays());
    for vert in P.vertices():
        vertvec=vert.vector()
        if vertvec!=v:
            newrays=newrays+[vertvec-v]
    return Polyhedron(vertices=[v],rays=newrays,lines=P.lines())

# like p_lambda_int, but returns the (often larger) cone at lambda that contains p_lambda_int
# I guessed that this would be faster in situations where we think we are getting a point, but
# it seems about the same
def p_lambda_cone(B, la, seqs):
    P = cone_at(p_lambda_faster(B, la, seqs[0]),la)
    for s in seqs:
        P = P.intersection(cone_at(p_lambda_faster(B, la, s),la))
        if P.dimension()==0:
            break
    return P


def B(A,c):  # Cartan matrix (assumes nonpositive off-diagonal entries) and Coxeter element (a list)
    n=A.nrows()
    out=Matrix([[0]*n]*n)
    #print(out)
    for i in range(n):
        out[i,i]=0
        for j in range(i+1,n):
            out[c[i],c[j]]=-A[c[i],c[j]]
            out[c[j],c[i]]=A[c[j],c[i]]
    return out

def K(c,v):  # A set of sequences depending on a sortable element v, given by its c-sorting word
    out = []
    if c==[]:   # The terminal case (empty sequence for the rank-0 Coxeter group)
        out=[[]]
    elif v==[] or c[0]!=v[0]:  # The "induction on rank" case.
        oldK=K(c[1:],v)
        for k in oldK:
            out=out+[k+[c[0]]]
        out=oldK+out
    else:  # The "induction on length" case.  This is provably the right thing to do.
        for k in K(c[1:]+[c[0]],v[1:]):
            if k!=[] and k[0]==v[0]:
                out=out+[k[1:]]
            else:
                out=out+[[v[0]]+k]
    return out

# a set of sequences, all prefixes of c^infty.  The induction-on-rank step is inserting c[0] into
# the sequences in the parabolic in the right places to make them prefixes of c^infty

def K_insert(c,v):  
    if c==[]:   # The terminal case (empty sequence for the rank-0 Coxeter group)
        out=[[]]
    elif v==[] or c[0]!=v[0]:  # The "induction on rank" case.
        out=[]
        oldK=K_insert(c[1:],v)
        for k in oldK:
            newk = insert(c,k)
            if newk==[] or newk[-1]==c[-1]:
                out=out+[newk,newk+[c[0]]]
            else:
                out=out+[newk]
    else:  # The "induction on length" case.  This is provably the right thing to do.
        out=[]
        for k in K_insert(c[1:]+[c[0]],v[1:]):
            if k!=[] and k[0]==v[0]:
                out=out+[k[1:]]
            else:
                out=out+[[v[0]]+k]
    return out

def insert(c,k):
    if k==[]:
        return []
    elif k[0]==c[1]:
        return [c[0],k[0]]+insert(c,k[1:])
    else:
        return [k[0]]+insert(c,k[1:])

def Kc(c,coxnum):  # powers of c
    return list(c*i for i in range(coxnum+2)) 

def Kc_even(c,coxnum):  # even powers of c
    return list(2*c*i for i in range(coxnum+1)) 

def Kprefix(c,coxnum):  # prefixes of powers of c
    return list(c_inf_prefix(c,i) for i in range(len(c)*coxnum)) 

def c_inf_prefix(c,i):
    return list(c[j%len(c)] for j in range(i))

def Krev_prefix(c,coxnum):  # prefixes of powers of c
    return list(c_inf_prefix(c,i) for i in range(len(c)*(coxnum+2))) 

def c_inf_rev_prefix(c,i):
    return list(c[(-j)%len(c)] for j in range(i))

# Trying to find an index I such that the intersection of the P_{lambda,k} for sequences
# k=c_inf_prefix(c,I),...,c_inf_prefix(c,I+n) is a point.  (Takes a sortable element v, and 
# we would take lambda to be in the seed for v.
def prefix_start_index(c,v):
    if v==[]:
        return 0
    elif v[0]==c[0]:
        return prefix_start_index(c[1:]+[c[0]],v[1:])+1
    else:
        p=prefix_start_index(c[1:],v)
        return p//(len(c)-1)*len(c)+p%(len(c)-1)



def Kbip(cplus,cminus,Coxnum):
    out = [[]]
    for i in range(1,Coxnum):
        out=out+[plusminus(cplus,cminus,-i/float(2)),plusminus(cplus,cminus,i/float(2))]
    out=out+[plusminus(cplus,cminus,Coxnum/float(2))]
    return out


# c=cplus+cminus
# cinv=cminus+cplus
# returns [] if n is 0
# returns c^n if n is a positive whole number
# returns c^(n-0.5)+cplus if n is a positive number with fractional part 0.5
# returns cinv^(-n) if n is a negative whole number
# returns cinv^(-n-0.5)+cminus if n is a positive number with fractional part 0.5
def plusminus(cplus,cminus,n):
    if abs(n)<0.1:  # trying to be careful of the "is a float zero" issue.
        return []
    elif n>0:
        return cplus+plusminus(cminus,cplus,n-0.5)
    else:
        return cminus+plusminus(cminus,cplus,n+0.5)



#The reflection $s_i$ acts on $V^*$ by fixing $\rho_j$ for $j\neq i$ and sending $\rho_i$ to $\rho_i-\sum_{k=1}^na_{ki}\rho_k$.


# reflects a vector by a simple reflection
def reflect(A,vec,s):
    return vector(vec)-vec[s]*vector(A.column(s))

# a vector in the interior of the Coxeter cone given by the sortable element
# (and therefore in the cone of the g-vector fan given by the sortable element)
def lam(A,v):
    n=A.nrows()
    out=vector([1]*n)
    for s in reversed(v):
        out=reflect(A,out,s)
    return out

# finds the intersection of the P_lambda,k for lambda=lam(A,v) and all k in K(c,v)
def p_lambda_v(A,c,v):
    return p_lambda_int(B(A,c),lam(A,v),K(c,v))

def sortables(A,c):
    yield from sortables_remaining(A,c,[],vector([1]*A.nrows()))

def sortables_remaining(A,c,sort,vec):
    if c==[]:
        yield sort
    else:
        if vec[c[0]]>0:
            yield from sortables_remaining(A,c[1:]+[c[0]],sort+[c[0]],reflect(A,vec,c[0]))
        yield from sortables_remaining(A,c[1:],sort,vec)

def parasortables(A,c):  # c-sortable elements not starting with c, i.e. in some proper parabolic
    for v in sortables(A,c):
        if len(v)<len(c) or v[0:len(c)]!=c:
            yield v

def longest_sortable(A,c):
    return longest_sortable_remaining(A,c,[],vector([1]*A.nrows()))

def longest_sortable_remaining(A,c,sort,vec):
    if c==[]:
        return sort
    else:
        if vec[c[0]]>0:
            return longest_sortable_remaining(A,c[1:]+[c[0]],sort+[c[0]],reflect(A,vec,c[0]))
        else:
            return longest_sortable_remaining(A,c[1:],sort,vec)

# Is this the dumbest way to find the Coxeter number or the smartest?
def coxeter_number(A):
    n=A.nrows()
    return len(longest_sortable(A,list(range(n))))*2/n

# Each c-sortable element v encodes a cone in the c-Cambrian fan

# This finds the c^{-1}-sortable element that encodes the same cone in the antipodal c^{-1}-Cambrian fan.
def inv_sortable(A,c,v):
    return inv_sortable_remaining(A,c,[],lam(A,v))

def inv_sortable_remaining(A,c,sort,vec):
    if c==[]:
        return sort
    else:
        if vec[c[-1]]>0:
            return inv_sortable_remaining(A,[c[-1]]+c[:-1],sort+[c[-1]],reflect(A,vec,c[-1]))
        else:
            return inv_sortable_remaining(A,c[:-1],sort,vec)


'''
#This was just a test of how python works
def even(n):
    for i in range(n):
        if i==(i//2)*2:
            yield(i)
'''





'''
# An old version (failed in B3):
def K(c,v):  # A set of sequences depending on a sortable element v, given by its c-sorting word
    out = []
    if c==[]:   # The terminal case (empty sequence for the rank-0 Coxeter group)
        out=[[]]
    elif v==[] or c[0]!=v[0]:  # The "induction on rank" case.
        oldK=K(c[1:],v)
        #print([c],oldK)
        for k in oldK:
            out=out+[[c[0]]+k]
        out=oldK+out
    else:  # The "induction on length" case.  This is provably the right thing to do.
        for k in K(c[1:]+[c[0]],v[1:]):
            if k!=[] and k[0]==v[0]:
                out=out+[k[1:]]
            else:
                out=out+[[v[0]]+k]
    return out

# Another old version (also failed in B3)
def K(c,v):  # A set of sequences depending on a sortable element v, given by its c-sorting word
    out = []
    if c==[]:   
        out=[[]]
    elif v==[] or c[0]!=v[0]:
        oldK=K(c[1:],v)
        for k in oldK:
            out=out+[k+[c[0]]]
        out=oldK+out
    else:
        for k in K(c[1:]+[c[0]],v[1:]):
            if k!=[] and k[0]==v[0]:
                out=out+[k[1:]]
            else:
                out=out+[[v[0]]+k]
    return out

#Another
def K(c,v):  # A set of sequences depending on a sortable element v, given by its c-sorting word
    out = []
    if c==[]:   # The terminal case (empty sequence for the rank-0 Coxeter group)
        out=[[]]
    elif v==[] or c[0]!=v[0]:  # The "induction on rank" case.
        for i in range(len(c)):
            out=out+[c[0:i+1]]
        out=K(c[1:],v)+out
    else:  # The "induction on length" case.  This is provably the right thing to do.
        for k in K(c[1:]+[c[0]],v[1:]):
            if k!=[] and k[0]==v[0]:
                out=out+[k[1:]]
            else:
                out=out+[[v[0]]+k]
    return out



# An experiment (a pretty dumb one... yields the list [v])
def K(c,v):  # A set of sequences depending on a sortable element v, given by its c-sorting word
    out = []
    if c==[]:   # The terminal case (empty sequence for the rank-0 Coxeter group)
        out=[[]]
    elif v==[] or c[0]!=v[0]:  # The "induction on rank" case.
        out=K(c[1:],v)
    else:  # The "induction on length" case.  This is provably the right thing to do.
        for k in K(c[1:]+[c[0]],v[1:]):
            if k!=[] and k[0]==v[0]:
                out=out+[k[1:]]
            else:
                out=out+[[v[0]]+k]
    return out

# old
def Kbip(cplus,cminus,Coxnum):
    out = [cminus+cplus,cminus,[]]
    for i in range(Coxnum-1):
        if (i//2)*2==i:  #even
            out=out+[out[-1]+cplus]
        else:
            out=out+[out[-1]+cminus]
    return out


'''
