# walking around a volcano

PHI2 = [[3, 0, 1], [0, 3, 1], [2, 0, -162000], [0, 2, -162000], [2, 1, 1488], [1, 2, 1488], [2, 2, -1], [1, 0, 8748000000], [0, 1, 8748000000], [1, 1, 40773375], [0, 0, -157464000000000]]
PHI3 = [[1, 0, 1855425871872000000000L], [0, 1, 1855425871872000000000L], [1, 1, -770845966336000000], [2, 0, 452984832000000], [0, 2, 452984832000000], [2, 1, 8900222976000], [1, 2, 8900222976000], [2, 2, 2587918086], [3, 0, 36864000], [0, 3, 36864000], [3, 1, -1069956], [1, 3, -1069956], [3, 2, 2232], [2, 3, 2232], [3, 3, -1], [4, 0, 1], [0, 4, 1]]
PHI4 = [[0, 0, 280949374722195372109640625000000000000L], [1, 0, -364936327796757658404375000000000000L], [0, 1, -364936327796757658404375000000000000L], [1, 1, -94266583063223403127324218750000L], [2, 0, 158010236947953767724187500000000L], [0, 2, 158010236947953767724187500000000L], [2, 1, 188656639464998455284287109375L], [1, 2, 188656639464998455284287109375L], [2, 2, 26402314839969410496000000L], [3, 0, -22805180351548032195000000000L], [0, 3, -22805180351548032195000000000L], [3, 1, 12519806366846423598750000L], [1, 3, 12519806366846423598750000L], [3, 2, -914362550706103200000L], [2, 3, -914362550706103200000L], [3, 3, 2729942049541120], [4, 0, 24125474716854750000L], [0, 4, 24125474716854750000L], [4, 1, 1194227244109980000], [1, 4, 1194227244109980000], [4, 2, 1425220456750080], [2, 4, 1425220456750080], [4, 3, 80967606480], [3, 4, 80967606480], [4, 4, 7440], [5, 0, -8507430000], [0, 5, -8507430000], [5, 1, 561444609], [1, 5, 561444609], [5, 2, -2533680], [2, 5, -2533680], [5, 3, 2976], [3, 5, 2976], [5, 4, -1], [4, 5, -1], [6, 0, 1], [0, 6, 1]]
PHI5 = [[0, 0, 141359947154721358697753474691071362751004672000L], [1, 0, 53274330803424425450420160273356509151232000L], [0, 1, 53274330803424425450420160273356509151232000L], [1, 1, -264073457076620596259715790247978782949376L], [2, 0, 6692500042627997708487149415015068467200L], [0, 2, 6692500042627997708487149415015068467200L], [2, 1, 36554736583949629295706472332656640000L], [1, 2, 36554736583949629295706472332656640000L], [2, 2, 5110941777552418083110765199360000L], [3, 0, 280244777828439527804321565297868800L], [0, 3, 280244777828439527804321565297868800L], [3, 1, -192457934618928299655108231168000L], [1, 3, -192457934618928299655108231168000L], [3, 2, 26898488858380731577417728000L], [2, 3, 26898488858380731577417728000L], [3, 3, -441206965512914835246100L], [4, 0, 1284733132841424456253440L], [0, 4, 1284733132841424456253440L], [4, 1, 128541798906828816384000L], [1, 4, 128541798906828816384000L], [4, 2, 383083609779811215375L], [2, 4, 383083609779811215375L], [4, 3, 107878928185336800], [3, 4, 107878928185336800], [4, 4, 1665999364600], [5, 0, 1963211489280], [0, 5, 1963211489280], [5, 1, -246683410950], [1, 5, -246683410950], [5, 2, 2028551200], [2, 5, 2028551200], [5, 3, -4550940], [3, 5, -4550940], [5, 4, 3720], [4, 5, 3720], [5, 5, -1], [6, 0, 1], [0, 6, 1]]

PHI = [ [], [], PHI2, PHI3, PHI4, PHI5  ]

TOOLONG = 200

q = next_prime(121043222281)  # q=121043222311301, 192, 861, 3
FF = FiniteField(q)
RR.<x> = PolynomialRing(FF)
print q

def Phieval(N, i, j, m):
    """ eval PHI_N(i,j) mod m  """
    v = 0
    for monomial in PHI[N]:
        v = ( v + (i**monomial[0]) * (j**monomial[1]) * (monomial[2]%m)  )%m
    return v

def onesidephi(N, i):
    """ input, N, i, output the polynomial f(y) = PHI_N(i,y)  """
    f = 0
    for monomial in PHI[N]:
        f = f + (i**monomial[0])*monomial[2]* x^monomial[1]
    return f.roots()

FC1 = 3  # the degree of isogeny

def buildedge():
    ls = []
    for i in range(q):
        g = onesidephi(FC1, i)
        ls.append( [ i,g ] )
    return ls

#LIST_JV = buildedge()

def moveahead(path):
    """ input: the existing path e.g. [a, b, c];  output [a, b, c, d]   """
    if len(path)>TOOLONG:
        return TOOLONG
#    edges = LIST_JV[ path[-1] ]  # all the edges from the last vertex
    edges = onesidephi(FC1, path[-1])
    for e in edges:
        v = e[0]
        if v not in path:
            longerpath = path + [v]
            #print longerpath
            return moveahead( longerpath )
        elif (v in path) and len(path)>1 and v!=path[-2]: # the normal way of forming a loop of length >2
            #print path
            return path + [v]
        elif len(path)>1 and v==path[-2] and e[1]==1: # walking backwards, skip and continue
            continue
        elif len(path)>1 and v==path[-2] and e[1]>1:  # a loop of length 2
            return path + [v]
        else:
#            print len(path), path, edge, v
            continue

def probe():
    """ probe the forwarding depth, excluding loops of length 2  """
    for source in range(1,1000): #LIST_JV[1:50]:
        if len( onesidephi(FC1, source) )>0:
            #print source
            path = [ source ]
            route = moveahead( path )
            if route!=TOOLONG and route!=None:
                print len(route)-1, source
                print route

probe()
