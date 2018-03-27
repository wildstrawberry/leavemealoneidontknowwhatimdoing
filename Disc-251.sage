# class group of K = QQ(\sqrt(-251))
# References:
# http://doc.sagemath.org/html/en/reference/curves/sage/schemes/elliptic_curves/ell_curve_isogeny.html
# http://www.math.uwaterloo.ca/~mrubinst/modularpolynomials/phi_l.html
# https://math.mit.edu/~drew/ClassicalModPolys.html
# https://eprint.iacr.org/2006/145.pdf

PHI2 = [[3, 0, 1], [0, 3, 1], [2, 0, -162000], [0, 2, -162000], [2, 1, 1488], [1, 2, 1488], [2, 2, -1], [1, 0, 8748000000], [0, 1, 8748000000], [1, 1, 40773375], [0, 0, -157464000000000]]
PHI3 = [[1, 0, 1855425871872000000000L], [0, 1, 1855425871872000000000L], [1, 1, -770845966336000000], [2, 0, 452984832000000], [0, 2, 452984832000000], [2, 1, 8900222976000], [1, 2, 8900222976000], [2, 2, 2587918086], [3, 0, 36864000], [0, 3, 36864000], [3, 1, -1069956], [1, 3, -1069956], [3, 2, 2232], [2, 3, 2232], [3, 3, -1], [4, 0, 1], [0, 4, 1]]
PHI4 = [[0, 0, 280949374722195372109640625000000000000L], [1, 0, -364936327796757658404375000000000000L], [0, 1, -364936327796757658404375000000000000L], [1, 1, -94266583063223403127324218750000L], [2, 0, 158010236947953767724187500000000L], [0, 2, 158010236947953767724187500000000L], [2, 1, 188656639464998455284287109375L], [1, 2, 188656639464998455284287109375L], [2, 2, 26402314839969410496000000L], [3, 0, -22805180351548032195000000000L], [0, 3, -22805180351548032195000000000L], [3, 1, 12519806366846423598750000L], [1, 3, 12519806366846423598750000L], [3, 2, -914362550706103200000L], [2, 3, -914362550706103200000L], [3, 3, 2729942049541120], [4, 0, 24125474716854750000L], [0, 4, 24125474716854750000L], [4, 1, 1194227244109980000], [1, 4, 1194227244109980000], [4, 2, 1425220456750080], [2, 4, 1425220456750080], [4, 3, 80967606480], [3, 4, 80967606480], [4, 4, 7440], [5, 0, -8507430000], [0, 5, -8507430000], [5, 1, 561444609], [1, 5, 561444609], [5, 2, -2533680], [2, 5, -2533680], [5, 3, 2976], [3, 5, 2976], [5, 4, -1], [4, 5, -1], [6, 0, 1], [0, 6, 1]]
PHI5 = [[0, 0, 141359947154721358697753474691071362751004672000L], [1, 0, 53274330803424425450420160273356509151232000L], [0, 1, 53274330803424425450420160273356509151232000L], [1, 1, -264073457076620596259715790247978782949376L], [2, 0, 6692500042627997708487149415015068467200L], [0, 2, 6692500042627997708487149415015068467200L], [2, 1, 36554736583949629295706472332656640000L], [1, 2, 36554736583949629295706472332656640000L], [2, 2, 5110941777552418083110765199360000L], [3, 0, 280244777828439527804321565297868800L], [0, 3, 280244777828439527804321565297868800L], [3, 1, -192457934618928299655108231168000L], [1, 3, -192457934618928299655108231168000L], [3, 2, 26898488858380731577417728000L], [2, 3, 26898488858380731577417728000L], [3, 3, -441206965512914835246100L], [4, 0, 1284733132841424456253440L], [0, 4, 1284733132841424456253440L], [4, 1, 128541798906828816384000L], [1, 4, 128541798906828816384000L], [4, 2, 383083609779811215375L], [2, 4, 383083609779811215375L], [4, 3, 107878928185336800], [3, 4, 107878928185336800], [4, 4, 1665999364600], [5, 0, 1963211489280], [0, 5, 1963211489280], [5, 1, -246683410950], [1, 5, -246683410950], [5, 2, 2028551200], [2, 5, 2028551200], [5, 3, -4550940], [3, 5, -4550940], [5, 4, 3720], [4, 5, 3720], [5, 5, -1], [6, 0, 1], [0, 6, 1]]
PHI9 = [[0, 0, 69044898341264155384039982812013450021641280641122350612180014989312000000000000000000000000000000000000L], [1, 0, -232683783909812135091180664408097870003379370875184432208292127703040000000000000000000000000000000000L], [0, 1, -232683783909812135091180664408097870003379370875184432208292127703040000000000000000000000000000000000L], [1, 1, -750477769359323701346534735710985354748615023238392518722583107993600000000000000000000000000000000L], [2, 0, 376501703401756663167829622430300312634214504891079031582441139601408000000000000000000000000000000L], [0, 2, 376501703401756663167829622430300312634214504891079031582441139601408000000000000000000000000000000L], [2, 1, 2180180133682453618437602678049589323532982602942504524433011926630400000000000000000000000000000L], [1, 2, 2180180133682453618437602678049589323532982602942504524433011926630400000000000000000000000000000L], [2, 2, -195104180586505976656170267263663673342537332095867469769673831088128000000000000000000000000L], [3, 0, -356498882808745744675952813139106931526561504326407178437871617966080000000000000000000000000000L], [0, 3, -356498882808745744675952813139106931526561504326407178437871617966080000000000000000000000000000L], [3, 1, -734745143710253598299658819584221973541388698584653375578080857292800000000000000000000000000L], [1, 3, -734745143710253598299658819584221973541388698584653375578080857292800000000000000000000000000L], [3, 2, -2174126850109403409612089909266604830941666101520869114646967552573440000000000000000000000L], [2, 3, -2174126850109403409612089909266604830941666101520869114646967552573440000000000000000000000L], [3, 3, 1337407428761009589661412615877939747792960714552280876322431791267840000000000000000000L], [4, 0, 209224067686996826205300712689558941749016250509807218578387753762816000000000000000000000000L], [0, 4, 209224067686996826205300712689558941749016250509807218578387753762816000000000000000000000000L], [4, 1, -552474127304059622933147654782708735487934827715739874235214567833600000000000000000000000L], [1, 4, -552474127304059622933147654782708735487934827715739874235214567833600000000000000000000000L], [4, 2, 1072273000502181967970881238481789529001284473312660386993634558345216000000000000000000L], [2, 4, 1072273000502181967970881238481789529001284473312660386993634558345216000000000000000000L], [4, 3, -66711604798252702203240814195082538435645647967276610575498938941440000000000000000L], [3, 4, -66711604798252702203240814195082538435645647967276610575498938941440000000000000000L], [4, 4, -18479552284868738690246708627837190388988153031268327611159810146304000000000000L], [5, 0, -71846975897804771463705349756964908464196982078331911573168361308160000000000000000000000L], [0, 5, -71846975897804771463705349756964908464196982078331911573168361308160000000000000000000000L], [5, 1, 243396300874543727943730913944913649785290815733841200853312530808832000000000000000000L], [1, 5, 243396300874543727943730913944913649785290815733841200853312530808832000000000000000000L], [5, 2, -59901709382329805978633635413667280365706226187233481397070735605760000000000000000L], [2, 5, -59901709382329805978633635413667280365706226187233481397070735605760000000000000000L], [5, 3, -18835180011463843231694024714228261766240954982932133312407483187200000000000000L], [3, 5, -18835180011463843231694024714228261766240954982932133312407483187200000000000000L], [5, 4, 658003637405245442250095582476112855345312171647903310493473832960000000000L], [4, 5, 658003637405245442250095582476112855345312171647903310493473832960000000000L], [5, 5, 81733165849983664671370721136696538031126996011599298838686334976000000L], [6, 0, 11840440831010017485683642912266874224079473751966006057244794290176000000000000000000L], [0, 6, 11840440831010017485683642912266874224079473751966006057244794290176000000000000000000L], [6, 1, -9317132816976477132813130144500667240517074482239080939579936604160000000000000000L], [1, 6, -9317132816976477132813130144500667240517074482239080939579936604160000000000000000L], [6, 2, 715289164764285891371697249337632347917412105898043988469312651264000000000000L], [2, 6, 715289164764285891371697249337632347917412105898043988469312651264000000000000L], [6, 3, 510602989311796036944510442921534158938486209261514004104006860800000000000L], [3, 6, 510602989311796036944510442921534158938486209261514004104006860800000000000L], [6, 4, 99649379447299227710013261724885046276640125416907777760194199552000000L], [4, 6, 99649379447299227710013261724885046276640125416907777760194199552000000L], [6, 5, 1287403158247585064152497034973072570105141673574412698991984640000L], [5, 6, 1287403158247585064152497034973072570105141673574412698991984640000L], [6, 6, -2003361981538502515786431353709655012603350293310522697051236L], [7, 0, 2893588694944278522177764102620581575000495540984238180610867200000000000000000L], [0, 7, 2893588694944278522177764102620581575000495540984238180610867200000000000000000L], [7, 1, -60495652821744993287486932849986202311614543306205957452432670720000000000000L], [1, 7, -60495652821744993287486932849986202311614543306205957452432670720000000000000L], [7, 2, 64735274511100762564459248803228674900735227692774967911928299520000000000L], [2, 7, 64735274511100762564459248803228674900735227692774967911928299520000000000L], [7, 3, -16474618115158934218802257796780294354859440459742602945869381632000000L], [3, 7, -16474618115158934218802257796780294354859440459742602945869381632000000L], [7, 4, 642148740736630948004496052690013708092074404011516770532392960000L], [4, 7, 642148740736630948004496052690013708092074404011516770532392960000L], [7, 5, -2064039409580538723682303423606644781559896017127825973607192L], [5, 7, -2064039409580538723682303423606644781559896017127825973607192L], [7, 6, 589799947861930821558414471517379624168655297448978560L], [6, 7, 589799947861930821558414471517379624168655297448978560L], [7, 7, 114208634776799415892612599617569985637849232080L], [8, 0, 235558341987396949911428581593686239877979341954146161917952000000000000L], [0, 8, 235558341987396949911428581593686239877979341954146161917952000000000000L], [8, 1, 6546983950323272622474592016942583982874006866244655026012160000000000L], [1, 8, 6546983950323272622474592016942583982874006866244655026012160000000000L], [8, 2, 8776475257933895754391596525840222413136747280735390811553792000000L], [2, 8, 8776475257933895754391596525840222413136747280735390811553792000000L], [8, 3, 1773899530331128588088312119839752844385302827595855128985600000L], [3, 8, 1773899530331128588088312119839752844385302827595855128985600000L], [8, 4, 61334584124724430829107139003055411562736252229635281920495L], [4, 8, 61334584124724430829107139003055411562736252229635281920495L], [8, 5, 295043918779312988862524085979914335281706535935143680L], [5, 8, 295043918779312988862524085979914335281706535935143680L], [8, 6, 114663187991652512877947591740259351959099261248L], [6, 8, 114663187991652512877947591740259351959099261248L], [8, 7, 2196008982690079369308616861054207114320L], [7, 8, 2196008982690079369308616861054207114320L], [8, 8, -169084890908576041794892443742536L], [9, 0, 6390980147531295015493344616502870354075036858198261760000000000L], [0, 9, 6390980147531295015493344616502870354075036858198261760000000000L], [9, 1, -7900333936192849023918427261965278932265209355223171072000000L], [1, 9, -7900333936192849023918427261965278932265209355223171072000000L], [9, 2, 3273266810212629480595452963053694318464393523934986240000L], [2, 9, 3273266810212629480595452963053694318464393523934986240000L], [9, 3, -527782836316123418691170962447078429119508813357952220L], [3, 9, -527782836316123418691170962447078429119508813357952220L], [9, 4, 29938980095729674278837381908388909886666835116800L], [4, 9, 29938980095729674278837381908388909886666835116800L], [9, 5, -452102708759835815999184660653014461675230688L], [5, 9, -452102708759835815999184660653014461675230688L], [9, 6, 1097815847178520649575574301039075207792L], [6, 9, 1097815847178520649575574301039075207792L], [9, 7, -169096306433121398819742262191810L], [7, 9, -169096306433121398819742262191810L], [9, 8, 205874310760628521421376L], [8, 9, 205874310760628521421376L], [9, 9, 28587961955951784], [10, 0, 10331567886902497628770879898357071872000000L], [0, 10, 10331567886902497628770879898357071872000000L], [10, 1, 6381231899147017430314467070087302021120000L], [1, 10, 6381231899147017430314467070087302021120000L], [10, 2, 155705417634012907024266501589913689446466L], [2, 10, 155705417634012907024266501589913689446466L], [10, 3, 655424730501203626951599797646911785920L], [3, 10, 655424730501203626951599797646911785920L], [10, 4, 680444811295518681180723971143182528L], [4, 10, 680444811295518681180723971143182528L], [10, 5, 186204831778242651626938540276560L], [5, 10, 186204831778242651626938540276560L], [10, 6, 11645320898401795868144158404L], [6, 10, 11645320898401795868144158404L], [10, 7, 102969059545961636573088L], [7, 10, 102969059545961636573088L], [10, 8, 28587961990122552], [8, 10, 28587961990122552], [10, 9, 15624], [9, 10, 15624], [10, 10, -1], [11, 0, 5567288717204029440000L], [0, 11, 5567288717204029440000L], [11, 1, -3894864835363363281932L], [1, 11, -3894864835363363281932L], [11, 2, 160958016085240175040L], [2, 11, 160958016085240175040L], [11, 3, -1807128632206069128], [3, 11, -1807128632206069128], [11, 4, 8462621974879728], [4, 11, 8462621974879728], [11, 5, -19911358807902], [5, 11, -19911358807902], [11, 6, 25558882848], [6, 11, 25558882848], [11, 7, -18155340], [7, 11, -18155340], [11, 8, 6696], [8, 11, 6696], [11, 9, -1], [9, 11, -1], [12, 0, 1], [0, 12, 1]]

PHI = [ [], [], PHI2, PHI3, PHI4, PHI5, [], [], [], PHI9 ]

PP = 83
QQ = 173
NN = PP*QQ
FFP = FiniteField(PP)
FFQ = FiniteField(QQ)
RingN = Integers(NN)
PolyRingN.<x> = PolynomialRing(RingN)
print PP, QQ, NN

db = HilbertClassPolynomialDatabase()
f251 = db[-251]

def PhiEval(N, i, j, m):
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
    return f #.roots()

S83 = [15, 48, 23, 29, 34, 55, 71]
S173 = [2, 162, 36, 117, 134, 116, 167]

for i in range(7):
    CRT1 = crt(S83[i], S173[i], PP, QQ)
    CRT2 = crt(S83[i-5], S173[i-5], PP, QQ)
    f1 = onesidephi(3, CRT1)
    f2 = onesidephi(5, CRT2)
    print i, CRT1, CRT2, f2.gcd(f1)

def buildedge():
    ls = []
    fc1 = 4
    for i in range(q):
        g = onesidephi(fc1, i)
        if len(g)>0:
            print i,fc1,g
        ls.append( [ i,g ] )
    return

#buildedge()

def finddots(q, FF):
    """ identify curves in specific isogeny classes """
    possiblejv = []
    for A in xrange(0,q):
        for B in xrange(0,q):
            try:
                E1 = EllipticCurve(FF, [0,0,0,A,B])
                #print E1.j_invariant()
                if E1.count_points(1)==153:
                    if E1.j_invariant() not in possiblejv:
                        possiblejv.append(E1.j_invariant())
                    #print "E1:", E1, "j(E1):", E1.j_invariant(), "#(E1)=", E1.count_points(1)
            except (ArithmeticError):
                continue
    print possiblejv

#finddots()
