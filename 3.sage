
PP = [2,3,5,7,11]#,13]#, 17, 19]#, 23, 29, 31, 37]

prod = 1
for p in PP:
    prod = p*prod
    
print prod

for D in range(-1800, 0):
    fl = 1
    for p in PP:
        if kronecker(D,p)<>1:
            fl = 0
    if fl ==1:
        print D, factor(D)
        print D-1, factor(D-1)
        print D-4, factor(D-4)
        #print D+1, factor(D+1)
