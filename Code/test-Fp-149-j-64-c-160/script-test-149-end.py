load('/data/cyril/github/Code/isogeny-computing-testless.sage')
load('/data/cyril/github/Code/tate-module-extension-testless.sage')
load('/data/cyril/github/Code/extension_corps.sage')
#load('~/Documents/github-3/Code/isogeny-computing-testless.sage')
#load('~/Documents/github-3/Code/tate-module-extension-testless.sage')
#load('~/Documents/github-3/Code/extension_corps.sage')

F=FiniteField(149)
K=Tower_two(F,1)
print K._top
E=EllipticCurve(j=F(64))
if E.cardinality()!=160:
	E=E.quadratic_twist()
print E
for end in range(2,24):
	deg=end**2
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	if calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)[0]==True:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init_bis(b,c,b2,c2,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init_bis(b,c,b2,c2,2,d,deg,e,f,g)',number=10,repeat=5,seconds=True,preparse=True)
		A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=10,repeat=5,seconds=True,preparse=True)		
	else:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)',number=10,repeat=5,seconds=True,preparse=True)		
		if calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)[0]==True:
			A=timeit('calcul_isogenie_step_bis(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=10,repeat=5,seconds=True,preparse=True)
		else:
			A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=10,repeat=5,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=10,repeat=5,seconds=True,preparse=True)
	Fichier = open('/data/cyril/github/Code/test-script-149.txt','a')
	#Fichier=open('~/Documents/github-3/Code/test-script-149.txt','a')
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()

L=[31,47,67]

for l in L:
	deg=l
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	if calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)[0]==True:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init_bis(b,c,b2,c2,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init_bis(b,c,b2,c2,2,d,deg,e,f,g)',number=10,repeat=5,seconds=True,preparse=True)
		A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=10,repeat=5,seconds=True,preparse=True)		
	else:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)',number=10,repeat=5,seconds=True,preparse=True)		
		if calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)[0]==True:
			A=timeit('calcul_isogenie_step_bis(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=10,repeat=5,seconds=True,preparse=True)
		else:
			A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=10,repeat=5,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=10,repeat=5,seconds=True,preparse=True)
	Fichier = open('/data/cyril/github/Code/test-script-149.txt','a')
	#Fichier=open('~/Documents/github-3/Code/test-script-149.txt','a')
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()


L=[5,7,19,41,59,71]
E3=EllipticCurve(j=F(63))
if E3.cardinality()!=160:
	E3=E3.quadratic_twist()

for l in L:
	deg=l
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E3,((16.0)/3)*deg,K,2,conservation=True)
	if calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)[0]==True:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init_bis(b,c,b2,c2,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init_bis(b,c,b2,c2,2,d,deg,e,f,g)',number=10,repeat=5,seconds=True,preparse=True)
		A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=10,repeat=5,seconds=True,preparse=True)		
	else:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)',number=10,repeat=5,seconds=True,preparse=True)		
		if calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)[0]==True:
			A=timeit('calcul_isogenie_step_bis(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=10,repeat=5,seconds=True,preparse=True)
		else:
			A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=10,repeat=5,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=10,repeat=5,seconds=True,preparse=True)
	Fichier = open('/data/cyril/github/Code/test-script-149.txt','a')
	#Fichier=open('~/Documents/github-3/Code/test-script-149.txt','a')
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()

