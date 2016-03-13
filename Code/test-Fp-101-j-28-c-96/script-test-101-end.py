load('/data/cyril/github/Code/isogeny-computing-testless.sage')
load('/data/cyril/github/Code/tate-module-extension-testless.sage')
load('/data/cyril/github/Code/extension_corps.sage')
#load('~/Documents/github-3/Code/isogeny-computing-testless.sage')
#load('~/Documents/github-3/Code/tate-module-extension-testless.sage')
#load('~/Documents/github-3/Code/extension_corps.sage')

F=FiniteField(101)
K=Tower_two(F,1)
print K._top
E=EllipticCurve(j=F(28))
if E.cardinality()!=96:
	E=E.quadratic_twist()
print E
for end in range(2,24):
	deg=end**2
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)
	A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=10,repeat=5,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=10,repeat=5,seconds=True,preparse=True)
	C=timeit('calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)',number=10,repeat=5,seconds=True,preparse=True)
	Fichier = open('/data/cyril/github/Code/test-script-101.txt','a')
	#Fichier=open('~/Documents/github-3/Code/test-script-101.txt','a')
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()
L=[3,31,41,47,71]
E3=EllipticCurve(j=F(65))
if E3.cardinality()!=96:
	E3=E3.quadratic_twist()
for l in L:
	deg=l
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E3,((16.0)/3)*deg,K,2,conservation=True)
	M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)
	A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=10,repeat=5,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=10,repeat=5,seconds=True,preparse=True)
	C=timeit('calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)',number=10,repeat=5,seconds=True,preparse=True)
	Fichier = open('/data/cyril/github/Code/test-script-101.txt','a')
	#Fichier=open('~/Documents/github-3/Code/test-script-101.txt','a')
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()
