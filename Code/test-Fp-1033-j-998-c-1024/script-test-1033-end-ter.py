load('../isogeny-computing-testless.sage')
load('../tate-module-extension-testless.sage')
load('../extension_corps.sage')

F=FiniteField(1033)
K=Tower_two(F,1)
print K._top
E=EllipticCurve(j=F(998))
if E.cardinality()!=1024:
	E=E.quadratic_twist()
print E
for end in range(2,443):
	deg=end**2
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	if calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)[0]==True:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init_bis(b,c,b,c,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init_bis(b,c,b,c,2,d,deg,e,f,g)',number=3,repeat=5,seconds=True,preparse=True)
		A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=3,repeat=5,seconds=True,preparse=True)		
	else:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)',number=3,repeat=5,seconds=True,preparse=True)		
		if calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)[0]==True:
			A=timeit('calcul_isogenie_step_bis(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=3,repeat=5,seconds=True,preparse=True)
		else:
			A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=3,repeat=5,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=3,repeat=5,seconds=True,preparse=True)
	Fichier = open('../../benchmarks/test-script-1033-quater.tsv','a')
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()

E3=EllipticCurve(j=F(462))
if E3.cardinality()!=1024:
	E3=E3.quadratic_twist()
L=[7,37,43]
for l in L:
	deg=l
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E3,((16.0)/3)*deg,K,2,conservation=True)
	if calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)[0]==True:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init_bis(b,c,b2,c2,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init_bis(b,c,b2,c2,2,d,deg,e,f,g)',number=3,repeat=5,seconds=True,preparse=True)
		A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=3,repeat=5,seconds=True,preparse=True)		
	else:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)',number=3,repeat=5,seconds=True,preparse=True)		
		if calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)[0]==True:
			A=timeit('calcul_isogenie_step_bis(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=3,repeat=5,seconds=True,preparse=True)
		else:
			A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=3,repeat=5,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=3,repeat=5,seconds=True,preparse=True)
	Fichier = open('../../benchmarks/test-script-1033-quater.tsv','a')
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()

E3=EllipticCurve(j=F(613))
if E3.cardinality()!=1024:
	E3=E3.quadratic_twist()
L=[11,23,29,71]
for l in L:
	deg=l
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E3,((16.0)/3)*deg,K,2,conservation=True)
	if calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)[0]==True:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init_bis(b,c,b2,c2,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init_bis(b,c,b2,c2,2,d,deg,e,f,g)',number=3,repeat=5,seconds=True,preparse=True)
		A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=3,repeat=5,seconds=True,preparse=True)		
	else:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)',number=3,repeat=5,seconds=True,preparse=True)		
		if calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)[0]==True:
			A=timeit('calcul_isogenie_step_bis(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=3,repeat=5,seconds=True,preparse=True)
		else:
			A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=3,repeat=5,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=3,repeat=5,seconds=True,preparse=True)
	Fichier = open('../../benchmarks/test-script-1033-quater.tsv','a')
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()

E3=EllipticCurve(j=F(757))
if E3.cardinality()!=1024:
	E3=E3.quadratic_twist()
L=[3]
for l in L:
	deg=l
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E3,((16.0)/3)*deg,K,2,conservation=True)
	if calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)[0]==True:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init_bis(b,c,b2,c2,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init_bis(b,c,b2,c2,2,d,deg,e,f,g)',number=3,repeat=5,seconds=True,preparse=True)
		A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=3,repeat=5,seconds=True,preparse=True)		
	else:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init(b,c,b2,c2,2,d,deg,e,f,g)',number=3,repeat=5,seconds=True,preparse=True)		
		if calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)[0]==True:
			A=timeit('calcul_isogenie_step_bis(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=3,repeat=5,seconds=True,preparse=True)
		else:
			A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=3,repeat=5,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=3,repeat=5,seconds=True,preparse=True)
	Fichier = open('../../benchmarks/test-script-1033-quater.tsv','a')
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()