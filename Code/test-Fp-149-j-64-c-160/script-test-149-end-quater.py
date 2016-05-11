load('../isogeny-computing-testless.sage')
load('../tate-module-extension-testless.sage')
load('../extension_corps.sage')

n=15
r=10
fichier='../../benchmarks/test-script-149-5.tsv'


F=FiniteField(149)
K=Tower_two(F,1)
print K._top
E=EllipticCurve(j=F(64))
if E.cardinality()!=160:
	E=E.quadratic_twist()
print E
for end in range(2,32):
	deg=end**2
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	if calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)[0]==True:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init_bis(b,c,b,c,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init_bis(b,c,b,c,2,d,deg,e,f,g)',number=n,repeat=r,seconds=True,preparse=True)
		A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)		
	else:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)',number=n,repeat=r,seconds=True,preparse=True)		
		if calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)[0]==True:
			A=timeit('calcul_isogenie_step_bis(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)
		else:
			A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=n,repeat=r,seconds=True,preparse=True)
	Fichier = open(fichier,'a')	
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()


for end in range(51,60):
	deg=end**2
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	if calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)[0]==True:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init_bis(b,c,b,c,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init_bis(b,c,b,c,2,d,deg,e,f,g)',number=n,repeat=r,seconds=True,preparse=True)
		A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)		
	else:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)',number=n,repeat=r,seconds=True,preparse=True)		
		if calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)[0]==True:
			A=timeit('calcul_isogenie_step_bis(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)
		else:
			A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=n,repeat=r,seconds=True,preparse=True)
	Fichier = open(fichier,'a')	
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()

for end in range(106,115):
	deg=end**2
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	if calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)[0]==True:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init_bis(b,c,b,c,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init_bis(b,c,b,c,2,d,deg,e,f,g)',number=n,repeat=r,seconds=True,preparse=True)
		A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)		
	else:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)',number=n,repeat=r,seconds=True,preparse=True)		
		if calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)[0]==True:
			A=timeit('calcul_isogenie_step_bis(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)
		else:
			A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=n,repeat=r,seconds=True,preparse=True)
	Fichier = open(fichier,'a')	
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()

for end in range(217,226):
	deg=end**2
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	if calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)[0]==True:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init_bis(b,c,b,c,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init_bis(b,c,b,c,2,d,deg,e,f,g)',number=n,repeat=r,seconds=True,preparse=True)
		A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)		
	else:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)',number=n,repeat=r,seconds=True,preparse=True)		
		if calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)[0]==True:
			A=timeit('calcul_isogenie_step_bis(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)
		else:
			A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=n,repeat=r,seconds=True,preparse=True)
	Fichier = open(fichier,'a')	
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()

