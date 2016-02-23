load('/data/cyril/github/Code/isogeny-computing-testless.sage')
load('/data/cyril/github/Code/tate-module-extension-testless.sage')
load('/data/cyril/github/Code/extension_corps.sage')
#load('/home/cyril/Documents/github-3/Code/isogeny-computing-testless.sage')
#load('/home/cyril/Documents/github-3/Code/tate-module-extension-testless.sage')
#load('/home/cyril/Documents/github-3/Code/extension_corps.sage')
F=FiniteField(269)
K=Tower_two(F,1)
print K._top
E=EllipticCurve(j=F(70))
if E.cardinality()!=256:
	E=E.quadratic_twist()
print E
deg=71
a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
a,b2,c2,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
R=PolynomialRing(b[0].parent(),'x')
P2,Q2,Tower,L,M,Lc,Vr,o2,o1=calcul_isogenie_initialisation(b,c,b2,c2,R,2,d,d,deg,e,f,g)
#print Couveignes_algorithme(E,E,deg,Tower)
