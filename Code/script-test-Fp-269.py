#some script to test on some elliptic curves defined on Fp 269 located on a cyclic  with p=1 mod 4
load('/data/cyril/github/Code/isogeny-computing-testless.sage')
load('/data/cyril/github/Code/tate-module-extension-testless.sage')
load('/data/cyril/github/Code/extension_corps.sage')
F=FiniteField(269)
K=Tower_two(F,1)
print K._top
E=EllipticCurve(j=F(70))
E3=EllipticCurve(j=F(195))
if E.cardinality()!=E3.cardinality():
	E3=E3.quadratic_twist()
print E,E3
print Couveignes_algorithme(E,E3,83,K)
