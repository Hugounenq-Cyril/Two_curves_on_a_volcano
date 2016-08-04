import sage.schemes.elliptic_curves.weierstrass_morphism as wm


def antecedent_2_mutliplication(E,Q,Tower):
	'''
	Input:
	-E an elliptic curve defined which posseses at least 3 2-torsion 
	primitive points
	-Q a point which possesses a pre image by the multiplication by 2 in Tower
	-Tower the 2-adic tower under which Q and E are defined

	Output:
	-R a point of E such that 2*R=Q
	'''
	#we compute the 2 2isogenies from which we will compute pre image and
	#also the weierstrass ismorphism which permit us to have the same 
	#representation of points
	phi1,phi2,Pl=calcul_2_isogenie_antecedent(E,Tower)
	#we compute a pre image of the image of Q by the weierstrass 
	#isomorphism
	R=calcul_antecedent_2_isogenie(Pl(Q),phi2,Tower)
	if phi2(R)!=Pl(Q):
		R=-R
	#we store the value of R to be able to check later if we got the good
	#pre image
	R2=R
	R= calcul_antecedent_2_isogenie(R,phi1,Tower)
	if phi1(R)!=R2:
		R=-R
	if 2*R!=Q:
		R=E(R[0],-R[1])
	return  R

def calcul_2_isogenie_antecedent(E,Tower):
	'''
	Input:
	-E an elliptic curve defined which posseses at least 3 2-torsion 
	primitive points
	-Tower the 2-adic tower under which E is defined	

	Output:
	-ph1 a 2-isogeny that has E1 for domain
	-phi2 a 2-isogeny that has E2 for codomain
	-Pl a weierstrass isomorphism that goes from E to the codomain of phi2
	'''
	K=Tower._base
	R=PolynomialRing(E.base_field(),'x')
	A=K(E.a4()[0])
	B=K(E.a6()[0])
	#we work only on the base field for the construction of the isogeny
	#because Sage is not able to construct them on an other level of the
	#2 adic tower
	Ebis=EllipticCurve([K(A),K(B)])
	L=Ebis(0).division_points(2)
	P1=L[1]
	P2=L[2]
	phi1=E.isogeny(R([K(-P1[0]),1]),degree=2)
	Q2=phi1(P2)
	E1=phi1.codomain()	
	phi2=E1.isogeny(R([K(-Q2[0][0]),1]),degree=2)	
	E2=phi2.codomain()
	A2=K(E2.a4()[0])
	B2=K(E2.a6()[0])
	#we construct E2bis for the same reason as for Ebis
	E2bis=EllipticCurve([K(A2),K(B2)])
	#the next 2 lines are only for the computation of the Weierstass 
	#isomorphism
	Pm=wm.isomorphisms(Ebis,E2bis,True) 
	Pl=wm.WeierstrassIsomorphism(E, Pm ,phi2.codomain())
	return phi1,phi2,Pl
	

def calcul_antecedent_2_isogenie(Q,phi,Tower):
	'''
	Input:
	-Q the abscissa of a 2-torsion point
	-phi a 2-isogeny with Q that belongs to the codomain of phi
	-Tower the 2-adic tower on which the curve is defined	

	Output:
	Preimage of Q by phi
	'''
	E=phi.domain()
	#we get the rational maps for the x, since we are only intersted in 
	#the abscissas of the points		
	f=phi.rational_maps()[0]
	PR2=PolynomialRing(Q[0].parent(),'x')	
	a2=E.a4()
	b2=E.a6()
	K=Tower._base
	num=f.numerator()
	den=f.denominator()
	#we have to tell sage that this is a polynomial and not a rational function
	num=PR2([num.coefficients()[2],num.coefficients()[1],num.coefficients()[0]])
	den=PR2([den.coefficients()[1],den.coefficients()[0]])
	solv=num-den*Q[0]
	#now we can use the implemented fonction to compute square roots to solve 
	#this quadratic equation
	delta_sqrt=Tower.root_computing(solv[1]**2-4*solv[2]*solv[0])
	#we make sure they are defined in the same level of the 2-adic tower to 
	#be able to compare them	
	delta_sqrt,b=Tower.meeting(delta_sqrt,solv[1]**2-4*solv[2]*solv[0])
	while delta_sqrt**2!=b:
		delta_sqrt=Tower.root_computing(solv[1]**2-4*solv[2]*solv[0])
		delta_sqrt,b=Tower.meeting(delta_sqrt,solv[1]**2-4*solv[2]*solv[0])
	#to make sure that the square root is at the good level in the tower	
	delta_sqrt,solv1=Tower.meeting(delta_sqrt,solv[1])
	delta_sqrt,solv2=Tower.meeting(delta_sqrt,solv[2])
	#to make sure that we work in the good polynomial ring
	if delta_sqrt.parent()!=solv[2].parent():
		PR2=PolynomialRing(delta_sqrt.parent(),'x')	
	A2=PR2((-solv1+delta_sqrt)/(2*solv2))[0]
	#to make sure that the square root is at the good level in the tower		
	A2,a2=Tower.meeting(A2,a2)
	A2,b2=Tower.meeting(A2,b2)
	B=Tower.root_computing(A2**3+a2*A2+b2)
	A2,B=Tower.meeting(A2,B)
	if E.base_field()!=A2.parent():
		a2,A2=Tower.meeting(a2,A2)
		b2,A2=Tower.meeting(b2,A2)
		E=EllipticCurve([a2,b2])
	#we test if we have the good point or its oppposite
	if phi(E(A2,B))!=Q :
		return E(A2,-B)
	else :
		return E(A2,B)

def division_point_qbis(l,P,E,Tower):
	#this function is made to work on point with coordinates on the tower
	'''
	This function computes an l-division point of P if it exists on the 
	field where the elliptic curve E is defined. If there do not exist 
	l-division point of P the function returns the same point P. The choice
	of the l-division is arbitrary, the function returns the first 
	l-division point computed.

	INPUT:

	- l an integer
	- P a point that belongs to E
	- E an elliptic curve defined over a field
	- Tower the 2-adic tower under which E is defined
	
	OUTPUT:
	A l-division point of P if it exists, if it does not exist it returns 
	the point P.

	EXAMPLE:
	sage: k=FiniteField(101)
	sage: E1=EllipticCurve(j=k(65)).quadratic_twist();
	sage: P=E1.random_point(); P=2*P; P
	(52 : 74 : 1)
	sage: P.order()
	12
	sage: Q=division_point_q(2,P,E1); print Q,Q.order()
	(43 : 14 : 1) 24
	'''
	if l==2 :
		return antecedent_2_mutliplication(E,P,Tower)
	#not useful since we work only with 2 torsion	
	else :
        	ans=[]
        	nP=-P
        	P_is_2_torsion = (P == -P)
        	g = E._multiple_x_numerator(l) - P[0]*E._multiple_x_denominator(l)
		if P_is_2_torsion:
			if l % 2 == 0:
                		g = g.gcd(g.derivative())*Tower.root_computing(g.leading_coefficient())
			
            		else: #sert à rien en pratique....
                		g0 = g.variables()[0] - P[0]

                		g = g // g0
                		g = g.gcd(g.derivative())*Tower.root_computing(g.leading_coefficient())
                		g = g0*g       
		g=g/g.leading_coefficient()	
		L=cantor_zassenhauss_4(g,g.parent().gen(),Tower); i=0;
		a1, a2, a3, a4, a6 = E.ainvs()
		#print 'L',L
        	while (len(ans)==0 and i<len(L)):
            		x=L[i]; i=i+1;
           	#Q = E.lift_x(x) #a implanter soi meme pour eviter un test de racine carree qui ne marchera pas
		#x=E.base_field()(x)
		#x,j=Tower.determination_level(x)
			x=Tower.meeting2(x,a2)
        		f = ((x + a2) * x + a4) * x + a6
			ys = Tower.root_computing(f)
			x,ys=Tower.meeting(x,ys)
			while x.parent()!=E.base_field():
				x=Tower.lift1l(x)
				ys=Tower.lift1l(ys)
			f = ((x + a2) * x + a4) * x + a6
			Q=E.point([x, ys, x.parent()(1)])
            		nQ = -Q
       			lQ = l*Q
       			if P_is_2_torsion:
            			if lQ == P:
            				ans.append(Q)
                		if nQ != Q:
                	        	ans.append(nQ)
			else:
                		if lQ == P:
                	        	ans.append(Q)
                	    	elif lQ == nP:
					ans.append(nQ)
		if (i==len(L) and len(ans)==0):
			print 'pas bouger'
			return P
        	return ans[0]

def calcul_torsion_max(E,l):
	'''
	This function computes the maximal rationnal l-torsion and returns 
	the rational generator point(s) of this l-torsion with the power 
	of the l rationnal torsion of each point returned.

	INPUTS:

	-E an elliptic curve defined over a finite field and on the cyclic
	crater of the volcano of l-isogeny
	-l an integer

	OUTPUTS:
	The function returns the l-th power of the eventually two subgroups 
	with generator points.

	EXAMPLE:
	sage: k=FiniteField(101)
	sage: E1=EllipticCurve(j=k(65)).quadratic_twist();
	sage: [k1,P,k2,Q]=calcul_torsion_max(E1,2); print k1,P,P.order().factor(),"k2",k2,Q,Q.order().factor()
	3 (95 : 91 : 1) 2^3 k2 2 (82 : 47 : 1) 2^2
	'''
	#we first define the rank of the rational l-torsion
	M=[]
        g = E.division_polynomial(l)
        for x in g.roots(multiplicities=False):
      	    if E.is_x_coord(x):
      	        Q = E.lift_x(x)
                lQ = l*Q
                if lQ == E(0):
                    M.append(Q)
        L=M; j=1;
	if len(M) != 1 :
	#if the rational torsion is not of rank 1 then we determine two points
	# that do not generate the same group, for this we use the weil pairing	 
		if l==2:
			L=[L[0],L[j]]
		else :
			while  L[j].weil_pairing(L[0],l).multiplicative_order()!=l:
	            		j=j+1
	        	L=[L[0],L[j]]
		k1=1
		P=division_point_q(l,L[0],E)
		while ( P !=L[0]):
			L[0]=P
			P=division_point_q(l,L[0],E)
			k1=k1+1	
		k2=1
		Q=division_point_q(l,L[1],E)
		while ( Q !=L[1]):
			L[1]=Q
			Q=division_point_q(l,L[1],E)
			k2=k2+1
		if k1==k2:
		#we look here if we can go further, it is possible that one of 
		#the two generator points is the generator of a descending 
		#isogeny, then we have to get this point to one who is not a 
		#generator of a descending isogeny 
			R=P+Q
			j=1
			P=division_point_q(l,R,E)		
			while (l*P!=R and j<=l-1):
				R=R+Q
				P=division_point_q(l,R,E)
				j=j+1
			#here we assume that we are on a crater
			while ( l*P==R):
				R=P
				P=division_point_q(l,R,E)
				k1=k1+1
		if k2>k1:
			[k1,P,k2,Q]=[k2,Q,k1,P]
		return k1,P,k2,Q
	#in the following case we have made a bad choice for the curve since it
	#only has cyclic l-torsion
	else : 
		k1=1
		P=division_point_q(l,L[0],E)
		while (P!=L[0] ):
			L[0]=P
			k1=k1+1
			P=division_point_q(l,L[0],E)
		return k1,P

def fonction_ent_diago(P,Q,Tower,rs,Lambda_1, Lambda_2,power,h):
	'''
	Input
	-P a point of 2**rs torsion
	-Q a point of 2**rs torsion
	-Tower the 2-adic tower under which P and Q are defined
	-rs the 2th-power of the order of P and Q
	-Lambda_1 , Lambda_2 the values of the eigenvalues of the Frobenius
	for the 2**(r-1) torsion points.
	-power the power of the frobenius to compute
	-h the 2 adic value of Lambda_1 - Lambda_2

	Output:
	P, Q, Lambda_1, Lambda_1 such that \pi(P)=Lambda_1P \pi(Q)=Lambda_2Q	
	'''
	#we compute the coefficients of the matrix of the Frobenius in the 
	#basis P,Q 
	a,b,c,d=calcul_coeff_diagonalisation(P,Q,Tower,rs,Lambda_1,Lambda_2,power)
	#we return then a diagonalized basis P,Q for the Frobenius and the new
	#values of the Lambda
	return fonction_diagonalisation(P,Q,a,b,c,d,Lambda_1,Lambda_2 , rs,h);

def calcul_coeff_diagonalisation(P,Q,Tower,rs,Lambda_1, Lambda_2,power):
	'''
	Input:
	-P a point of 2**rs torsion
	-Q a point of 2**rs torsion
	-Tower the 2-adic tower under which P and Q are defined
	-r the 2th-power of the order of P and Q
	-Lambda_1 , Lambda_2 the values of the eigenvalues of the Frobenius
	for the 2**(rs-1) torsion points.
	-power the power of the frobenius to compute
	-h the 2-adic value of Lambda1 - Lambda2

	Output:
	-a,b,c,d boolean such that pi(P)=(Lambda_1 + a2**(r-1))P + b2**(r-1)Q
	pi(Q)= c2**(r-1)P + (d2**(r-1) + Lambda_1)Q 
	'''
	A=Tower.frobenius_computation(P[0],power)
	B=Tower.frobenius_computation(P[1],power)
	R=Lambda_1*P
	if (R[0]==A and R[1]==B):
		a=0; b=0;
	elif (A==((R+2**(rs-1)*P)[0]) and B==((R+2**(rs-1)*P)[1])):
		a=1; b=0;
	elif (A==((R+2**(rs-1)*Q)[0]) and B==((R+2**(rs-1)*Q)[1])):
		a=0; b=1;
	else:
		a=1; b=1;
	A=Tower.frobenius_computation(Q[0],power)
	B=Tower.frobenius_computation(Q[1],power)
	R=Lambda_2*Q	
	if (R[0]==A and R[1]==B):
		c=0; d=0;
	elif (A==((R+2**(rs-1)*Q)[0]) and B==((R+2**(rs-1)*Q)[1])):
		d=1; c=0;
	elif (A==((R+2**(rs-1)*P)[0]) and B==((R+2**(rs-1)*P)[1])):
		c=1; d=0;
	else:
		c=1; d=1;
	return a,b,c,d

def fonction_diagonalisation(P,Q,a,b,c,d,Lambda_1, Lambda_2 , rs,h):
	'''
	Input:
	-a,b,c,d are integers in {0,1}
	-P, Q are 2**r-primitive torsion points
	-Lambda_1 is an eigenvalue associated to P
	-Lambda_2 is an eigenvalue associated to Q
	-h the 2 adic value of Lambda_1 - Lambda_2

	Output:
	P, Q, Lambda_1, Lambda_1 such that \pi(P)=Lambda_1P \pi(Q)=Lambda_2Q
	'''
	
	if b==0 and c==0:
	#easiest case, the matrix is already diagonalized
		return P,Lambda_1+2**(rs-1)*a, Q, Lambda_2+2**(rs-1)*d
	else:
		if b*c!=1:
		#the matrix is already trigonalized
			if b==0:
				Q=2**(rs-1-h)*P+Q
				return P,Lambda_1+2**(rs-1)*a, Q, Lambda_2+2**(rs-1)*d
			if c==0:
				P=2**(rs-1-h)*Q+P
				return P,Lambda_1+2**(rs-1)*a, Q, Lambda_2+2**(rs-1)*d
		else:
		#the matrix is nor trigonalized nor diagonalized
			P=2**(rs-1-h)*Q+P
			#first we trigonalized the matrix
			Q=2**(rs-1-h)*P+Q
			#then we diagonalized the matrix
			return P,Lambda_1+2**(rs-1)*a, Q, Lambda_2+2**(rs-1)*d

def suite_calcul_torsion_max(E,P,Q,k1,k2,l,Tower,i):
	'''

	This function computes on an elliptic curve E defined over a finite 
	fied K the l rationnal torsion of E on K from the knowledge of two 
	generator points P,Q of the l rationnal torsion of E over a subfield of
	K, with P and Q two points of l^k1 , l^k2 torsion respectively. The 
	functions returns the l-th power (k1 and k2) of the order of each 
	generator point returned P , Q of the rational l-torsion on K.

	INPUTS
	- E an elliptic curve defined over a finite field K and on the cyclic 
	crater of a volcano of l-isogeny
	- P,Q two point defined over E and generators of the l rational 
	torsion of E over a subfield of K
	- l an integer
	- k1 an integer such that P is a point of l^k1 torsion
	- k2 an integer such that Q is a point of l^k1 torsion
	- i is an integer and is for the number of steps we want to do

	OUTPUT:
	The l-th power of the generator points with the generator points		

	EXAMPLE:
	sage: k=FiniteField(101)
	sage: E1=EllipticCurve(j=k(65)).quadratic_twist();
	sage: [k1,P,k2,Q]=calcul_torsion_max(E1,2); print k1,P,P.order().factor(),"k2",k2,Q,Q.order().factor()
	3 (95 : 91 : 1) 2^3 k2 2 (82 : 47 : 1) 2^2
	sage: J.<a>=k.extension(2,conway=True, prefix='z'); J
	Finite Field in a of size 101^2
	sage: E1b,P,Q=construction_lift_better(E1,k,J,P,Q)
	sage: suite_calcul_torsion_max(E1b,P,Q,k1,k2,2)
	(4, (85*a + 7 : 78*a + 29 : 1), 3, (83*a + 7 : 91*a + 39 : 1))
	'''
	#here we make sure to have P and Q of same torsion
	if k1!=k2:
		if k1>k2:
			P=2**(k1-k2)*P
			k1=k2
		else :
			Q=2**(k2-k1)*P
			k2=k1
	h=k1
	#we initialise the eigenvalue for the Frobenius
	Lambda_1=1
	Lambda_2=1
	for j in range(i):
		R=division_point_qbis(l,P,E,Tower)
		P=R
		R=division_point_qbis(l,Q,E,Tower)	
		Q=R
		#at each step we diagonlaized the basis and we update the value
		#of the eigenvalues
		P,Lambda_1,Q,Lambda_2=fonction_ent_diago(P,Q,Tower,j+k2+1,Lambda_1, Lambda_2,Tower._base.cardinality(),h)
		if Lambda_1==Lambda_2:
			h+=1
	return P,Lambda_1,Q,Lambda_2,k2+i,k2+i,h

def suite_calcul_torsion_max_2(E,P,Q,k1,k2,l,Tower,i,h,Lambda_1,Lambda_2):
	'''

	This function computes on an elliptic curve E defined over a finite 
	fied K the l rationnal torsion of E on K from the knowledge of two 
	generator points P,Q of the l rationnal torsion of E over a subfield of
	K, with P and Q two points of l^k1 , l^k2 torsion respectively. The 
	functions returns the l-th power (k1 and k2) of the order of each 
	generator point returned P , Q of the rational l-torsion on K.

	INPUTS
	- E an elliptic curve defined over a finite field K and on the cyclic 
	crater of a volcano of l-isogeny
	- P,Q two point defined over E and generators of the l rational 
	torsion of E over a subfield of K
	- l an integer
	- k1 an integer such that P is a point of l^k1 torsion
	- k2 an integer such that Q is a point of l^k1 torsion
	- i is an integer and is for the number of steps we want to do
	- h is the 2 adic valuation of the difference between Lambda_1 and Lambda_2

	OUTPUT:
	The l-th power of the generator points with the generator points		

	EXAMPLE:
	sage: k=FiniteField(101)
	sage: E1=EllipticCurve(j=k(65)).quadratic_twist();
	sage: [k1,P,k2,Q]=calcul_torsion_max(E1,2); print k1,P,P.order().factor(),"k2",k2,Q,Q.order().factor()
	3 (95 : 91 : 1) 2^3 k2 2 (82 : 47 : 1) 2^2
	sage: J.<a>=k.extension(2,conway=True, prefix='z'); J
	Finite Field in a of size 101^2
	sage: E1b,P,Q=construction_lift_better(E1,k,J,P,Q)
	sage: suite_calcul_torsion_max(E1b,P,Q,k1,k2,2)
	(4, (85*a + 7 : 78*a + 29 : 1), 3, (83*a + 7 : 91*a + 39 : 1))
	'''
	R=division_point_qbis(l,P,E,Tower)
	P=R
	R=division_point_qbis(l,Q,E,Tower)	
	Q=R
	#il faut redresser la base et actualiser les valeurs propres
	P,Lambda_1,Q,Lambda_2=fonction_ent_diago(P,Q,Tower,k2+1,Lambda_1, Lambda_2,Tower._base.cardinality(),h)
	if Lambda_1==Lambda_2:
		h+=1
	return P,Lambda_1,Q,Lambda_2,k2+i,k2+i,h

def construction_lift_better(E1,k,K,P,Q): 
	'''
	This function lifts an elliptic curve E1 defined over the field k into 
	an elliptic curve E2 defined over the extension field K of k. The 
	function lifts also the point P, Q defined on E1 into points defined on
	E2. This function returns the elliptic curve E2 and the points P,Q who 
	belong to E2.


	INPUTS:
	-E1 an elliptic curve defined over k
	-k a subfield of K
	-K a field
	-P a point defined over E1
	-Q a point defined over E1

	OUTPUTS:
	A lift E2 of the elliptic curve E1 defined over K and the lift of the 
	points P,Q belonging to E2.

	EXAMPLE:
	sage: k=FiniteField(101)
	sage: E1=EllipticCurve(j=k(65)).quadratic_twist(); E1
	Elliptic Curve defined by y^2 = x^3 + 99*x + 54 over Finite Field of si
	ze 101
	sage: J.<a>=k.extension(2,conway=True, prefix='z'); J
	Finite Field in a of size 101^2
	sage: E1b,R,S=construction_lift_better(E1,k,J,R,S); E1b
	Elliptic Curve defined by y^2 = x^3 + 94*x + 6 over Finite Field in a o
	f size 101^2
	sage: R=E1.random_point(); R
	(92 : 59 : 1)
	sage: S=E1.random_point(); S
	(2 : 0 : 1)
	sage: E1b,R,S=construction_lift_better(E1,k,J,R,S); E1b
	Elliptic Curve defined by y^2 = x^3 + 94*x + 6 over Finite Field in a o
	f size 101^2
	sage: print R,S, R.curve(), R.curve()==S.curve()
	(92 : 59 : 1) (2 : 0 : 1) Elliptic Curve defined by y^2 = x^3 + 94*x + 
	6 over Finite Field in a of size 101^2 True	
	'''

        a=E1.a4()
	b=E1.a6()
        a=k(a)
        b=k(b)
	E2=EllipticCurve(K,[a,b])
	P=E2(K(P[0]),K(P[1]))
	Q=E2(K(Q[0]),K(Q[1]))
	return E2,P,Q

def construction_lift_betterbis(E1,k,K,P,Q,Tower): 
	'''
	This function lifts an elliptic curve E1 defined over the field k into 
	an elliptic curve E2 defined over the extension field K of k. The 
	function lifts also the point P, Q defined on E1 into points defined on
	E2. This function returns the elliptic curve E2 and the points P,Q who 
	belong to E2.


	INPUTS:
	-E1 an elliptic curve defined over k
	-k a subfield of K
	-K a field
	-P a point defined over E1
	-Q a point defined over E1

	OUTPUTS:
	A lift E2 of the elliptic curve E1 defined over K and the lift of the 
	points P,Q belonging to E2.

	EXAMPLE:
	sage: k=FiniteField(101)
	sage: E1=EllipticCurve(j=k(65)).quadratic_twist(); E1
	Elliptic Curve defined by y^2 = x^3 + 99*x + 54 over Finite Field of si
	ze 101
	sage: J.<a>=k.extension(2,conway=True, prefix='z'); J
	Finite Field in a of size 101^2
	sage: E1b,R,S=construction_lift_better(E1,k,J,R,S); E1b
	Elliptic Curve defined by y^2 = x^3 + 94*x + 6 over Finite Field in a o

	f size 101^2
	sage: R=E1.random_point(); R
	(92 : 59 : 1)
	sage: S=E1.random_point(); S
	(2 : 0 : 1)
	sage: E1b,R,S=construction_lift_better(E1,k,J,R,S); E1b
	Elliptic Curve defined by y^2 = x^3 + 94*x + 6 over Finite Field in a o
	f size 101^2
	sage: print R,S, R.curve(), R.curve()==S.curve()
	(92 : 59 : 1) (2 : 0 : 1) Elliptic Curve defined by y^2 = x^3 + 94*x + 
	6 over Finite Field in a of size 101^2 True	
	'''
	if k.cardinality().is_prime():
		a=E1.a4()
		b=E1.a6()
	else :
        	a=E1.a4()[0]
		b=E1.a6()[0]
	r=K.random_element()
	E2=EllipticCurve(K,[Tower._base(a),Tower._base(b)])
	P0=Tower.meeting2(P[0],r)
	Q0=Tower.meeting2(Q[0],r)
	P1=Tower.meeting2(P[1],r)
	Q1=Tower.meeting2(Q[1],r)
	P=E2(P0,P1)
	Q=E2(Q0,Q1)
	return E2,P,Q

def centering_frob(E,P,Q,Tower,l,o,Lambda_1,h,stair=None):
	'''	
	From a basis P,Q of the l^k torsion of the elliptic curve E it returns
	a point P ensuring that P is associated to the first eigenvalue of the 
	Frobenius endomorphism.


	INPUTS: 
	-E an elliptic curve defined over an extension field of K
	-P a point of E of order l**o
	-Q a point of E of order l**o
	-Tower a 2-adic tower under which P,Q,E are defined
	-l an integer	
	-o the power of a prime number denoted l	
	-Lambda_1 the eigenvalue for the frobenius associated to P
	-h the 2-adic valuation of Lambda_1-Lambda_2

	OUTPUTS:
	The function returns a point P such that P is an eigenvector for the 
	first eigenvalue of the Frobenius endomorphism.

	EXAMPLE:
	sage: k=FiniteField(101)
	sage: J.<a>=k.extension(2,conway=True, prefix='z'); J
	Finite Field in a of size 101^2
	sage: E1=EllipticCurve(j=k(65)).quadratic_twist();
	sage: [k1,P,k2,Q]=calcul_torsion_max(E1,2)
	sage: E1b,P,Q=construction_lift_better(E1,k,J,P,Q)
	sage: [k1,P,k2,Q]=suite_calcul_torsion_max(E1b,P,Q,3,2,2)
	sage: P
	(74 : 56 : 1)
	sage: P=P+Q; P
	(22*a + 51 : 65*a + 68 : 1)
	sage: (P).xy()[0]==(P.xy()[0])^101
	False
	sage: P=centering_left(E1b,P,Q)
	sage: P.xy()[0] in k
	True
	sage: P
	(50 : 27 : 1)
	sage: P.order()
	8
	'''
	Pfx=Tower.frobenius_computation(P[0],Tower._base.cardinality())
	Pmx=Lambda_1*P
	Pfy=Tower.frobenius_computation(P[1],Tower._base.cardinality())
	if (Pfx != Pmx[0]) or (Pfy != Pmx[1])  :
		P=P+l**(o-1-h)*Q
	return P

def straightening_step(E,P,Q,l,k,Tower,Lambda_1,h,stair=None):
	'''
	Returns an elliptic curve E which is the codomain of the isogeny phi
	generated by the point l^(k-1)*P. It returns also two points P,Q of l^k 
	order on the elliptic curve E, whith P associated to the first eigenva-
	lue of the Frobenius endomorphism and Q associated to the second eigen-
	value of the Frobenius endomorphism.
	
	INPUTS:
	-E an elliptic curve defined over an extension field of K
	-P a point of E of order l^k associated to the first eigenvalue 
	Lambda_1
	-Q a point of E of order l^k associated to the second eigenvalue of the
	Frobenius endomorphism
	-l an integer
	-k an integer
	-Tower, the 2 adic tower under P,Q and E are defined
	-Lambda_1 the eigenvalue of the frobenius endomorphism modulo l**k 
	associated to the point P
	-h the 2 adic value of Lambda_1-Lambda_2 (eigenvalues of the forbenius
	endomorphism)

	OUTPUTS:
	An elliptic curve isogeneous to the input curve E on the crater of the 
	volcano of l-isogeny. This isogeneous curve is the codomain of the 
	l^k-isogeny associated to the eigenvector Lambda_1  of the Frobenius 
	endomorphism. It also outputs the list of the dual of the l-isogenies 
	we used between the input and the output curve, and a l**k torsion 
	points P on the new curve each associated to the eigenvalue Lambda_1 of
	the Frobenius endomorphism.

	EXAMPLE:
	sage: k=FiniteField(101)
	sage: J.<a>=k.extension(2,conway=True, prefix='z'); J
	Finite Field in a of size 101^2
	sage: E1=EllipticCurve(j=k(65)).quadratic_twist();
	sage: [k1,P,k2,Q]=calcul_torsion_max(E1,2)
	sage: E1b,P,Q=construction_lift_better(E1,k,J,P,Q)
	sage: [k1,P,k2,Q]=suite_calcul_torsion_max(E1b,P,Q,3,2,2)
	sage: P=centering_left(E1b,2*P,Q)
	sage: E,P,Q=right_step(E1b,P,Q,2,3,k)
	sage: E
	Elliptic Curve defined by y^2 = x^3 + 13*x + 34 over Finite Field in a of size 101^2
	sage: print P,P.order(),Q,Q.order(), Q.xy()[0]^101==(5*Q).xy()[0]
	(67 : 53 : 1) 8 (41*a + 33 : 46*a + 68 : 1) 8 True
	'''
	#we compute an l-torsion point which will allow us to compute later 
	#horizontal l-iosgeny
	Q2=l**(k-1)*Q
	L=[]
	R=PolynomialRing(E.base_field(),'x')
	for i in range(k-1):
		#we compute an l-torsion point which will allow us to compute 
		#later horizontal l-iosgeny
		S=l**(k-1)*P
		phi=E.isogeny(R([-S[0],1]),degree=2)
		E=phi.codomain()
		f=phi.rational_maps()[0]
		g=phi.rational_maps()[1]
		Q0=f(Q[0],Q[1])
		Q1=g(Q[0],Q[1])
		Q=E(Q0,Q1)
		Q2=phi(Q2)
		#here we compute the dual isogeny of phi
		psi=E.isogeny(R([-Q2[0],1]),degree=2)
		K=Tower._base
		#we compute also the weierstrass isomorphism to make all those
		#isogeny constructed compatible
		E2bis=EllipticCurve([K(phi.domain().a4()[0]),K(phi.domain().a6()[0])])	
		Ebis=EllipticCurve([K(psi.codomain().a4()[0]),K(psi.codomain().a6()[0])])
		Pm=wm.isomorphisms(Ebis,E2bis,True)
		E1=psi.codomain()
		E2=phi.domain()
		Pl=wm.WeierstrassIsomorphism(E1,Pm,E2)
		#here we stock the isogeny and it's weierstrass isomorphism
		L.append([psi,Pl])
		#we push P
		P=E(f(P[0],P[1]),g(P[0],P[1])) #or P=phi(P) but this writing don t work with sage
		#we compute a pre image of P		
		R2=division_point_qbis(l,P,E,Tower)		
		if 2*R2!=P:
			print 	'2*R2,phi(P)',2*R2,P
		P=R2
		#we diagonalise the basis with P		
		P=centering_frob(E,P,Q,Tower,l,k,Lambda_1,h,stair) 
	return E,P,L

def way_back(P,L):
	'''
	Returns the image of the point P by the successive isogenies contained 
	in the list L

	INPUTS:
	-P a point of E of order l^k associated to one of the two "eigenvalues"
	-L a list of horizontal 2-isogenies associated to the second eigenvalue

	OUTPUT:
	A point P on the elliptic curve who is the codomain of the isogeny
	generated by Q.

	EXAMPLE:
	sage: k=FiniteField(101)
	sage: J.<a>=k.extension(2,conway=True, prefix='z'); J
	Finite Field in a of size 101^2
	sage: E1=EllipticCurve(j=k(65)).quadratic_twist();
	sage: [k1,P,k2,Q]=calcul_torsion_max(E1,2)
	sage: E1b,P,Q=construction_lift_better(E1,k,J,P,Q)
	sage: [k1,P,k2,Q]=suite_calcul_torsion_max(E1b,P,Q,3,2,2)
	sage: Pm=2^(k1-k2)*P
	sage: Er,Pr,Qr=right_step(E1b,Pm,Q,2,3,k)
	sage: El,Pl,Ql=left_step(E1b,Pm,Q,2,3,k)
	sage: for i in range(k2-1):
	    Er,Pr,Qr=right_step(Er,Pr,Qr,2,3,k)
 	    El,Pl,Ql=left_step(El,Pl,Ql,2,3,k)
	sage: P=way_back_left(Er,Pr,Qr)
	sage: print P,P.order(),P.curve(),P.curve().j_invariant()==E1b.j_invari
	ant()
	(69 : 88 : 1) 8 Elliptic Curve defined by y^2 = x^3 + 12*x + 92 over Fi
	nite Field in a of size 101^2 True	
	'''
	n=len(L)
	for i in range(n):
		phi,Pl=L[-i-1]
		E=phi.codomain()
		f=phi.rational_maps()[0]
		g=phi.rational_maps()[1]
		P0=f(P[0],P[1])
		P1=g(P[0],P[1])
		P=E(P0,P1)#or P=phi(P) but this don t work with sage
		P=Pl(P)
	return P

def tate_module(E,b,Tower,l,conservation=None):
	'''
	Return an horizontal basis (P,Q) of E with P,Q of order l**k such that
	l**{2k}>b. 
	
	INPUTS:
	-E an ellptic curve defined over a finite field K and such that E is 
	on the cyclic crater of a volcano of l-isogeny
	-b an integer bound 
	-Tower a 2-adic tower over a prime field
	-l a prime number
	-conservation a boolean variable to know if the algorithm outputs the 
	2-adic tower under the algorithm has worked, useful to have compatible
	representations of elements of a same finite field

	OUTPUTS:
	The same elliptic curve as the input curve E with an horizontal basis 
	(P,Q) of E with P,Q of order l**k such that l**{2k}>b.  and if 
	conservation is true it returns the tower under which are defined 
	points of the elliptic curve E.
	
	EXAMPLE:
	sage: k=FiniteField(101)
	sage: J.<a>=k.extension(2,conway=True, prefix='z'); J
	Finite Field in a of size 101^2
	sage: E1=EllipticCurve(j=k(65)).quadratic_twist();
	sage: E1b=construction_lift_bis(E1,k,J); E1b
	Elliptic Curve defined by y^2 = x^3 + 26*x + 68 over Finite Field in a 
	of size 101^2
	sage: [E,P,Q,k2]=tate_modulus(E1b,68,J,2); print E,P,Q, k2;
	Elliptic Curve defined by y^2 = x^3 + 42*x + 1 over Finite Field in c o
	f size 101^8 (21*c^7 + 41*c^6 + 23*c^5 + 19*c^4 + 87*c^3 + c^2 + 46*c +
	 75 : 21*c^7 + 60*c^6 + 40*c^5 + 14*c^4 + 6*c^3 + 44*c^2 + 68*c + 41 : 
	1) (31*c^7 + 89*c^6 + 67*c^5 + 11*c^4 + 46*c^2 + 57*c + 17 : 82*c^7 + 1
	9*c^6 + 80*c^5 + 73*c^4 + 56*c^3 + 99*c^2 + 87*c + 61 : 1) 5
	sage: P.curve()==Q.curve()
	True
	sage: print P.order().factor(), Q.order().factor()
	2^5 2^5
	sage: P.weil_pairing(Q,2^k2).multiplicative_order()
	32

	sage: F=FiniteField(101)
	sage: E=EllipticCurve(j=F(28))
	sage: K=Tower_two(F,1,'x') #you will need the other file with 2-adic tower construction
	sage: tate_modulus(E,73,K,2)
	(Elliptic Curve defined by y^2 = x^3 + 87*x + 77 over 
	Univariate Quotient Polynomial Ring in x5 over Finite Field of 		
	size 101 with modulus x5^8 + 15,
	 (99*x5^6 + 56*x5^4 + 94*x5^2 + 24 : 84*x5^6 + 29*x5^4 + 		
	30*x5^2 + 58 : 1),
	 (74*x5^4 + 19 : 42*x5^4 + 33 : 1),
 	4,
 	3,
 	7)


	'''
	K=E.base_field()
	if Tower.cardinality_field(Tower._levels[1])!=Tower.cardinality_field(Tower._levels[0])**2:
		#we made a Tower compatible with the one in 
		#input, with all the intermediary levels
		if (K!=Tower._base):
			K=FiniteField(Tower._base.cardinality())
			E=EllipticCurve(K,[E.a4()[0],E.a6()[0]])		
		i=valuation(valuation(Tower._top_cardinal,Tower._base.cardinality()),2)		
		Tower=Tower_two(E.base_field(),1,Tower._name,Tower._alpha)
		for r in range(i):
			Tower.add_one_level()
	
	#we set the eigenvalues of the Frobenius to 1	
	Lambda_1,Lambda_2=1,1
	#we compute in the 2 following blocs the rational torsion of the 
	#elliptic curve, it will be useful to know how many levels above
	#in the tower we will be working and it will permit us to have a basis
	# of the rational torsion of the elliptic curve that we will use to 
	#compute first diagonal basis and then horizontal basis 
	if (K==Tower._base):
		k1,P,k2,Q=calcul_torsion_max(E,l)
		ind=-2
		h=k2
	else :
		#we compute which level we are working in the tower
		ind=Tower.floor(K.random_element())
		#we compute a basis in the prime finite field 		
		K1=FiniteField(Tower._base.cardinality())
		E3=EllipticCurve(K1,[E.a4()[0],E.a6()[0]])
		k1,P,k2,Q=calcul_torsion_max(E3,l)
		#we convert the basis in the initial base field of the input 
		#curve 
		E,P,Q=construction_lift_better(E3,K1,E.base_field(),P,Q)
		#we evaluate how many times we will have to increase the 2-adic
		# value of the basis P,Q just before doing it
		i=valuation(valuation(Tower.cardinality_field(K),Tower._base.cardinality()),2)		
		P,Lambda_1,Q,Lambda_2,k1,k2,h=suite_calcul_torsion_max(E,P,Q,k1,k2,l,Tower,i)
		
	if k2==0:
		print("probleme choix courbe")
		raise TypeError('the curve must have some rational 2 torsion to work on!')
	#here we evaluate if the rational torsion is sufficient relatively to
	#the bound b
	if l**(2*k2)<b:
		i=1
		r=l**(2*k2)
		while (r*l**2<b):
			i=i+1
			r=r*l**2
		#K2=K.extension(2**i,conway=True, prefix='z',name='b')
		c=i
		if (E.base_field()==Tower._base):
			N=Tower._top_cardinal
			M=E.base_field().cardinality()
			O=N.valuation(M)
			M=O.valuation(2)
			if i>M:
				#we add the necessary levels in the tower 
				#relatively to the 2-torsion we want to work on
				while M<i:
					Tower.add_one_level()
					M+=1
			elif M>i:
				if Tower.cardinality_field(Tower._levels[1])!=Tower.cardinality_field(Tower._levels[0])**2:
				#we made a Tower compatible with the one in 
				#input, with all the intermediary levels		
					Tower=Tower_two(E.base_field(),1,Tower._name,Tower._alpha)
					for r in range(i):
						Tower.add_one_level()
			ind=-1
			#we store the highest level in the tower, where we will
			#work
			K2=Tower._levels[-1]
		else:
			while ind<-2 and i>0:
				i-=1 
				ind+=2
			if i>0:
				for j in range(i):
					Tower.add_one_level()
			#we store the highest level in the tower, where we will
			#work			
			K2=Tower._levels[ind]		
		#we express the base point and the elliptic curve in the level
		# in the tower where we will work accoring to the bound		
		E2,P,Q=construction_lift_betterbis(E,K,K2,P,Q,Tower)
		#print P,Q
		#we compute the diagonal basis of the frobenius with the 
		#eigenvalues
		P,Lambda_1,Q,Lambda_2,k1,k2,h=suite_calcul_torsion_max(E2,P,Q,k1,k2,l,Tower,c)
	else :
		P=l**(k1-k2)*P
		k1=k2
		K2=E.base_field()
		E2=E
	#if the eigenvalues are identical then we need to increase the torsion
	# we work on to reach a 2**(h+1) torsion
	while h>=k2:
		
		K=E2.base_field()
		if K==Tower._base:
			K2=Tower._levels[1]
			ind=Tower.floorfield(K2)		
		elif ind<-2:
			ind+=2
			K2=Tower._levels[ind]
		else:
			Tower.add_one_level()
			K2=Tower._levels[ind]		
		E2,P,Q=construction_lift_betterbis(E2,K,K2,P,Q,Tower)
		P,Lambda_1,Q,Lambda_2,k1,k2,h=suite_calcul_torsion_max_2(E2,P,Q,k1,k2,l,Tower,1,h,Lambda_1,Lambda_2)
	#we compute here diagonal points Pr,Pl on curves on the crater with 
	#list of horizontal isogenies that permits to go back to E and do not 
	#change the order of Pr and Pl but make their images horizontal	
	Er,Pr,Lr=straightening_step(E2,P,Q,l,k2,Tower,Lambda_1,h)
	El,Pl,Ll=straightening_step(E2,Q,P,l,k2,Tower,Lambda_2,h)
	#we compute the horizontal images of Pr and Pl
	P=way_back(Pr,Lr)
	Q=way_back(Pl,Ll)
	if conservation==None:
		return E2,P,Q,k2,Lambda_1,Lambda_2
	else:	
		return E2,P,Q,k2,Lambda_1,Lambda_2,Tower


def patch_not_on_volcano(E,b,Tower,l,conservation=None,path,ascending=True)
	'''
	Return a basis (P,Q) of E with P,Q of order l**k such that
	P and Q generates asending isogenies until they reach the crater, where
	they generate two horizontal isogenies of distincts directions. We have 
	also l**{2k}>b. 
	
	INPUTS:
	-E an ellptic curve defined over a finite field K and such that E is 
	on a volcano of l-isogeny with cyclic crater
	-b an integer bound 
	-Tower a 2-adic tower over a prime field
	-l a prime number
	-conservation a boolean variable to know if the algorithm outputs the 
	2-adic tower under the algorithm has worked, useful to have compatible
	representations of elements of a same finite field
	-path is a series 

	OUTPUTS:
	The same elliptic curve as the input curve E with an horizontal basis 
	(P,Q) of E with P,Q of order l**k such that l**{2k}>b.  and if 
	conservation is true it returns the tower under which are defined 
	points of the elliptic curve E.
	'''
	


'''
----------------------------------------------------------------
The rest of the functions were designed for l!=2 but not implemented yet 
efficiently with l-adic tower they were working with sage construction
'''
def division_point_q(l,P,E):
	'''
	This function is not use for the case l=2, it was made for the eventual
	cases of l different from 2 

	This function computes an l-division point of P if it exists on the 
	field where the elliptic curve E is defined. If there do not exist 
	l-division point of P the function returns the same point P. The choice
	of the l-division is arbitrary, the function returns the first 
	l-division point computed.

	INPUT:

	- l an integer
	- P a point that belongs to E
	- E an elliptic curve defined over a field
	
	OUTPUT:
	A l-division point of P if it exists, if it does not exist it returns 
	the point P.

	EXAMPLE:
	sage: k=FiniteField(101)
	sage: E1=EllipticCurve(j=k(65)).quadratic_twist();
	sage: P=E1.random_point(); P=2*P; P
	(52 : 74 : 1)
	sage: P.order()
	12
	sage: Q=division_point_q(2,P,E1); print Q,Q.order()
	(43 : 14 : 1) 24
	'''
	ans=[]
	nP=-P
        P_is_2_torsion = (P == -P)
        g = E._multiple_x_numerator(l) - P[0]*E._multiple_x_denominator(l)
        if P_is_2_torsion:
            if l % 2 == 0:
                g = g.gcd(g.derivative())*g.leading_coefficient().sqrt()

            else: #sert à rien en pratique....
                g0 = g.variables()[0] - P[0]
                g = g // g0
                g = g.gcd(g.derivative())*g.leading_coefficient().sqrt()
                g = g0*g
        L=g.roots(multiplicities=False); i=0;
        while (len(ans)==0 and i<len(L)):
            x=L[i]; i=i+1;
            if E.is_x_coord(x):
                Q = E.lift_x(x)
                nQ = -Q

                lQ = l*Q
                if P_is_2_torsion:
                    if lQ == P:
                        ans.append(Q)
                        if nQ != Q:
                            ans.append(nQ)
                else:
                    if lQ == P:
                        ans.append(Q)
                    elif lQ == nP:
                        ans.append(nQ)
	if (i==len(L) and len(ans)==0):
		return P
        return ans[0]


def fonction_test_point_redresse(P,l,k):
	'''
	Just a function to check if the point obtained of order l**k is a 
	generator of a l**k isogeny is horizontal
	
	Input:
	-P a point on an elliptic curve of order l**k
	-l a prime integer
	-k an integer

	Ouput:
	The list of the j-invariants obtained by computing the successive l**r
	isogenies with 1<=r<=k-1

	Example:
		
	'''
	
	E=P.curve()
	R=PolynomialRing(E.base_field(),x)	
	for i in range(k-1):
		S=2**(k-1-i)*P
		phi=E.isogeny(R([-S[0],1]),degree=2)
		P=phi(P)
		E=phi.codomain()
		print E.j_invariant()
	return E.j_invariant()
