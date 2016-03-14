def creation_list_interpolation(P,Q,o1,o2,Lambda_1,Lambda_2,k2,a,b,Pb,Qb,Tower):
	'''
	Function that tabulates the abscissas we want to interpolate from the 
	diagonal basis of the two curves, the eigenvalues of the Frobenius and
	the mapping coefficient

	Input:
	-P,Q two points that generate horizontal isogenies each one 
	associated to the eigenvalue Lambda_1 (resp. Lambda_2)
	-Pb,Qb form also an horizontal basis
	-a and b integers such that we want \phi(P)=a*Pb and \phi(Q)=b*Qb
	-Lambda_1, Lambda_2 the 2 eigenvalues associated to P and Q for 
	the Frobenius
	-o1 and o2 the multiplicative order of Lambda_1 (resp.Lambda_2)
	modulo the order of P and Q
	-k2 is the 2 valuation of the order of P and Q

	Output:
	The list of the abscissas we want to interpolate with the abscissas of
	the desired image and the coefficients of them in the diagonal basis.
	The representatives points are in different rows of the table according
	to the length of their orbit under the Frobenius action
	'''
	T=[]
	C=[]
	M=[0]*(o2-o1+1)
	#print 'k2',k2
	for i in range(2**k2): 
		for j in range(2**k2):
			if (i%2==1)or(j%2==1 ):
				C.append((i*P+j*Q)[0])
			else:
				C.append(False)
		T.append(C)

		C=[]
	#print "o1,o2",o1,o2
	#we consider all the orders o1<oi<=o2
	for h in range(o2-o1):
		for i in range(2**k2):
			for j in range(2**(k2-h)):
				if T[i][2**h*j]!=False and j%2==1:
					#we consider the elements in the lowest level tower where they are defined
					i1=Tower.pushD(T[i][2**h*j],h)
					j1=Tower.pushD((a*i*Pb+b*2**h*j*Qb)[0],h)
					#we store the abscissas of the points 
					#and their image, and also their 
					#coefficients in the basis to obtain
					#new images later
					C.append([i1,j1,i,2**h*j])
					T[i][2**h*j]=False
					#we remove also the opposite
					T[-i%2**k2][-2**h*j%2**k2]=False
					i1=i
					j1=2**h*j
					#we remove also the other points of
					#the orbit
					for r in range(2**(o2-h)):
						i1=(i1*Lambda_1)%2**k2
						j1=(j1*Lambda_2)%2**k2
						T[i1][j1]=False
						T[-i1%2**k2][-j1%2**k2]=False
			j=0
		M[h]=C		
		#each row of corrrespond to representatives of orbits of same 
		#length for the Frobenius
		C=[]
	#in case of we have the orders o1 and o2 equals
	if o1==o2:
		for i in range(2**k2):
			for j in range(2**k2):
				if T[i][j]!=False and (i%2==1 or j%2==1):
					i1=Tower.pushD(T[i][j],o2)#we consider the elements in the lowest level of the tower where they are defined	
					j1=Tower.pushD((a*i*Pb+b*j*Qb)[0],o2)					
					C.append([i1,j1,i,j])
					T[i][j]=False
					T[-i%2**k2][-j%2**k2]=False
					i1=i
					j1=j
					for r in range(2**o2):
						i1=(i1*Lambda_1)%2**k2
						j1=(j1*Lambda_2)%2**k2
						T[i1][j1]=False
						T[-i1%2**k2][-j1%2**k2]=False
			j=0
		M[0]=C
		C=[]	
	else:
	#we make the order o1<o2
		for i in range(2**k2):
			for j in range(2**(k2-o2+o1)):
				if T[i][2**(o2-o1)*j]!=False and (i%2==1):
					i1=Tower.pushD(T[i][2**(o2-o1)*j],o2-o1)
					j1=Tower.pushD((a*i*Pb+b*2**(o2-o1)*j*Qb)[0],o2-o1)
					C.append([i1,j1,i,2**(o2-o1)*j])
					T[i][2**(o2-o1)*j]=False
					T[-i%2**k2][-2**(o2-o1)*j%2**k2]=False
					i1=i
					j1=2**(o2-o1)*j
					for r in range(2**o1):
						i1=(i1*Lambda_1)%2**k2
						j1=(j1*Lambda_2)%2**k2
						T[i1][j1]=False
						T[-i1%2**k2][-j1%2**k2]=False
		M[o2-o1]=C	
	return	M

def initialisation_T(L,o2,Tower):
	'''
	Input:
	-L a list of abscissas without their image but with their coefficients
	-o2 the 2 adic order of the greatest eigenvalue for the Frobenius
	-Tower, the 2-adic tower we are working on 
	
	Output:
	-A list with the abscissas of points ot interpolate with their 
	polynomials T^{i} for the interpolation computation
	
	'''
	#each row of L corresponds to representatives points of orbits of same
	#length
	for i in range(len(L)):
		for l in L[i]:
			C=[]
			#we work in the lowest level possible of the 2-adic
			#tower
			l[0]=Tower.determination_level(l[0])[0]
			R=PolynomialRing(l[0].parent(),'x')
			P=R.gen()-l[0]
			for r in range(o2-i):
				P2=P
				#we make act the frobenius at the right order 
				#according to the row of L
				P3=frobenius_polynomial(P,Tower,Tower._base.cardinality()**(2**(o2-i-1-r)))
				P=P2*P3
				#at each time we combine a polynomial with its 
				#frobenius conjugate we reduce its base ring of
				#one level in the 2-adic tower
				P=red_pol_basis(P,Tower)
				C.append([P2,P3])
			P=red_pol_basis(P,Tower)
			C.append([P])
			P=red_pol_basis(P,Tower,True)
			DP=P.diff()(l[0])
			DP2=R(DP**(-1))[0]
			DP2=Tower.determination_level(DP2)[0]
			C.append(DP2)
			l.append(C)
	#in the end we return the Polynomials T^(i,j) associated to 
	#representatives vof points of differents orbits, with also T'(v). They 
	#are stored in the i-th row according to the length of their orbit for 
	#the Frobenius
	return L

def interpolation_L(L,o2,o1,Tower,Pb,Qb,a,b):
	'''
	Input:
	-L a list of abscissas with their image and their coefficients
	-o2 the 2 adic order of the greatest eigenvalue for the Frobenius
	-Tower, the 2-adic tower we are working on 
	-Pb,Qb two points that forms an horizontal basis of the curve
	-a,b two mapping coefficients for the interpolation

	Output:
	-A table of interpolation polynomials defined on F_q
	'''
	N=[]
	P=[]
	for i in range(len(L)-1):
		for l in L[i]:
			w=(l[2]*a*Pb+l[3]*b*Qb)[0]
			#we make sure that w is defined on the lowest level in 
			#the 2-adic tower
			w=Tower.determination_level(w)[0]
			w=w*l[-1][-1]
			w2=Tower.frobenius_computation(w,Tower._base.cardinality()**(2**(o2-1-i)))
			M=[w,w2]
			for r in range(o2-i-1):
				#since we might have l[-1][r+1][0] and 
				#l[-1][r][0] not defined on the same polynomial
				#ring we have to make sure of this with this test
				if l[-1][r+1][0].parent()!=l[-1][r][0].parent():
					A=lift_pol_basis(l[-1][r+1][0],Tower)
				else:
					A=l[-1][r+1][0]
				L1=M[0]*(A/l[-1][r][0])
				L2=M[1]*(A/l[-1][r][1])
				L0=(L1+L2).numerator()
				#at each time we recombine a polynomial with  
				#its frobenius conjugate we reduce its base 
				#field				
				L0=red_pol_basis(L0,Tower)
				#we make act again the Frobenius ont the 
				#polynomial obtained
				L01=frobenius_polynomial(L0,Tower,Tower._base.cardinality()**(2**(o2-i-2-r)))
				M=[L0,L01]
			if l[-1][o2-i][0].parent()!=l[-1][o2-i-1][0].parent():
				A=lift_pol_basis(l[-1][o2-i][0],Tower)
			else:
				A=l[-1][o2-i][0]		
			L1=M[0]*A/l[-1][o2-i-1][0]
			L2=M[1]*A/l[-1][o2-i-1][1]
			#we make sure that sage do not represent it as a 
			#rational fraction
			L0=(L1+L2).numerator()
			L0=red_pol_basis(L0,Tower)
			P.append([L0,l[-1][-2][0]])
			#il faudrait peut etre rajouter le polynome T
		N.append(P)
		P=[]
	i=len(L)-1
	#we check here if we need to make act the frobenius again
	if o1!=0:
		for l in L[len(L)-1]:
			w=(l[2]*a*Pb+l[3]*b*Qb)[0]
			w=Tower.determination_level(w)[0]
			w=w*l[-1][-1]
			w2=Tower.frobenius_computation(w,Tower._base.cardinality()**(2**(o2-1-i)))
			M=[w,w2]
			for r in range(o2-i-1):
				if l[-1][r+1][0].parent()!=l[-1][r][0].parent():
					A=lift_pol_basis(l[-1][r+1][0],Tower)
				else:
					A=l[-1][r+1][0]				
				L1=M[0]*A/l[-1][r][0]
				L2=M[1]*A/l[-1][r][1]
				L0=(L1+L2).numerator()
				L0=red_pol_basis(L0,Tower)
				L01=frobenius_polynomial(L0,Tower,Tower._base.cardinality()**(2**(o2-i-2-r)))
				M=[L0,L01]
			if l[-1][o2-i][0].parent()!=l[-1][o2-i-1][0].parent():
				A=lift_pol_basis(l[-1][o2-i][0],Tower)
			else:
				A=l[-1][o2-i][0]
			L1=M[0]*A/l[-1][o2-i-1][0]
			L2=M[1]*A/l[-1][o2-i-1][1]
			L0=(L1+L2).numerator()
			L0=red_pol_basis(L0,Tower)
			P.append([L0,l[-1][-2][0]])
	else:
		#in this loop we work without the frobenius since we consider
		#elements associated to o1==0
		for l in L[len(L)-1]:		
			w=(l[2]*a*Pb+l[3]*b*Qb)[0]
			w=w*l[-1][-1]
			P.append([w,l[-1][-2][0]])
	N.append(P)
	P=[]
	return N	

def mult_tableau_interpola(N,Tower):
	'''
	Input:
	-N a table of interpolation polynomials associated to the different 
	points of different orbits for the Frobenius

	Output:
	-An interpolation polynomial L which is the combination of all the
	computed ones in N, with the the T associated which is a product of all
	the computed ones in N
	-Lc a pre computed list for the Chinese Remainder techniques
	'''
	Lc=[]
	c=0
	#each row of N correspond to points of differents orbits of same 
	#order for the Frobenius, we start by the ones with the smallest orbits
	#thus by the last one of N with the bigger one the first one of N
	while len(N)>1:
		l=len(N[0])		
		for i in range(len(N[-1])):
			#we reduce eventually the base field of the polynomial 
			#ring to have the polynomials defined on the base field
			#of the 2-adic tower
			A=red_pol_basis(N[-1][i][0],Tower,True)
			B=red_pol_basis(N[-1][i][1],Tower,True)
			C=red_pol_basis(N[0][c%l][0],Tower,True)
			D=red_pol_basis(N[0][c%l][1],Tower,True)
			A,B,U,V=CRT(A,B,C,D)
			#we store the results in the first row since we wont 
			#remove it
			N[0][c%l][0]=A
			N[0][c%l][1]=B
			#we store the pre computations for the CRT in a list
			Lc.append([U,V])			
			c+=1
		#as we have worked with elements in the last row once we have 
		#taken them into account we delete them
		N=N[0:-1]
	#now we have only one row, we just have to combine them with chinese 
	#remainder techniques
	k=0	
	while 2**k< len(N[0]):
		k+=1
	for j in range(k):		 
		for i in range(2**(k-j-1)):
			A=red_pol_basis(N[0][i][0],Tower,True)
			B=red_pol_basis(N[0][i][1],Tower,True)
			C=red_pol_basis(N[0][-i-1][0],Tower,True)
			D=red_pol_basis(N[0][-i-1][1],Tower,True)
			A,B,U,V=CRT(A,B,C,D)
			N[0][i][0]=A
			N[0][i][1]=B			
			Lc.append([U,V])
			c+=1
		N[0]=N[0][0:2**(k-j-1)]
		#we reduce each iteration by 2 the size of the row
	return N[0][0][0],N[0][0][1],Lc
	#we return here the result of combinations of all L modulo T obtained
	#with chines remainder techniques
	
def mult_tableau_interpola_pre(N,Tower,Lc):
	'''
	Function that computes the interpolation polynomial and return also a
	precomputed list.

	Input:
	-N a table of interpolation polynomials associated to the different 
	points of different orbits for the Frobenius
	-Lc a pre computed list for the Chinese Remainder techniques

	
	Output:
	-An interpolation polynomial L which is the combination of all the
	computed ones in N, with the the T associated which is a product of all
	the computed ones in N
	'''	
	#this function do the same as "mult_tableau_interpola" it just uses
	#pre computations to speed up CRT 
	c=0	
	while len(N)>1:
		l=len(N[0])		
		for i in range(len(N[-1])):
			A=red_pol_basis(N[-1][i][0],Tower,True)
			B=red_pol_basis(N[-1][i][1],Tower,True)
			C=red_pol_basis(N[0][c%l][0],Tower,True)
			D=red_pol_basis(N[0][c%l][1],Tower,True)			
			A,B=CRTm(A,B,C,D,Lc[c][0],Lc[c][1])
			N[0][c%l][0]=A
			N[0][c%l][1]=B
			c+=1
		N=N[0:-1]
	k=0	
	while 2**k< len(N[0]):
		k+=1	
	for j in range(k):		 
		for i in range(2**(k-j-1)):
			A=red_pol_basis(N[0][i][0],Tower,True)
			B=red_pol_basis(N[0][i][1],Tower,True)
			C=red_pol_basis(N[0][-i-1][0],Tower,True)
			D=red_pol_basis(N[0][-i-1][1],Tower,True)			
			A,B=CRTm(A,B,C,D,Lc[c][0],Lc[c][1])
			N[0][i][0]=A
			N[0][i][1]=B			
			c+=1
		N[0]=N[0][0:2**(k-j-1)]
	return N[0][0][0],N[0][0][1]

def frobenius_polynomial(P,Tower,power):
	'''
	Input:
	-P a polynomial with coefficients in Tower in the polynomial ring PR
	-PR a polynomial ring with coefficients in tower
	-power an integer which is a power of the Frobenius

	Output:
	-P the polynomial P with the frobenius power applied on all
	 of his coefficients
	'''
	PR=PolynomialRing(P[0].parent(),'x')
	L=P.list()
	M=[]
	for l in L:
		M.append(Tower.frobenius_computation(l,power))
	P=PR(M)
	return P

def red_pol_basis(P,Tower,Total=False):
	'''
	Input:
	-P a polynomial with coefficients in Tower
	-Tower a 2-adic tower
	-Total is a boolean value set by default to False, it says if the 
	reduction has to be total id est until coefficients of the polynomial
	are elements of the base field of the 2-adic Tower

	Ouput:
	-P the same polynomial but with his coefficients defined one level down in the tower

	Example:
	sage: F=FiniteField(101)
	sage: K=Tower_two(F,3,'x')
	sage: K1=K._top
	sage: K.add_one_level()
	sage: K2=K._top
	sage: Pol2=PolynomialRing(K2,'x')
	sage: Pol1=PolynomialRing(K1,'y')
	sage: L=[]
	sage: L2=[]
	sage: for i in range(5):
	....:     a=K1.random_element()
	....:     L2.append(a)
	....:     L.append(K.lift1l(a))
	....:     
	sage: print L2,L
	[82*x1^7 + 41*x1^6 + 57*x1^5 + 54*x1^4 + 4*x1^2 + 68*x1 + 27, 46*x1^7
	 + 15*x1^6 + 86*x1^5 + 42*x1^4 + 91*x1^3 + 89*x1^2 + 3*x1 + 4, 10*x1^7
	 + 7*x1^6 + 30*x1^5 + 45*x1^4 + 23*x1^3 + 2*x1^2 + 93*x1 + 33, 61*x1^7
	 + 89*x1^6 + 90*x1^5 + 12*x1^4 + 8*x1^3 + 18*x1^2 + 6*x1 + 4, 20*x1^7 
	+ 44*x1^6 + 15*x1^5 + x1^4 + 95*x1^3 + 78*x1^2 + 52*x1 + 98] [82*x3^14
	 + 41*x3^12 + 57*x3^10 + 54*x3^8 + 4*x3^4 + 68*x3^2 + 27, 46*x3^14 + 
	15*x3^12 + 86*x3^10 + 42*x3^8 + 91*x3^6 + 89*x3^4 + 3*x3^2 + 4, 
	10*x3^14 + 7*x3^12 + 30*x3^10 + 45*x3^8 + 23*x3^6 + 2*x3^4 + 93*x3^2 + 
	33, 61*x3^14 + 89*x3^12 + 90*x3^10 + 12*x3^8 + 8*x3^6 + 18*x3^4 + 
	6*x3^2 + 4, 20*x3^14 + 44*x3^12 + 15*x3^10 + x3^8 + 95*x3^6 + 78*x3^4 +
	 52*x3^2 + 98]
	sage: P2=Pol2(L)
	sage: P2
	(20*x3^14 + 44*x3^12 + 15*x3^10 + x3^8 + 95*x3^6 + 78*x3^4 + 52*x3^2 + 
	98)*x^4 + (61*x3^14 + 89*x3^12 + 90*x3^10 + 12*x3^8 + 8*x3^6 + 18*x3^4 
	+ 6*x3^2 + 4)*x^3 + (10*x3^14 + 7*x3^12 + 30*x3^10 + 45*x3^8 + 23*x3^6 
	+ 2*x3^4 + 93*x3^2 + 33)*x^2 + (46*x3^14 + 15*x3^12 + 86*x3^10 + 
	42*x3^8 + 91*x3^6 + 89*x3^4 + 3*x3^2 + 4)*x + 82*x3^14 + 41*x3^12 + 
	57*x3^10 + 54*x3^8 + 4*x3^4 + 68*x3^2 + 27
	sage: red_pol_basis(P2,K)
	(20*x1^7 + 44*x1^6 + 15*x1^5 + x1^4 + 95*x1^3 + 78*x1^2 + 52*x1 + 
	98)*x^4 + (61*x1^7 + 89*x1^6 + 90*x1^5 + 12*x1^4 + 8*x1^3 + 18*x1^2 + 
	6*x1 + 4)*x^3 + (10*x1^7 + 7*x1^6 + 30*x1^5 + 45*x1^4 + 23*x1^3 + 
	2*x1^2 + 93*x1 + 33)*x^2 + (46*x1^7 + 15*x1^6 + 86*x1^5 + 42*x1^4 + 
	91*x1^3 + 89*x1^2 + 3*x1 + 4)*x + 82*x1^7 + 41*x1^6 + 57*x1^5 + 54*x1^4 
	+ 4*x1^2 + 68*x1 + 27
	sage: P1=Pol1(L2)
	sage: P1
	(20*x1^7 + 44*x1^6 + 15*x1^5 + x1^4 + 95*x1^3 + 78*x1^2 + 52*x1 + 98)*
	y^4 + (61*x1^7 + 89*x1^6 + 90*x1^5 + 12*x1^4 + 8*x1^3 + 18*x1^2 + 6*x1 
	+ 4)*y^3 + (10*x1^7 + 7*x1^6 + 30*x1^5 + 45*x1^4 + 23*x1^3 + 2*x1^2 + 
	93*x1 + 33)*y^2 + (46*x1^7 + 15*x1^6 + 86*x1^5 + 42*x1^4 + 91*x1^3 + 
	89*x1^2 + 3*x1 + 4)*y + 82*x1^7 + 41*x1^6 + 57*x1^5 + 54*x1^4 + 4*x1^2 
	+ 68*x1 + 27
	'''
	L=P.list()
	if Total==False:
		if (P[0].parent()==Tower._levels[1] and P.parent()==PolynomialRing(Tower._levels[1],'x')) or (P[0].parent()==Tower._levels[0] and P.parent()==PolynomialRing(Tower._levels[0],'x')):
			return P	
	else: 
		if P[0].parent()==Tower._levels[0]:
			if  P.parent()==PolynomialRing(Tower._levels[0],'x'):		
				return P
			else:
				Pol=PolynomialRing(Tower._levels[0],'x')
				L2=[]
				for l in L:
					L2.append(l)
				return Pol(L2)
	Pol=PolynomialRing(Tower.push1l(P[0]).parent(),'x')
	L2=[]
	for l in L:
		L2.append(Tower.push1l(l))
	return Pol(L2)

def lift_pol_basis(P,Tower):
	'''
	Input:
	-P a polynomial with coefficients in Tower
	-Tower a 2-adic tower
	
	Ouput:
	-P the same polynomial but with his coefficients defined one level up in the tower

	Example:
	sage: K=Tower_two(F,1)
	sage: F=FiniteField(101)
	sage: K.add_one_level()
	sage: PR=PolynomialRing(K._levels[1],'x')
	sage: P2=PR.random_element()
	sage: P2
	(61*x1 + 14)*x^2 + (29*x1 + 94)*x + 31*x1 + 11
	sage: lift_pol_basis(P2,K)
	(61*x3^2 + 14)*x^2 + (29*x3^2 + 94)*x + 31*x3^2 + 11
	'''
	Pol=PolynomialRing(Tower.lift1l(P[0]).parent(),'x')
	L2=[]
	L=P.list()
	for l in L:
		L2.append(Tower.lift1l(l))
	return Pol(L2)

def CRT(A,TA,B,TB):

	'''
	Input:
	-A an interpolation polynomial, a remainder for the CRT
	-TA a modulus polynomial
	-B an interpolation polynomial, a remainder for the CRT
	-TB a modulus polynomial

	Output:
	A polynomial P such that P=AmodTA et P=B mod TB and his mudlus S
	and also the Bezout coefficients
	'''
	#chinese remainder theorem
	#1st step find the modulus prime with each other
	PR=TA.parent()
	r0=TA; u0=1; v0=0; r1=TB; u1=0; v1=1;
	while (r1!=0):
		q,r2=r0.quo_rem(r1)
		(r0,u0,v0,r1,u1,v1)=(r1,u1,v1,r2,u0-q*u1,v0-q*v1)
	if r0.degree()>0: 
		#we have a problem normaly the modulus are prime with each other
		print 'problem r0.degree()>0, r0.degree(), TA.degree(), TB.degree()', r0.degree(), TA.degree(), TB.degree(), r0		
		return False	 
	else:#we already have modulus prime with each other
		u0=u0.quo_rem(r0)[0]; v0=v0.quo_rem(r0)[0];
		S1=TA*TB
		R1=(TA*u0*B+v0*TB*A).quo_rem(S1)[1]
		#we return also pre computed value for the next occurrence
		#since the modulus do not change at each iteration try
		return (R1,S1,TA*u0,TB*v0)

def CRTm(A,TA,B,TB,u0,v0):

	'''
	Input:
	-A an interpolation polynomial, a remainder for the CRT
	-TA a modulus polynomial
	-B an interpolation polynomial, a remainder for the CRT
	-TB a modulus polynomial
	-u0 the precomputed bezout coefficient multiplied already by TA
	-v0 the precomputed bezout coefficient multiplied already by TB

	Output:
	A polynomial P such that P=AmodTA et P=B mod TB and his mudlus S
	and also the Bezout coefficients
	'''
	S=TA*TB
	R=(A*v0+B*u0)%S
	return (R,S)

def calcul_isogenie(P1,Q1,P2,Q2,l,order,d,Lambda_1,Lambda_2,Tower,interpol=None):
	'''
	INPUT:
	-P1 a primitive l^order-torsion point generating an horizontal l^order-isogeny on the crater
	-P2 a primitive l^order-torsion point generating an horizontal l^order-isogeny on the crater
	-Q1 a primitive l^order-torsion point generating an horizontal l^order-isogeny on the crater
	-Q2 a primitive l^order-torsion point generating an horizontal l^order-isogeny on the crater
	-l an integer
	-order an integer
	-T the l^order division polynomial on the domain curve by the isogeny
	-d the degree of the isogeny we want to compute	
	
	OUTPUT:
	An isogeny between the curve which P1 belogns and the curve which P2 belongs
	'''
	#we set the interpolation coefficient mapping both to 1	
	i=1
	j=1
	#we set the boolean value to test if we have found the isogeny to 
	#false
	Test=False
	#we compute the 2-adic valuation of the eigenvalues for the Frobenius
	o1=valuation(mod(Lambda_1,l**order).multiplicative_order(),l)
	o2=valuation(mod(Lambda_2,l**order).multiplicative_order(),l)
	power=Tower._base.cardinality()
	#we just swap the eigenvalues to have the convention that o2>o1
	if o1>o2:
		(o2,o1,Lambda_2,Lambda_1,Q1,P1,Q2,P2)=(o1,o2,Lambda_1,Lambda_2,P1,Q1,P2,Q2)
	#we compute the abscissas we want to interpolate
	L=creation_list_interpolation(P1,Q1,o1,o2,Lambda_1,Lambda_2,order,i,j,P2,Q2,Tower)
	#we tabulate all the different T associated to different 
	#representatives of orbits	
	M=initialisation_T(L,o2,Tower)
	#we compute  all the different interpolation polynomials L associated
	#to the different representatives	
	N=interpolation_L(M,o2,o1,Tower,P2,Q2,i,j)
	#we combine all the different interpolation polynomials and we get the
	#polynomials A and TA to do a rational reconstruction, Lc is a 
	#precomputed list for the next occurrences
	A,TA,Lc=mult_tableau_interpola(N,Tower)
	R2=PolynomialRing(A[0].parent(),'x')
	#we test if the polynomials obtained are good one if so the function 
	#return the rational function otherwise Test stays equal to False
	Test=fonction_test_iso(A,TA,R2,d,Tower)
	#we change the value of the coefficients mapping
	i+=1
	while (Test==False and j<l**order):
		while(Test==False and i<l**order and j%l!=0):
			if i%l==0:
				i+=1
			#we compute new interpolation polynomials according to 
			#the new mapping coefficients
			N=interpolation_L(M,o2,o1,Tower,P2,Q2,i,j)
			A,TA=mult_tableau_interpola_pre(N,Tower,Lc)
			Test=fonction_test_iso(A,TA,R2,d,Tower)
			#print 'i,j',i,j		
			i+=1
		i=1
		j+=1
	if Test!=False:
		r=Test[1]
		Test=Test[0]
		phi=Test/Test.leading_coefficient()
		return r,phi,r/(phi**2)
		#The rest of the code is made to return an isogeny using Sage 
		#constructor for isogeny
		K=Tower._base	
		PR=PolynomialRing(Tower._base,'x')
		#we are obliged to defined a curve on the lowest level in 
		#the 2-adic tower to be able to use Sage contructor
		E22=EllipticCurve([K(P1.curve().a4()[0]),K(P1.curve().a6()[0])])
		L=[]
		#we just coerce our list of coefficients of phi in the base field				
		for l in phi.list():
			L.append(l[0])
		phi2=E22.isogeny(PR([phi.list()[0][0],1]),degree=d)
		Test=phi2
		return Test
	#if no isogeny has been found we return the of Test which is 
	#equal to False in that case
	return 	Test	

def test_square(P,Tower):
		'''
		Input:
		-P a polynomial defined over a 2-adic tower
		-Tower a 2-adic tower

		Output:
		A boolean value saying if P is a square and also the root 
		polynomial of P
		'''
		#we get rid off the leading coefficient to not work in a tower
		#extension for the root computing, we return then the leading
		#coefficient to take this modification into account
		Lc=P.leading_coefficient()
		P=P/Lc
		DP=P.diff()
		P2=P.gcd(DP)
		P3=P.quo_rem(P2**2)
		if P3[0].degree()==0 and P3[1]==0 and P3[0][0]==Tower._levels[0](1) and P2.degree()==(P.degree()/2):
			return [True,P2,P3[0][0]*Lc]
		else:
			return [False]

def fonction_test_iso(A,T,R,d,Tower):
		'''
		INPUT
		-A a polynomial which is the interpolation of the primitive 
		torsion points on the others primitive torsion points of the 
		other curve
		-T the division polynomial of the torsion points
		-d the degree of the isogeny
		-Tower a 2-adic tower under which the polynomials are defined
		
		OUTPUT
		If the iosgeny exists it returns the  rational reconstruction 
		of the 	d-isogeny that respects A, otherwise it returns the 
		boolean value False.
		'''
		q0,r0=T,A
		q1,r1=T.quo_rem(A)		
		u1=R(0)
		u0=R(1)
		v1=R(1)
		v0=R(0)
		test=False
		deg=-1
		while (deg<d-1 and r1!=0):
			i=u1			
			u1=u0-q1*u1
			u0=i
			i=v1
			v1=v0-q1*v1
			v0=i
			i=r1
			r2=r0
			(q1,r1)=r0.quo_rem(r1)
			r0=i
			deg=R(v0).degree()
		if (deg==d-1):
			#we test if we have a square
			test=test_square(v0,Tower)
			#if we have a square we return the elements of the
			#rational function
			if test[0]==True:
				#we take here into account account the leading
				#coefficient that we put aside in test_square		
				r2=r2/test[2]
				return test[1],r2 
			else :
				return test[0]
		else :
			return test

def Couveignes_algorithme(E1,E2,r,Tower):
		'''
		INPUT:
		-E1 an elliptic curve defined on a level of Tower r-isogenous 
		to E2
		-E2 an elliptic curve defined on a level of Tower r-isogenous 
		to E1
		-r an integer
		-Tower a 2-adic tower

		Output:
		The rational map on the abscissas of the isogeny of degree r 
		that goes from E1 to E2

		Example:
		sage: F=FiniteField(101)
		sage: E=EllipticCurve(j=F(28))
		sage: G3=E(0).division_points(3); G3
		[(0 : 1 : 0), (12 : 18 : 1), (12 : 83 : 1)]
		sage: E2=E.isogeny(G3[1]).codomain(); E2
		Elliptic Curve defined by y^2 = x^3 + 48*x + 66 over Finite 
		Field of size 101
		sage: K=Tower_two(F,1,'x')
		sage: Couveignes_algorithme(E,E2,3,K)
		Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 
		87*x + 77 over Finite Field of size 101 to Elliptic Curve 
		defined by y^2 = x^3 + 48*x + 66 over Finite Field of size 101


 		'''
		#wew fix a same base field for the 2 curves
		K=E1.base_field()
		#we compute the bound necessary to reach to have enough points
		#for the interpolation 		
		B=(16.0/3)*(r)
		if E2.base_field()!=K:
			raise TypeError('the curves must be defined on the same field')
		else :
			#we compute horizontal basis on the 2 curves
			E1,P1,Q1,k1,Lambda_1,Lambda_2,Tower=tate_module(E1,B,Tower,2,conservation=True)
			E2,P2,Q2,k2,Lambda_12,Lambda_22=tate_module(E2,B,Tower,2)
			#we make sure to have the points associated to the same eigenvalues
			if Lambda_1==Lambda_22 and Lambda_2==Lambda_12:
				P2,Q2=Q2,P2	
			elif Lambda_1 != Lambda_12 or Lambda_2 != Lambda_22 :
				raise TypeError("probleme valeur propres,Lambda_1, Lambda_12, Lambda_2, Lambda_22" ,Lambda_1, Lambda_12, Lambda_2, Lambda_22) 
			return calcul_isogenie(P1,Q1,P2,Q2,2,k1,r,Lambda_1,Lambda_2,Tower)
'''
---------------------------------------
The rest of the code is just made for Timing, it is just some code cut at some points to make timing at precise moments
'''

def calcul_isogenie_init(P1,Q1,P2,Q2,l,order,d,Lambda_1,Lambda_2,Tower):
	'''
	A truncated function of calcul_isogenie just made to do some timing on
	computations necessary for the initialisation of the entire computation.
	The function returns all the necessary data to pursue the computations.
	'''
		#we set the interpolation coefficient mapping both to 1	
	i=1
	j=1
	#we set the boolean value to test if we have found the isogeny to 
	#false
	Test=False
	#we compute the 2-adic valuation of the eigenvalues for the Frobenius
	o1=valuation(mod(Lambda_1,l**order).multiplicative_order(),l)
	o2=valuation(mod(Lambda_2,l**order).multiplicative_order(),l)
	power=Tower._base.cardinality()
	#we just swap the eigenvalues to have the convention that o2>o1
	if o1>o2:
		(o2,o1,Lambda_2,Lambda_1,Q1,P1,Q2,P2)=(o1,o2,Lambda_1,Lambda_2,P1,Q1,P2,Q2)
	#we compute the abscissas we want to interpolate
	L=creation_list_interpolation(P1,Q1,o1,o2,Lambda_1,Lambda_2,order,i,j,P2,Q2,Tower)
	#we tabulate all the different T associated to different 
	#representatives of orbits	
	M=initialisation_T(L,o2,Tower)
	#we compute  all the different interpolation polynomials L associated
	#to the different representatives	
	N=interpolation_L(M,o2,o1,Tower,P2,Q2,i,j)
	#we combine all the different interpolation polynomials and we get the
	#polynomials A and TA to do a rational reconstruction, Lc is a 
	#precomputed list for the next occurrences
	A,TA,Lc=mult_tableau_interpola(N,Tower)
	R2=PolynomialRing(A[0].parent(),'x')
	#we test if the polynomials obtained are good one if so the function 
	#return the rational function otherwise Test stays equal to False
	Test=fonction_test_iso(A,TA,R2,d,Tower)
	if Test!=False:
		print "computation of the isogeny almost done"
		return [True,Test]
	else:
		return M,Lc,P2,Q2,R2,o2,o1,Tower,d

def calcul_isogenie_init_bis(P1,Q1,P2,Q2,l,order,d,Lambda_1,Lambda_2,Tower):
	'''
	A truncated function of calcul_isogenie just made to do some timing on
	computations necessary for the initialisation of the entire computation.
	The function returns all the necessary data to pursue the computations.
	'''
		#we set the interpolation coefficient mapping both to 1	
	i=3
	j=1
	#we set the boolean value to test if we have found the isogeny to 
	#false
	Test=False
	#we compute the 2-adic valuation of the eigenvalues for the Frobenius
	o1=valuation(mod(Lambda_1,l**order).multiplicative_order(),l)
	o2=valuation(mod(Lambda_2,l**order).multiplicative_order(),l)
	power=Tower._base.cardinality()
	#we just swap the eigenvalues to have the convention that o2>o1
	if o1>o2:
		(o2,o1,Lambda_2,Lambda_1,Q1,P1,Q2,P2)=(o1,o2,Lambda_1,Lambda_2,P1,Q1,P2,Q2)
	#we compute the abscissas we want to interpolate
	L=creation_list_interpolation(P1,Q1,o1,o2,Lambda_1,Lambda_2,order,i,j,P2,Q2,Tower)
	#we tabulate all the different T associated to different 
	#representatives of orbits	
	M=initialisation_T(L,o2,Tower)
	#we compute  all the different interpolation polynomials L associated
	#to the different representatives	
	N=interpolation_L(M,o2,o1,Tower,P2,Q2,i,j)
	#we combine all the different interpolation polynomials and we get the
	#polynomials A and TA to do a rational reconstruction, Lc is a 
	#precomputed list for the next occurrences
	A,TA,Lc=mult_tableau_interpola(N,Tower)
	R2=PolynomialRing(A[0].parent(),'x')
	#we test if the polynomials obtained are good one if so the function 
	#return the rational function otherwise Test stays equal to False
	Test=fonction_test_iso(A,TA,R2,d,Tower)
	if Test!=False:
		print "computation of the isogeny almost done"
		return [True,Test]
	else:
		return M,Lc,P2,Q2,R2,o2,o1,Tower,d


def calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,d):
	'''
	A truncated part of calcul_isogenie that computes the step which is 
	repeated for the interpolation tries. Made only for timing here.
	'''
	i=3
	j=1
	Test=False
	#we compute new interpolation polynomials according to 
	#the new mapping coefficients
	N=interpolation_L(M,o2,o1,Tower,P2,Q2,i,j)
	A,TA=mult_tableau_interpola_pre(N,Tower,Lc)
	Test=fonction_test_iso(A,TA,R2,d,Tower)
	if Test!=False:
		print "computation of the isogeny almost done"
		return [True,Test]
	else:
		return [Test]

def calcul_isogenie_step_bis(M,Lc,P2,Q2,R2,o2,o1,Tower,d):
	'''
	A truncated part of calcul_isogenie that computes the step which is 
	repeated for the interpolation tries. Made only for timing here.
	'''
	i=1
	j=1
	Test=False
	#we compute new interpolation polynomials according to 
	#the new mapping coefficients
	N=interpolation_L(M,o2,o1,Tower,P2,Q2,i,j)
	A,TA=mult_tableau_interpola_pre(N,Tower,Lc)
	Test=fonction_test_iso(A,TA,R2,d,Tower)
	if Test!=False:
		print "computation of the isogeny almost done"
		return [True,Test]
	else:
		return [Test]
