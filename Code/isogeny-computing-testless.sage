def creation_matrice(P,Q,l,order,PR):
		#not useful since we have implemented an other function
		'''
		Functions that tabulates linear combinations of P and Q

		INPUT
		-P a primitive l^order torsion point on an elliptic curve E
		-Q a primitive l^order torsion point on an elliptic curve E
		-PR the polynomial ring under which further calculations will be made

		OUTPUT:
		A list of abscissas of multiples of the point P and Q
		'''
		L=[]
		R=P
		S=P.curve()(0)
		M=[]
		M.append(0)
		for j in range(0,l**(order-1)+1):			
			R+=S			
			for i in range(0,l**order-1):
				M.append(PR(R[0]))
				R+=P	
			L.append(M)
			R=P
			S+=Q
			M=[]
			M.append(S[0])
		return L
'''
def creation_list_interpolation(M1,M2,b,d,R):
		#not useful since we have implemented other fonctions to do this according to the frobenius
		
		INPUT:
		-M1 a 2 dimension list of abscissas of points of an elliptic curve
		-M2 a 2 dimension list of abscissas of points of an elliptic curve
		-b an integer prime to the order of the primitive points who spanned M1 and M2
		-d an integer prime to the order of the primitive points who spanned M1 and M2

		OUPUT:
		An interpolation list who has for first entry a point of M1 and second entry a
		point of M2 according to their index and b and d
		
		L=[]
		L0=len(M1[0])
		i=0
		for j in range(1,L0/2+1):
				L.append([M1[0][j],M2[0][s]])
		for i in range(1,L0/2):
			for j in range(0,L0):
				r=i*b%L0
				s=j*d%L0
				if r>L0/2:
					L.append([M1[i][j],M2[L0-r][-s]])
				else :
					L.append([M1[i][j],M2[r][s]])
		r=(L0/2)*b%L0
		for j in range(0,L0/2+1):
			s=j*d%L0
			if r>L0/2:
				L.append([M1[L0/2][j],M2[L0-r][-s]])
			else:
				L.append([M1[L0/2][j],M2[r][s]])
		return L
'''
'''
def modif_list_interpolation(L,M1,M2,b,d,R):
		#fonction inutile une autre a ete implante

		INPUT:
		-M1 a 2 dimension list of abscissas of points of an elliptic curve
		-M2 a 2 dimension list of abscissas of points of an elliptic curve
		-b an integer prime to the order of the primitive points who spanned M1 and M2
		-d an integer prime to the order of the primitive points who spanned M1 and M2

		OUPUT:
		An interpolation polynomial who sents the entries of M1 on the entries of M2 according to their index and b and d

		L0=len(M1[0])
		i=0
		c=0
		for j in range(1,L0/2+1):
				r=i*b%L0
				s=j*d%L0
				if r>L0/2:
					L[c][1]=M2[L0-r][-s]
				else :
					L[c][1]=M2[r][s]
				c+=1
		for i in range(1,L0/2):
			for j in range(0,L0):
				r=i*b%L0
				s=j*d%L0
				if r>L0/2:
					L[c][1]=M2[L0-r][-s]
				else :
					L[c][1]=M2[r][s]
				c+=1
		i=L0/2
		for j in range(0,L0/2+1):
			r=i*b%L0
			s=j*d%L0
			if r>L0/2:
				L[c][1]=M2[L0-r][-s]
			else:
				L[c][1]=M2[r][s]
			c+=1
		return L

'''
'''
def creation_liste(P,l,order):#devenu obsolete puisque l on cree directement un tableau a deux dimensions
		
		INPUT:
		-P a primitive l^order torsion point on an elliptic curve E
		-l an integer
		-order an integer

		OUPUT:
		A list of abscissas of multiples of the point P
		
		L=[P.xy()[0]]
		Q=P
		for i in range(1,l**order):
			L.append(Q.xy()[0])
			Q+=P
		return L
'''
'''
def creation_poly_interpolation_old(L1,L2,b,L3,L4,d,R): #devenu obsolete vu que l on travaille avec des tableaux de dimension 2
		
		INPUT:
		-L1 a list of abscissas of multiples of point corresponding to it's abscissa in the first entry 
		-L2 a list of abscissas of multiples of point corresponding to it's abscissa in the first entry 
		-R the polynomial ring in which we work

		OUPTUT:
		An interpolation polynomial that sends the elements of the first list L1 on elements of the second list according to the index
 		
		L=[]
		L0=len(L1)
		for i in range(1,L0/2+1):
			L.append([L1[i],L2[(b*i)%L0]])
			L.append([L3[i],L4[(d*i)%L0]])
		return R.lagrange_polynomial(L)
'''

'''
def creation_poly_interpolation_bis(L1,L2,b,R,l):
	
	INPUT:
	-L1 a list of abscissas of multiples of point corresponding to it's abscissa in the first entry 
	-L2 a list of abscissas of multiples of point corresponding to it's abscissa in the first entry 
	-R the polynomial ring in which we work

	OUPTUT:
	An interpolation polynomial that sends the elements of the first list L1 on elements of the second list according to the index
 	
	L=[]
	L0=len(L1)
	for i in range(1,L0):
		L.append([L1[i],L2[(b*i)%L0]])
	return interpolation_global(L,R)	

def initialisation_P(L,R):
	
	Initialisation of ploynomials Pij and the list V

	INPUT:
	-L a list of points and images of length a power of 2
	-R the Polynomial Ring we are working on

	OUTPUT:
	The polynomials Pij and the list of coefficients V 
	
	n=len(L)
	k=0
	while 2**k<n :
		k+=1
	P=[[]]
	x=R.gen()
	V=[]
	for j in range(n):
		P[0].append(x-L[j][0])
	u=n
	kb=k
	while (u>2):
		P.append([])
		for j in range(0,2**(kb-1),2):
			P[-1].append(P[-2][j]*P[-2][j+1])
		j=2**(kb-1)		
		if (u-j)>j/2:
			while (j+1)<u:
				P[-1].append(P[-2][j]*P[-2][j+1])
				j+=2
			if j<u:
				P[-1].append(P[-2][j])
		elif kb>0 :
			while j<u and u>2:
				P[-1].append(P[-2][j])
				j+=1
		u=len(P[-1])
		kb-=1
	if u==2:
		P.append([])
		P[-1].append(P[-2][0]*P[-2][1])
		
	for i in range(n):
		V.append((P[-1][0].diff()(L[i][0]))**(-1)*L[i][1])
	return P,V,k
'''
'''
def collage(L1,L2):
	P=[]
	n1=len(L1)
	print L2,2**(k-1),n
	n2=len(L2)
	for i in range(1,n2+1):
		M=L1[-i]+L2[-i]
		P.append(M)
	for i in range(n2+1,n1+1):
		P.append(L1[-i])
	return P.reverse()	
'''
'''
def transformation_liste_bis(L,n):
	M=[]
	b=bin(n).count('1')
	for i in range(b-1):
		m=L.index(0)
		c=len(L)
		M.append(L[0:m])               
		L=L[m+1:c]
	M.append(L)
	L=M[0]
	for i in range(b-2):
		for j in range(len(M[i+1])):
			L[-i-2-j]=L[-i-2-j]+M[i+1][-j-1]
	for j in range(len(M[b-1])):
		L[-b-j+1]=L[-b-j+1]+M[b-1][-j-1]
	for i in range(b-1):
		L[-b+2+i]=L[-b+2+i]+[L[-b+1+i][-1]*L[-b+1+i][-2]]
	L.append([L[-1][-1]*L[-1][-2]])			
	return L
'''

#def creation_poly_interpolation(P,Q,o1,o2,Lambda_1,Lambda_2,k2,a,b,Pb,Qb):
	#'''
	#Fonction qui calcule le polynome d'interpolation avec a et b les coefficients d'interpolation
	#'''
	#L=creation_list_interpolation(P,Q,o1,o2,Lambda_1,Lambda_2,k2,a,b,Pb,Qb)
	#return pol_inter_creation(L) 
	#return L

def creation_list_interpolation(P,Q,o1,o2,Lambda_1,Lambda_2,k2,a,b,Pb,Qb,Tower):
	'''
	Fonction qui calcule la liste à interpoler avec a et b les coefficients d'interpolation

	Input:
	-P,Q two points that generate horizontal isogenies each one 
	associated to the eigenvalue Lambda_1 (resp. Lambda_2)
	-a and b integers such that we want \phi(P)=a*Pb and \phi(Q)=b*Qb
	-Lambda_1, Lambda_2 the 2 eigenvalues associated to P and Q for 
	the Frobenius
	-o1 and o2 the multiplicative order of Lambda_1 (resp.Lambda_2)
	modulo the order of P and Q
	-k2 is the 2 valuation of the order of P and Q

	Output:
	The list of the abscissas we want to interpolate with the abscissa
	of the desired image
	'''
	T=[]
	C=[]
	M=[0]*(o2-o1+1)
	#print 'k2',k2
	for i in range(2**k2): 
		for j in range(2**k2):
			if (i%2==1)or(j%2==1 ):
				C.append((i*P+j*Q)[0])
				#print 'i,j,i*P+j*Q',i,j,i*P+j*Q
			else:
				C.append(False)
				#C.append((i*P+j*Q)[0])
		T.append(C)

		C=[]
	for h in range(o2-o1):#on fait tous les ordres au dessus de o1, o2 y compris
		for i in range(2**k2):
			for j in range(2**(k2-h)):
				if T[i][2**h*j]!=False and j%2==1:
					i1=Tower.pushD(T[i][2**h*j],h)#we consider the elements in the lowest level tower where they are defined
					j1=Tower.pushD((a*i*Pb+b*2**h*j*Qb)[0],h)#we consider the elements in the lowest level tower where they are defined
					#print "h,i1.parent()",h,i1.parent()
					C.append([i1,j1,i,2**h*j])
					#il faut aussi enlever les opposés
					T[i][2**h*j]=False
					T[-i%2**k2][-2**h*j%2**k2]=False
					i1=i
					j1=2**h*j
					for r in range(2**(o2-h)):
						i1=(i1*Lambda_1)%2**k2
						j1=(j1*Lambda_2)%2**k2
						T[i1][j1]=False
						T[-i1%2**k2][-j1%2**k2]=False
			j=0
		M[h]=C		
		#M.append([C])#chaque colonne de la matrice correspondra a un ordre
		C=[]
	if o1==o2:#cas ou les ordres sont egaux ne devrait pas arriver en fait...
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
		for i in range(2**k2):#on fait l ordre o1<o2
			for j in range(2**(k2-o2+o1)):
				if T[i][2**(o2-o1)*j]!=False and (i%2==1):#faut pas recompter ceux qui ont deja ete compte
					i1=Tower.pushD(T[i][2**(o2-o1)*j],o2-o1)#we consider the elements in the lowest level of the tower where they are defined	
					j1=Tower.pushD((a*i*Pb+b*2**(o2-o1)*j*Qb)[0],o2-o1)
					C.append([i1,j1,i,2**(o2-o1)*j])
					#il faut aussi enlever les opposés
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
	#M.append([C])
	#on a alors fini les ordres
	#print 'len(C),o1,o2,k2,Lambda_1,Lambda_2',len(C),o1,o2,k2,Lambda_1,Lambda_2	
	#for c in C:	
		#print c[-1],c[-2]
	#print 'M[0]',len(M[0]),M[0]
	#print 'M[1]',len(M[1]),M[1]
	return	M

def modif_list_interpolation(C,a,b,Pb,Qb,Tower):
	'''
	Fonction qui modifie une liste d interpolation prexistante afin qu elle corresponde aux coefficients a et b actualisés 
	Input:
	-a,b integers coefficeints for the interpolation
	-Pb,Qb points of an elliptice curver
	-Tower the 2-adic tower under which points are defined
	-C a list of abscissas of points to interpolate
	'''	
	L=[]	
	i0=Tower.floor(C[0][0])
	i=Tower.floor(Pb[0])
	nb=(i-i0)/2
	for c in C:
		L.append([c[0],Tower.pushD((c[-2]*a*Pb+c[-1]*b*Qb)[0],nb),c[-2],c[-1]])
	return L


def initialisation_poly(L,R,Tower):
	'''
	Input:
	-L a list of abscissas without their image but with their coefficients
	-R a polynomial ring

	Output:
	-M a table of Pij
	-k the bit length of the length of L
	'''
	T=[]
	for l in L:
		n=len(l)
		#we compute the 2 adic valuation of n
		k=0
		while 2**k<n :
			k+=1
		#
		M=initialisation_recu(l,R,R.gen(),Tower)
		if 2**k!=n:
			M=transformation_liste(M,n)
		T.append([k,M])	
	return T

def initialisation_recu(L,R,x,Tower):
	'''
	Initialisation of ploynomials Pij and the list V

	INPUT:
	-L a list of points and images of length a power of 2
	-R the Polynomial Ring we are working on

	OUTPUT:
	The polynomials Pij and the list of coefficients V 
	'''
	n=len(L)
	#we compute the 2 adic valuation of n
	k=0
	while 2**k<n :
		k+=1
	if 2**k!=n :
		#L1=initialisation_recu(L[0:2**(k-1)],R)
		#L2=initialisation_recu(L[2**(k-1):n],R)
		return initialisation_recu(L[0:2**(k-1)],R,x,Tower)+[0]+initialisation_recu(L[2**(k-1):n],R,x,Tower)
		#we use this format for the parsing of the output
	else:
		P=[[]]
		#x=R.gen()
		if L[0][0].parent()==Tower._base :
			R=PolynomialRing(Tower._levels[1],'x')
		else:	
			R=PolynomialRing(L[0][0].parent(),'x')
		x=R.gen()
		V=[]
		for j in range(n):
			P[0].append(x-L[j][0])
		u=n
		kb=k
		while (u>1):
			P.append([])
			for j in range(0,u,2):
				#if P[-2][j][0].parent()!=P[-2][j+1][0].parent():
					#print 'P[-2][j]',P[-2][j]
					#print  'P[-2][j+1]',P[-2][j+1]
				P[-1].append(P[-2][j]*P[-2][j+1])
				#faire ici un push down sur le polynome si c'est possible...		
			u=u/2
		return P

def transformation_liste(L,n):
	'''
	Input:
	-L a table of Pij not compatible with the subproduct tree
	-n the number of points we want to interpolate

	Output:
	-M a table of Pij compatible with the subproduct tree
	'''
	b=bin(n).count('1')
	#we count how many times it is decomposed as power of 2
	m=L.index(0)
	#because the list is cut with 0
	c=len(L)
	M=L[0:m]             
	L=L[m+1:c]
	N=[]
	for i in range(b-2):
	#we do this loop only if n is decomposed of at least 3 power of 2
		m=L.index(0)
		c=len(L)
		N=L[0:m]		
		L=L[m+1:c]
		for j in range(len(N)):
			#M[-i-2-j]=M[-i-2-j]+N[-j-1]
			M[-i-2-j]+=N[-j-1]
	for j in range(len(L)):
		#M[-b-j+1]=M[-b-j+1]+L[-j-1]
		M[-b-j+1]+=L[-j-1]
	for i in range(b-1):
		#M[-b+2+i]=M[-b+2+i]+[M[-b+1+i][-1]*M[-b+1+i][-2]]
		M[-b+2+i]+=[M[-b+1+i][-1]*M[-b+1+i][-2]]
	M.append([M[-1][-1]*M[-1][-2]])			
	return M

def precalcul_v(M,L,R,Tower):
	'''
	Function that precomputes the division of the subproduct

	Input:
	-L the list of points and images we want to interpolate
	-R the Polynomial Ring we are working on
	-M the table Pij

	Output:
	-V a vector
	'''
	n=len(L)
	V=[]
	P=M[-1][0].diff()**(-1)
	if L[0][0].parent()==Tower._base :
		R=PolynomialRing(Tower._levels[1],'x')
	else:	
		R=PolynomialRing(L[0][0].parent(),'x')
	for i in range(n):
		N=P(L[i][0])
		V.append(R(N))
	return V


def interpolation_global(M,k,L,V,R,Tower):
	'''
	Function to interpolate the list L with a polynomial on the Polynomial Ring R
	
	INPUT:
	-L the list of points and images we want to interpolate
	-R the Polynomial Ring we are working on
	-k the bit length of L
	-M the table Pij
	-V a precomputed vector 

	OUTPUT:
	The polynomial on R  which interpolates the point according to the initial list L
	'''
	n=len(L)
	if L[0][1].parent()==Tower._base :
		R=PolynomialRing(Tower._levels[1],'x')
	else:	
		R=PolynomialRing(L[0][1].parent(),'x')
	P,V=initialisation_global(M,L,V,R,n)
	return interpolation_recursive(P,L,V,k,0,n)

def initialisation_global(M,L,V,R,n):
	'''
	Function that computes the vector V useful for the interpolation

	Input:
	-M a table of the subproduc tree of the abscissas of the points we want
	to interpolate
	-L a list of abscissas of representatives of the orbits with 
	corresponding asbscissas on the other curves with also indexes of them
	-R the polynomial ring on which we want to work
	-n the length of L[i]
	
	Output:
	-M
	-V1 list of coefficients for the interpolation 
	'''
	V1=[]
	for i in range(n):
		V1.append(R(V[i]*L[i][1]))
	return M,V1

def interpolation_recursive(P,L,V,i,j,n):
	'''
	Recursive function to compute the polynomial interpolation of the list 

	INPUT:
	-P a pre computed list of polynomials in function of L
	-L the list of points and images we want to interpolate
	-V a pre computed list of coefficients in function of L
	-i an integer equals to the 2 valuation of len(L)
	-j an integer

	OUTPUT:
	The polynomial which interpolates the point according to the initial list L
	'''
	if n<2:
		return V[0]
	else :
		k=0
		while 2**k<n :
			k+=1
		u=2**(k-1)
		return interpolation_recursive(P,L[0:u],V[0:u],i-1,2*j,u)*P[i-1][2*j+1]+interpolation_recursive(P,L[u:n],V[u:n],i-1,2*j+1,n-u)*P[i-1][2*j]

def frobenius_polynomial(P,Tower,power):
	'''
	Input:
	-P a polynomial with coefficients in Tower in the polynomial ring PR
	-Tower, the 2-adic tower on which coefficients of P are defined
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

def red_pol_basis(P,Tower):
	'''
	Input:
	-P a polynomial with coefficients in Tower
	-Tower a 2-adic tower

	Ouput:
	-P the same polynomial but with his coefficients defined one level down in the Tower

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
	if P[0].parent()==Tower._levels[1]:
		return P
	Pol=PolynomialRing(Tower.push1l(P[0]).parent(),'x')
	L2=[]
	for l in L:
		L2.append(Tower.push1l(l))
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
	if r0.degree()>0: #we have to calculate the 3 prime modulus		
		TB1=TB.quo_rem(r0)[0]
		TA1=TA.quo_rem(r0)[0]
		R0=r0
		Test=True
		while Test!=False:
		#we continue the reduction until TB1 and R0 are prime
			r0=TB1; u0=1; v0=0; r1=R0; u1=0; v1=1;
			while (r1!=0):	
				q,r2=r0.quo_rem(r1)
				(r0,u0,v0,r1,u1,v1)=(r1,u1,v1,r2,u0-q*u1,v0-q*v1)
			s0=1
			if r0.degree()>0:		
				while TB1.quo_rem(r0)[1]==0:
					s0*=r0
					TB1=TB1.quo_rem(r0)[0]
				R0*=s0

			else:	
				TB1=TB1.quo_rem(r0)[0]
				R0*=r0
				Test=False
		while Test!=True:
		#we continue the reduction until TA1 and R0 are prime
			r0=TA1; u0=1; v0=0; r1=R0; u1=0; v1=1;
			while (r1!=0):	
				q,r2=r0.quo_rem(r1)
				(r0,u0,v0,r1,u1,v1)=(r1,u1,v1,r2,u0-q*u1,v0-q*v1)
			if r0.degree()>0:
				s0=1
				while TA1.quo_rem(r0)[1]==0:
					s0*=r0
					TA1=TA1.quo_rem(r0)[0]
				R0*=s0	
			else:
				TA1=TA1.quo_rem(r0)[0]
				R0*=r0
				Test=True
		#Now we have the three prime modulus TA1,TB1 and R0	 
	else:#we already have modulus prime with each other
		u0=u0.quo_rem(r0)[0]; v0=v0.quo_rem(r0)[0];
		S1=TA*TB
		R1=(TA*u0*B+v0*TB*A).quo_rem(S1)[1]
		
		return (R1,S1,TA*u0,TB*v0)
		
#2nd step compute the bezout identity
	#Test de primalite
	r0=TA1; u0=1; v0=0; r1=TB1; u1=0; v1=1;
	while (r1!=0):
		q,r2=r0.quo_rem(r1)
		(r0,u0,v0,r1,u1,v1)=(r1,u1,v1,r2,u0-q*u1,v0-q*v1)
	r0=TB1; u0=1; v0=0; r1=R0; u1=0; v1=1;
	while (r1!=0):
		q,r2=r0.quo_rem(r1)
		(r0,u0,v0,r1,u1,v1)=(r1,u1,v1,r2,u0-q*u1,v0-q*v1)
	r0=TA1; u0=1; v0=0; r1=R0; u1=0; v1=1;
	while (r1!=0):
		q,r2=r0.quo_rem(r1)
		(r0,u0,v0,r1,u1,v1)=(r1,u1,v1,r2,u0-q*u1,v0-q*v1)
	A1=A.quo_rem(TA1)[1] #useless en pratique
	B1=B.quo_rem(TB1)[1] #useless en pratique
	A2=A.quo_rem(R0)[1]
	B2=B.quo_rem(R0)[1]
	r0=TA1; u0=1; v0=0; r1=R0; u1=0; v1=1;
	while (r1!=0):
		q,r2=r0.quo_rem(r1)
		(r0,u0,v0,r1,u1,v1)=(r1,u1,v1,r2,u0-q*u1,v0-q*v1)
	S1=TA1*R0
	u0=u0.quo_rem(r0)[0]; v0=v0.quo_rem(r0)[0];
	R1=(TA1*u0*A2+v0*R0*A1).quo_rem(S1)[1]
	r0=TB1; u0=1; v0=0; r1=S1; u1=0; v1=1;
	while (r1!=0):
		q,r2=r0.quo_rem(r1)
		(r0,u0,v0,r1,u1,v1)=(r1,u1,v1,r2,u0-q*u1,v0-q*v1)
	S2=TB1*S1
	u0=u0.quo_rem(r0)[0]; v0=v0.quo_rem(r0)[0];
	R2=(TB1*u0*R1+v0*S1*B1).quo_rem(S2)[1]
	return R2,S2

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

def calcul_isogenie(P1,Q1,P2,Q2,R,l,order,T,d,Lambda_1,Lambda_2,Tower,interpol=None):
	'''
	INPUT:
	-P1 a primitive l^order-torsion point generating an horizontal l^order-isogeny on the crater
	-P2 a primitive l^order-torsion point generating an horizontal l^order-isogeny on the crater
	-Q1 a primitive l^order-torsion point generating an horizontal l^order-isogeny on the crater
	-Q2 a primitive l^order-torsion point generating an horizontal l^order-isogeny on the crater
	-l an integer
	-order an integer
	-T the l^order division polynomial on the domain curve by the isogeny
	-d the degree of the isogeny		
	
	OUTPUT:
	An isogeny between the curve which P1 belogns and the curve which P2 belongs
	'''
	i=1
	j=1
	Test=False
	o1=valuation(mod(Lambda_1,l**order).multiplicative_order(),l)
	o2=valuation(mod(Lambda_2,l**order).multiplicative_order(),l)
	power=Tower._base.cardinality()
	if o1>o2:
		(o2,o1,Lambda_2,Lambda_1,Q1,P1,Q2,P2)=(o1,o2,Lambda_1,Lambda_2,P1,Q1,P2,Q2)
	L=creation_list_interpolation(P1,Q1,o1,o2,Lambda_1,Lambda_2,order,i,j,P2,Q2,Tower)
	M=initialisation_poly(L,R,Tower)#M ne dependent pas des images choisies il est calcule pour tous les ordres du frobenius
	TA=M[0][1][-1][0]
	#1ere etape pour l ordre maximal
	L0=L[0]#on ne considere que les points lies a l ordre maximal
	#print "juste avant le precalcul"
	V=precalcul_v(M[0][1],L0,R,Tower)
	Vr=[]
	Vr.append(V)
	A=interpolation_global(M[0][1],M[0][0],L0,V,R,Tower)
	TAT=M[0][1][-1][0]#a copy of TA for this try of interpolation
	B=frobenius_polynomial(A,Tower,power**(2**(o2-1)))
	TB=frobenius_polynomial(TAT,Tower,power**(2**(o2-1))) 
	A,TA,U,V=CRT(A,TAT,B,TB)
	Lc=[]
	Lc.append([U,V])#on stocke les coefficients pour le CRT
	A=red_pol_basis(A,Tower)
	TA=red_pol_basis(TA,Tower)
	for r in range(o2-o1-1):#on fait les ordres intermediaires entre o1 et o2 et o1 NON compris
		Le=[]
		B=frobenius_polynomial(A,Tower,power**(2**(o2-1-r-1)))
		TB=frobenius_polynomial(TA,Tower,power**(2**(o2-1-r-1)))
		A,TA,U,V=CRT(A,TA,B,TB)
		A=red_pol_basis(A,Tower)
		TA=red_pol_basis(TA,Tower)
		Le.append([U,V])
		L0=L[r+1]#on ne considere que les points lies a l ordre intermediaire
		V=precalcul_v(M[r+1][1],L0,R,Tower)
		Vr.append(V)
		Aj=interpolation_global(M[r+1][1],M[r+1][0],L0,V,R,Tower)#on cree le poly interpol associe
		B=frobenius_polynomial(Aj,Tower,power**(2**(o2-1-r-1))) #on calcule son conjugue
		TAT=M[r+1][1][-1][0] # le modulus associe au poly ajoute
		TB=frobenius_polynomial(TAT,Tower,power**(2**(o2-1-r-1)))
		B,TB,U,V=CRT(Aj,TAT,B,TB)#pour le polynome ajoute
		B=red_pol_basis(B,Tower)
		TB=red_pol_basis(TB,Tower)
		Le.append([U,V])
		A,TA,U,V=CRT(A,TA,B,TB)		
		Le.append([U,V])
		Lc.append(Le) # ou Lc.append([Le])
	#ordre o1
	if o1>0: #on refait la boucle précédente normalement 
		#print "on commence la boucle o1>0"
		Le=[]
		B=frobenius_polynomial(A,Tower,power**(2**(o1-1)))
		TB=frobenius_polynomial(TA,Tower,power**(2**(o1-1)))
		A,TA,U,V=CRT(A,TA,B,TB)
		A=red_pol_basis(A,Tower)
		TA=red_pol_basis(TA,Tower)
		Le.append([U,V])
		L0=L[o2-o1]#on ne considere que les points lies a l ordre intermediaire
		V=precalcul_v(M[o2-o1][1],L0,R,Tower)
		Vr.append(V)
		Aj=interpolation_global(M[o2-o1][1],M[o2-o1][0],L0,V,R,Tower)#on cree le poly interpol associe	
		B=frobenius_polynomial(Aj,Tower,power**(2**(o1-1))) #on calcule son conjugue
		TAT=M[o2-o1][1][-1][0] # le modulus associe au poly ajoute
		TB=frobenius_polynomial(TAT,Tower,power**(2**(o1-1)))
		B,TB,U,V=CRT(Aj,TAT,B,TB)#pour le polynome ajoute
		B=red_pol_basis(B,Tower)
		TB=red_pol_basis(TB,Tower)
		Le.append([U,V])
		A,TA,U,V=CRT(A,TA,B,TB)	
		Le.append([U,V])
		Lc.append(Le)
	else :#du coup on ne fait pas agir le Frobenius...
		Le=[]		
		L0=L[o2-o1]#on ne considere que les points lies a l ordre intermediaire
		V=precalcul_v(M[o2-o1][1],L0,R,Tower)
		Vr.append(V)
		B=interpolation_global(M[o2-o1][1],M[o2-o1][0],L0,V,R,Tower)#on cree le poly interpol associe
		TB=M[o2-o1][1][-1][0] # le modulus associe au poly ajoute
		A,TA,U,V=CRT(A,TA,B,TB)		
		Lc.append([U,V])			
	for r in range(o1-1):
		B=frobenius_polynomial(A,Tower,power**(2**(o1-2-r)))
		TB=frobenius_polynomial(TA,Tower,power**(2**(o1-2-r)))
		A,TA,U,V=CRT(A,TA,B,TB)
		A=red_pol_basis(A,Tower)
		TA=red_pol_basis(TA,Tower)
		Lc.append([U,V])
	R2=PolynomialRing(A[0].parent(),R.gen())#we redefine the polynomial ring to make it match with the field we are working on
	Test=fonction_test_iso(A,TA,R2,d,Tower)
	i+=1
	while (Test==False and j<l**order):
		while(Test==False and i<l**order and j%l!=0):
			if i%l==0:
				i+=1
				#computation with pre computations
				#1ere etape pour l ordre maximal
			L0=modif_list_interpolation(L[0],i,j,P2,Q2,Tower)#on ne considere que les points lies a l ordre maximal
			A=interpolation_global(M[0][1],M[0][0],L0,Vr[0],R,Tower)
			TAT=M[0][1][-1][0]#a copy of TA for this try of interpolation
			B=frobenius_polynomial(A,Tower,power**(2**(o2-1)))
			TB=frobenius_polynomial(TAT,Tower,power**(2**(o2-1)))
			A,TA=CRTm(A,TAT,B,TB,Lc[0][0],Lc[0][1]) # on fait le CRT avec les valeurs pre calculees
			A=red_pol_basis(A,Tower)			
			TA=red_pol_basis(TA,Tower)
			for r in range(o2-o1-1):#on fait les ordres intermediaires entre o1 et o2
				Le=[]
				B=frobenius_polynomial(A,Tower,power**(2**(o2-1-r-1)))
				TB=frobenius_polynomial(TA,Tower,power**(2**(o2-1-r-1)))				
				A,TA=CRTm(A,TA,B,TB,Lc[r+1][0][0],Lc[r+1][0][1])
				A=red_pol_basis(A,Tower)
				TA=red_pol_basis(TA,Tower)
				L0=modif_list_interpolation(L[r+1],i,j,P2,Q2,Tower)#on ne considere que les points lies a l ordre intermediaire
				Aj=interpolation_global(M[r+1][1],M[r+1][0],L0,Vr[r+1],R,Tower)#on cree le poly interpol associe
				B=frobenius_polynomial(Aj,Tower,power**(2**(o2-1-r-1))) #on calcule son conjugue
				TAT=M[r+1][1][-1][0] # le modulus associe au poly ajoute
				TB=frobenius_polynomial(TAT,Tower,power**(2**(o2-1-r-1)))				
				B,TB=CRTm(Aj,TAT,B,TB,Lc[r+1][1][0],Lc[r+1][1][1])#pour le polynome ajoute
				B=red_pol_basis(B,Tower)
				TB=red_pol_basis(TB,Tower)				
				A,TA=CRTm(A,TA,B,TB,Lc[r+1][2][0],Lc[r+1][2][1])
			if o1>0: #on refait la boucle précédente normalement dans ce cas 
				B=frobenius_polynomial(A,Tower,power**(2**(o1-1)))
				TB=frobenius_polynomial(TA,Tower,power**(2**(o1-1)))				
				A,TA=CRTm(A,TA,B,TB,Lc[o2-o1][0][0],Lc[o2-o1][0][1])
				A=red_pol_basis(A,Tower)
				TA=red_pol_basis(TA,Tower)
				L0=modif_list_interpolation(L[o1-o2],i,j,P2,Q2,Tower)#on ne considere que les points lies a l ordre intermediaire			
				Aj=interpolation_global(M[o2-o1][1],M[o2-o1][0],L0,Vr[o2-o1],R,Tower)#on cree le poly interpol associe
				B=frobenius_polynomial(Aj,Tower,power**(2**(o1-1))) #on calcule son conjugue
				TAT=M[o2-o1][1][-1][0] # le modulus associe au poly ajoute
				TB=frobenius_polynomial(TAT,Tower,power**(2**(o1-1)))				
				B,TB=CRTm(Aj,TAT,B,TB,Lc[o2-o1][1][0],Lc[o2-o1][1][1])#pour le polynome ajoute
				B=red_pol_basis(B,Tower)
				TB=red_pol_basis(TB,Tower)				
				A,TA=CRTm(A,TA,B,TB,Lc[o2-o1][2][0],Lc[o2-o1][2][1])
			else :#du coup on ne fait pas agir le Frobenius... CRT A FINIR !!!
				L0=modif_list_interpolation(L[o2-o1],i,j,P2,Q2,Tower)#on ne considere que les points lies a l ordre intermediaire				
				B=interpolation_global(M[o2-o1][1],M[o2-o1][0],L0,Vr[o2-o1],R,Tower)#on cree le poly interpol associe
				TB=M[o2-o1][1][-1][0] # le modulus associe au poly ajoute
				A,TA=CRTm(A,TA,B,TB,Lc[-1][0],Lc[-1][1])		
			for r in range(o1-1):
				B=frobenius_polynomial(A,Tower,power**(2**(o1-2-r)))
				TB=frobenius_polynomial(TA,Tower,power**(2**(o1-2-r)))				
				A,TA=CRTm(A,TA,B,TB,Lc[o2-o1+r+1][0],Lc[o2-o1+r+1][1])
				A=red_pol_basis(A,Tower)
				TA=red_pol_basis(TA,Tower)		
			Test=fonction_test_iso(A,TA,R2,d,Tower)		
			i+=1
		i=1
		j+=1
	if Test!=False:
		r=Test[1]
		Test=Test[0]
		phi=Test/Test.leading_coefficient()
		return r,phi,r/(phi**2)
		#le reste du code est pour les isogenies faites a base de sous groupe de la courbe
		K=Tower._base	
		PR=PolynomialRing(Tower._base,'x')
		E22=EllipticCurve([K(P1.curve().a4()[0]),K(P1.curve().a6()[0])]) #oblige de prendre une construction sur une courbe situee sur le plus bas niveau afin de pouvoir utiliser la fonction isogeny par la suite
		L=[]
		for l in phi.list():
			L.append(l[0])
		print 'L[0],L[0].parent()',L[0],L[0].parent()
		phi2=E22.isogeny(PR([phi.list()[0][0],1]),degree=d)
		print 'phi2',phi2	
		Test=phi2
	return 	Test	

def test_square(P,Tower):
		'''
		Input:
		-P a polynomial defined over a 2-adic tower
		-Tower a 2-adic tower

		Output:
		A boolean value saying if P is a square
		'''
		DP=P.diff()
		P2=P.gcd(DP)
		P3=P.quo_rem(P2**2)
		if P3[0].degree()==0 and P3[1]==0 and P2.degree()==(P.degree()/2):
			#print 'P3[0],P3[0][0]',P3[0],P3[0][0]
			return [True,Tower.root_computing(P3[0][0])*P2,P3[0][0]]
		else:
			return [False]
		#calcul du pgcd de P et DP pour savoir si c'est un carre

def fonction_test_iso(A,T,R,d,Tower):
		'''
		INPUT
		-A a polynomial which is the interpolation of the primitive torsion points on the others primitive torsion points of the other curve
		-T the division polynomial of the torsion points
		-l the degree of the isogeny
		-F the Fraction Field of the Polynomial Ring where A and T are defined
		-d the degree of the isogeny
		-Tower a 2-adic tower under which the polynomials are defined
		
		OUTPUT
		The "denominator" of the l-isogeny that respects A
		'''
		q0,r0=T,A
		#q1,r1=divmod(T,A)
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
			#q1,r1=divmod(r1,r0)
			(q1,r1)=r0.quo_rem(r1)
			#j=r0
			#print "r1,r2",r1,r2
			r0=i
			deg=R(v0).degree()
			#on peut eventuellement retourner r0
		if (deg==d-1):
			test=test_square(v0,Tower)#v0.is_square(True) not implemented for all the polynomials
			if test[0]==True:			
				#I=F(j/v0)
				#phi=test[1]/test[1].leading_coefficient()
				#phi=EllipticCurveIsogeny(E,phi.list())
				r2=r2/test[2]
				#print "test[1],r2",test[1],r2				
				return test[1],r2 
			#tester si r0 est un carre sage le fait betement par une factorisation
			#computation of I the possible isogeny
			else :
				return test[0]
		else :
			return test

def Couveignes_algorithme(E1,E2,r,Tower):
		'''
		INPUT:
		-E1 an elliptic curve r-isogenous to E2
		-E2 an elliptic curve r-isogenous to E1
		-l an integer

		Output:
		The isogeny of degree r that goes from E1 to E2

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
		K=E1.base_field()
		#B=sqrt(r).round()
		#if B**2<r:
			#B+=1
		B=5.4*(r)
		if E2.base_field()!=K:
			raise TypeError('the curves must be defined on the same field')
		else :
			E1,P1,Q1,k1,Lambda_1,Lambda_2,Tower=tate_module(E1,B,Tower,2,conservation=True)
			E2,P2,Q2,k2,Lambda_12,Lambda_22=tate_module(E2,B,Tower,2)
			if Lambda_1==Lambda_22 and Lambda_2==Lambda_12:
				P2,Q2=Q2,P2	
			elif Lambda_1 != Lambda_12 or Lambda_2 != Lambda_22 :
				print "probleme valeur propres,Lambda_1, Lambda_12, Lambda_2, Lambda_22" ,Lambda_1, Lambda_12, Lambda_2, Lambda_22 
			R=PolynomialRing(P1[0].parent(),name='x')
			#if (Lambda_1*P1)[0]!=Tower.frobenius_computation(P1[0],101) or (Lambda_1*P2)[0]!=Tower.frobenius_computation(P2[0],101) or (Lambda_2*Q1)[0]!=Tower.frobenius_computation(Q1[0],101) or (Lambda_2*Q2)[0]!=Tower.frobenius_computation(Q2[0],101) :
				#print (Lambda_1*P1)[0]!=Tower.frobenius_computation(P1[0],101), (Lambda_1*P2)[0]!=Tower.frobenius_computation(P2[0],101), (Lambda_2*Q1)[0]!=Tower.frobenius_computation(Q1[0],101), (Lambda_2*Q2)[0]!=Tower.frobenius_computation(Q2[0],101) 
			T=E1.division_polynomial(2**k1,R.gen())
			return calcul_isogenie(P1,Q1,P2,Q2,R,2,k1,T,r,Lambda_1,Lambda_2,Tower)

def calcul_isogenie_initialisation(P1,Q1,P2,Q2,R,l,order,T,d,Lambda_1,Lambda_2,Tower,interpol=None):
	'''
	INPUT:
	-P1 a primitive l^order-torsion point generating an horizontal l^order-isogeny on the crater
	-P2 a primitive l^order-torsion point generating an horizontal l^order-isogeny on the crater
	-Q1 a primitive l^order-torsion point generating an horizontal l^order-isogeny on the crater
	-Q2 a primitive l^order-torsion point generating an horizontal l^order-isogeny on the crater
	-l an integer
	-order an integer
	-T the l^order division polynomial on the domain curve by the isogeny
	-d the degree of the isogeny		
	
	OUTPUT:
	An isogeny between the curve which P1 belogns and the curve which P2 belongs
	'''
	i=1
	j=1
	#M1=creation_matrice(P1,Q1,l,order,R)
	#M2=creation_matrice(P2,Q2,l,order,R)
	Test=False
	#FractionField(R)
	#L=creation_list_interpolation(M1,M2,i,j,R)
	o1=valuation(mod(Lambda_1,l**order).multiplicative_order(),l)
	o2=valuation(mod(Lambda_2,l**order).multiplicative_order(),l)
	power=Tower._base.cardinality()
	if o1>o2:
		(o2,o1,Lambda_2,Lambda_1,Q1,P1,Q2,P2)=(o1,o2,Lambda_1,Lambda_2,P1,Q1,P2,Q2)
	L=creation_list_interpolation(P1,Q1,o1,o2,Lambda_1,Lambda_2,order,i,j,P2,Q2,Tower)
	#print 'L',L
	#print 'len(L[0]),len(L[1]),len(L[2]),L',len(L[0]),len(L[1]),L[0][0],L[1][0]
	M=initialisation_poly(L,R,Tower)#M ne dependent pas des images choisies il est calcule pour tous les ordres du frobenius
	#print 'len(L),order,o1,o2,Lambda_1,Lambda_2,power',len(L),order,o1,o2,Lambda_1,Lambda_2,power
	TA=M[0][1][-1][0]
	#for r in range(len(L[0])):
		#if TA(L[0][r][0])!=0:
			#print 'TA(L[r][0])',TA(L[r][0])
	#1ere etape pour l ordre maximal
	L0=L[0]#on ne considere que les points lies a l ordre maximal
	#print "juste avant le precalcul"
	V=precalcul_v(M[0][1],L0,R,Tower)
	Vr=[]
	Vr.append(V)
	#print "juste avant l interpolation globale"
	A=interpolation_global(M[0][1],M[0][0],L0,V,R,Tower)
	#for r in range(len(L0)):
		#if A(L0[r][0])!=L0[r][1]:
			#print 'A(L[0][r][0]),L[i]',r,A(L0[r][0]),L0[r]
	TAT=M[0][1][-1][0]#a copy of TA for this try of interpolation
	B=frobenius_polynomial(A,Tower,power**(2**(o2-1)))
	TB=frobenius_polynomial(TAT,Tower,power**(2**(o2-1))) 
	#print "juste avant le CRT"
	#print "B,TB",B[0].parent(),TB[0].parent()
	A,TA,U,V=CRT(A,TAT,B,TB)
	Lc=[]
	Lc.append([U,V])#on stocke les coefficients pour le CRT
	#print "fin de la premiere etape"
	#for r in range(len(L0)):
		#if A(L0[r][0])!=L0[r][1]:
			#print 'A(L[0][i][0]),L[i]',r,A(L0[r][0]),L0[r], 'test normalement deja verifie'
	#print "test du polynome d interpolation passe"
	A=red_pol_basis(A,Tower)
	TA=red_pol_basis(TA,Tower)
	for r in range(o2-o1-1):#on fait les ordres intermediaires entre o1 et o2 et o1 NON compris
		Le=[]
		#print "r,o2,o1,o2-1-r-1",r,o2,o1,o2-1-r-1,o2>o1
		B=frobenius_polynomial(A,Tower,power**(2**(o2-1-r-1)))
		TB=frobenius_polynomial(TA,Tower,power**(2**(o2-1-r-1)))
		A,TA,U,V=CRT(A,TA,B,TB)
		A=red_pol_basis(A,Tower)
		TA=red_pol_basis(TA,Tower)
		Le.append([U,V])
		L0=L[r+1]#on ne considere que les points lies a l ordre intermediaire
		#print 'L[r+1]',L[r+1],'M[r+1][1]',M[r+1][1]
		V=precalcul_v(M[r+1][1],L0,R,Tower)
		#print "deuxieme precalcul passe,r",r
		Vr.append(V)
		Aj=interpolation_global(M[r+1][1],M[r+1][0],L0,V,R,Tower)#on cree le poly interpol associe
		#for s in range(len(L0)):
			#if Aj(L0[s][0])!=L0[s][1]:
				#print 'A(L[i][0]),L[i]',s,Aj(L0[s][0]),L0[s]
		B=frobenius_polynomial(Aj,Tower,power**(2**(o2-1-r-1))) #on calcule son conjugue
		TAT=M[r+1][1][-1][0] # le modulus associe au poly ajoute
		TB=frobenius_polynomial(TAT,Tower,power**(2**(o2-1-r-1)))
		B,TB,U,V=CRT(Aj,TAT,B,TB)#pour le polynome ajoute
		B=red_pol_basis(B,Tower)
		TB=red_pol_basis(TB,Tower)
		Le.append([U,V])
		A,TA,U,V=CRT(A,TA,B,TB)		
		#for s in range(len(L[0])):
			#if A[0].parent()!=
			#if A(L[0][s][0])!=L[0][s][1]:
				#print 'A(L[i][0]),L[i]',s,A(L[0][s][0]),L[0][s]
		#for s in range(len(L[1])):
			#if A(L[1][s][0])!=L[1][s][1]:
				#print 'A(L[i][0]),L[i]',s,A(L[1][s][0]),L[1][s]
		Le.append([U,V])
		Lc.append(Le) # ou Lc.append([Le])
	#ordre o1
	if o1>0: #on refait la boucle précédente normalement dans ce cas A FINIR !!!
		#print "on commence la boucle o1>0"
		Le=[]
		B=frobenius_polynomial(A,Tower,power**(2**(o1-1)))
		TB=frobenius_polynomial(TA,Tower,power**(2**(o1-1)))
		A,TA,U,V=CRT(A,TA,B,TB)
		A=red_pol_basis(A,Tower)
		TA=red_pol_basis(TA,Tower)
		Le.append([U,V])
		L0=L[o2-o1]#on ne considere que les points lies a l ordre intermediaire
		#print 'L[o2-o1]',L[o2-o1],'M[o2-o1][1]',M[o2-o1][1],'M[o2-o1][1][-1][0]',M[o2-o1][1][-1][0]#, M[o2-o1][1][-1][0].factor()		
		V=precalcul_v(M[o2-o1][1],L0,R,Tower)
		Vr.append(V)
		Aj=interpolation_global(M[o2-o1][1],M[o2-o1][0],L0,V,R,Tower)#on cree le poly interpol associe	
		#for s in range(len(L0)):
			#if Aj(L0[s][0])!=L0[s][1]:
				#print 'A(L[i][0]),L[i]',s,Aj(L0[s][0]),L0[s]
		B=frobenius_polynomial(Aj,Tower,power**(2**(o1-1))) #on calcule son conjugue
		TAT=M[o2-o1][1][-1][0] # le modulus associe au poly ajoute
		TB=frobenius_polynomial(TAT,Tower,power**(2**(o1-1)))
		B,TB,U,V=CRT(Aj,TAT,B,TB)#pour le polynome ajoute
		B=red_pol_basis(B,Tower)
		TB=red_pol_basis(TB,Tower)
		Le.append([U,V])
		A,TA,U,V=CRT(A,TA,B,TB)	
		#print "A[0].parent()",A[0].parent()
		#print "L[0][0][0].parent()",L[0][0][0].parent()	
		#for s in range(len(L[0])):
			#if A(L[0][s][0])!=L[0][s][1]:
				#print 'A(L[i][0]),L[i]',s,A(L[0][s][0]),L[0][s]
		#for s in range(len(L[1])):
			#if A(L[1][s][0])!=L[1][s][1]:
				#print 'A(L[i][0]),L[i]',s,A(L[1][s][0]),L[1][s]
		Le.append([U,V])
		Lc.append(Le)
	else :#du coup on ne fait pas agir le Frobenius...
		Le=[]		
		#print "L[o2-o1]",L[o2-o1]		
		#B=frobenius_polynomial(A,Tower,power**(2**(o1-1)))
		#TB=frobenius_polynomial(TA,Tower,power**(2**(o1-1)))
		#A,TA,U,V=CRT(A,TA,B,TB)
		#A=red_pol_basis(A,Tower)
		#TA=red_pol_basis(TA,Tower)
		#Le.append([U,V])
		L0=L[o2-o1]#on ne considere que les points lies a l ordre intermediaire
		#for s in range(len(L[o2-o1])):
			#if A(L[o2-o1][s][0])!=L[o2-o1][s][1]:
				#print 's,A(L[o2-o1][s][0]),L[o2-o1][s]',s,A(L[o2-o1][s][0]),L[o2-o1][s],'o1=0','avant CRT'
		V=precalcul_v(M[o2-o1][1],L0,R,Tower)
		Vr.append(V)
		B=interpolation_global(M[o2-o1][1],M[o2-o1][0],L0,V,R,Tower)#on cree le poly interpol associe
		#for s in range(len(L0)):
			#if B(L0[s][0])!=L0[s][1]:
				#print 'B(L[o2-o1][0]),L[i]',s,B(L0[s][0]),L0[s]
		#B=frobenius_polynomial(Aj,Tower,power**(2**(o1-1))) #on calcule son conjugue
		TB=M[o2-o1][1][-1][0] # le modulus associe au poly ajoute
		#TB=frobenius_polynomial(TAT,Tower,power**(2**(o1-1)))
		#B,TB,U,V=CRT(Aj,TAT,B,TB)#pour le polynome ajoute
		#B=red_pol_basis(Aj,Tower)
		#TB=red_pol_basis(TAT,Tower)
		#Le.append([U,V])
		#print 'A.degree(),TA.degree(),B.degree(),TB.degree()',A.degree(),TA.degree(),B.degree(),TB.degree()
		#for s in range(len(L[o2-o1-1])):
			#if A(L[o2-o1-1][s][0])!=L[o2-o1-1][s][1]:
				#print 'A(L[o2-o1-1][0]),L[i]',s,A(L[o2-o1-1][s][0]),L[o2-o1-1][s],'o1=0','premier essai avant CRT'
		A,TA,U,V=CRT(A,TA,B,TB)		
		#for s in range(len(L[o2-o1])):
			#if A(L[o2-o1][s][0])!=L[o2-o1][s][1]:
				#print 'A(L[o2-o1][0]),L[i]',s,A(L[o2-o1][s][0]),L[0][s],'o1=0','deuxieme essai apres CRT'
		#for s in range(len(L[o2-o1-1])):
			#if A(L[o2-o1-1][s][0])!=L[o2-o1-1][s][1]:
				#print 'A(L[o2-o1-1][0]),L[o2-o1-1]',s,A(L[o2-o1-1][s][0]),L[o2-o1-1][s],'o1=0','deuxieme essai different apres CRT'
		Lc.append([U,V])				
	for r in range(o1-1):
		B=frobenius_polynomial(A,Tower,power**(2**(o1-2-r)))
		TB=frobenius_polynomial(TA,Tower,power**(2**(o1-2-r)))
		A,TA,U,V=CRT(A,TA,B,TB)
		A=red_pol_basis(A,Tower)
		TA=red_pol_basis(TA,Tower)
		Lc.append([U,V])
	return P2,Q2,Tower,L,M,Lc,Vr,o2,o1

def calcul_isogenie_etape(P2,Q2,Tower,L,M,Lc,Vr,o2,o1,l):
	i=3
	j=3
	Test=False
	R=PolynomialRing(P2[0].parent(),name='x')
	power=Tower._base.cardinality()	
	while (Test==False and j<16):
		while(Test==False and i<16 and j%l!=0):
			if i%l==0:
				i+=1
			L0=modif_list_interpolation(L[0],i,j,P2,Q2,Tower)#on ne considere que les points lies a l ordre maximal
			A=interpolation_global(M[0][1],M[0][0],L0,Vr[0],R,Tower)
			TAT=M[0][1][-1][0]#a copy of TA for this try of interpolation
			B=frobenius_polynomial(A,Tower,power**(2**(o2-1)))
			TB=frobenius_polynomial(TAT,Tower,power**(2**(o2-1)))
			A,TA=CRTm(A,TAT,B,TB,Lc[0][0],Lc[0][1]) # on fait le CRT avec les valeurs pre calculees
			A=red_pol_basis(A,Tower)			
			TA=red_pol_basis(TA,Tower)
			for r in range(o2-o1-1):#on fait les ordres intermediaires entre o1 et o2
				Le=[]
				B=frobenius_polynomial(A,Tower,power**(2**(o2-1-r-1)))
				TB=frobenius_polynomial(TA,Tower,power**(2**(o2-1-r-1)))				
				A,TA=CRTm(A,TA,B,TB,Lc[r+1][0][0],Lc[r+1][0][1])
				A=red_pol_basis(A,Tower)
				TA=red_pol_basis(TA,Tower)
				L0=modif_list_interpolation(L[r+1],i,j,P2,Q2,Tower)#on ne considere que les points lies a l ordre intermediaire
				Aj=interpolation_global(M[r+1][1],M[r+1][0],L0,Vr[r+1],R,Tower)#on cree le poly interpol associe
				B=frobenius_polynomial(Aj,Tower,power**(2**(o2-1-r-1))) #on calcule son conjugue
				TAT=M[r+1][1][-1][0] # le modulus associe au poly ajoute
				TB=frobenius_polynomial(TAT,Tower,power**(2**(o2-1-r-1)))			
				B,TB=CRTm(Aj,TAT,B,TB,Lc[r+1][1][0],Lc[r+1][1][1])#pour le polynome ajoute
				B=red_pol_basis(B,Tower)
				TB=red_pol_basis(TB,Tower)			
				A,TA=CRTm(A,TA,B,TB,Lc[r+1][2][0],Lc[r+1][2][1])
			if o1>0: #on refait la boucle précédente normalement dans ce cas 
				#Le=[]
				B=frobenius_polynomial(A,Tower,power**(2**(o1-1)))
				TB=frobenius_polynomial(TA,Tower,power**(2**(o1-1)))							
				A,TA=CRTm(A,TA,B,TB,Lc[o2-o1][0][0],Lc[o2-o1][0][1])
				A=red_pol_basis(A,Tower)
				TA=red_pol_basis(TA,Tower)
				L0=modif_list_interpolation(L[o1-o2],i,j,P2,Q2,Tower)#on ne considere que les points lies a l ordre intermediaire			
				Aj=interpolation_global(M[o2-o1][1],M[o2-o1][0],L0,Vr[o2-o1],R,Tower)#on cree le poly interpol associe
				B=frobenius_polynomial(Aj,Tower,power**(2**(o1-1))) #on calcule son conjugue
				TAT=M[o2-o1][1][-1][0] # le modulus associe au poly ajoute
				TB=frobenius_polynomial(TAT,Tower,power**(2**(o1-1)))				
				B,TB=CRTm(Aj,TAT,B,TB,Lc[o2-o1][1][0],Lc[o2-o1][1][1])#pour le polynome ajoute
				B=red_pol_basis(B,Tower)
				TB=red_pol_basis(TB,Tower)				
				A,TA=CRTm(A,TA,B,TB,Lc[o2-o1][2][0],Lc[o2-o1][2][1])
			else :#du coup on ne fait pas agir le Frobenius... CRT A FINIR !!!
				L0=modif_list_interpolation(L[o2-o1],i,j,P2,Q2,Tower)#on ne considere que les points lies a l ordre intermediaire				
				B=interpolation_global(M[o2-o1][1],M[o2-o1][0],L0,Vr[o2-o1],R,Tower)#on cree le poly interpol associe
				TB=M[o2-o1][1][-1][0] # le modulus associe au poly ajoute
				A,TA=CRTm(A,TA,B,TB,Lc[-1][0],Lc[-1][1])		
			for r in range(o1-1):
				B=frobenius_polynomial(A,Tower,power**(2**(o1-2-r)))
				TB=frobenius_polynomial(TA,Tower,power**(2**(o1-2-r)))
				#print "9 eme occurence"				
				A,TA=CRTm(A,TA,B,TB,Lc[o2-o1+r+1][0],Lc[o2-o1+r+1][1])
				A=red_pol_basis(A,Tower)
				TA=red_pol_basis(TA,Tower)
			R2=A.parent()			
			Test=fonction_test_iso(A,TA,R2,d,Tower)
			return Test
