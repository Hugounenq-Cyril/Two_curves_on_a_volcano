def Calcul_j_invariant_cratere_isogenes(p,r,j4,p1):
	'''
	Input:
	-j4 le j invariant d'une courbe définie sur Fp^r
	-p un nombre premier
	-r un entier
	-p1 un nombre premier
	
	Output:
	Le j-invariant d une courbe situee sur le volcan et une liste et courbes p* isogenes pour p* un des 5 nombres premiers superieur a p1 
	'''
	if volcan_cyclique(p,r,j4):
		j4=Construction_Montante(p,r,j4)
		return volcan_cyclique_et_courbes_isogenes(p,r,j4,p1),j4
	else: 
		print "volcan pas cyclique"
		return False
	

def Construction_Montante(p,r,j4):
	'''
	Input:
	-j4 le j invariant d'une courbe définie sur Fp^r
	-p un nombre premier
	-r un entier

	Output:
	La courbe la plus proche sur le cratère de la courbe de j-invariant j4
	'''
	L=Construction_Descente(p,r,True,j4)
	while (L[0][-1]!=L[1][-1] and L[2][-1]!=L[1][-1] and L[0][-1]!=L[2][-1]):
		if L[0][-1]>L[2][-1] and L[0][1]>L[1][-1]:
			j4=L[0][1]
		elif L[1][-1]>L[2][-1] and L[1][1]>L[0][-1]:
			j4=L[1][1]
		else :
			j4=L[2][1]
		L=Construction_Descente(p,r,True,j4)
	if (L[0][-1]<L[1][-1] and L[0][-1]<L[2][-1]):
		return L[0][0]
	elif (L[1][-1]<L[0][-1] and L[0][-1]<L[2][-1]):
		return L[1][0]
	else:
		return L[2][0]


def volcan_cyclique_et_courbes_isogenes(p,r,j4,p1):
	'''
	Input:
	-j4 le j invariant d'une courbe définie sur Fp^r
	-p un nombre premier
	-r un entier
	-p1 un nombre premier

	Output:
	Dit si le volcan est cyclique ou pas et si la courbe admet des courbes p*-isogénes avec p* un des 5 nombre premiers à partir de p1
	'''
	if volcan_cyclique(p,r,j4):
		F.<a>=FiniteField(p^r)
		E=EllipticCurve(j=F(j4))
		BRAP=ClassicalModularPolynomialDatabase()
		L=[]
		for i in range(5):
			f=BRAP[p1]
			f=f.subs(j1=x)
			g=f.subs(j0=j4)
			Li=g.roots(ring=F)
			if len(Li)>0:
				M=[p1,Li]
				L.append(M)
			p1=next_prime(p1)
		if len(L)==0:
			print "pas de r isogenie"
			return False
		else:
			return L
	else:
		print "volcan pas cyclique"
		return False

def volcan_cyclique(p,r,j0):
	'''
	Input:
	-j0 le j invariant d'une courbe définie sur Fp^r
	-p un nombre premier
	-r un entier

	Output:
	Dit si le volcan est cyclique ou pas 
	'''
	F.<a>=FiniteField(p^r)
	E=EllipticCurve(j=F(j0))
	t=E.trace_of_frobenius()
	D=t^2-4*(p^r)
	while (D%4)==0:
		D=D/4
	if D%8==1:
		print 'D%8==',D%8
		return True
	else:
		print 'D%8==',D%8
		return False


def Construction_Descente(p,r,b,u):
	k.<a>=FiniteField(p^r)
	if b : #b est donc un booléen pour indiquer si on le fait à partir d'une courbe aléatoire ou choisie d'après son j-invariant u
		E=EllipticCurve(j=k(u)) 
	else :
		v=randint(0,p^r-1)
		E=EllipticCurve(j=a^v)   # on prend une courbe elliptique au hasard en prenant une puissance au hasard du générateur a
        	print E
	if E.is_singular()==True :
		#print "echec ordinaire"
                return E
	j4=E.j_invariant()
	BRAP=ClassicalModularPolynomialDatabase()
	f=BRAP[2]
	f=f.subs(j1=x) # pour avoir le classical modular polynome exprimé directement en l'indéterminé
	g=f.subs(j0=j4) # pour calculer les courbes 2 isogenes à partir de celle de départ 
	F=g.roots(ring=k)
	L=[[],[],[]]
	c=[0,0,0]
	if (len(F)<=2 and E.two_torsion_rank()==1):
		print "echec longueur"
                return [[j4,1],[j4,1],[j4,1]]
	if (len(F)==2 and E.two_torsion_rank()==2):
		for i in range(0,3) :
			L[i].append(j4)
			if (F[0][0]==j4):
				j4=F[1][0]
			else :
				j4=F[0][0] 
		g=f.subs(j0=j4)
		F=g.roots(ring=k)
		if (len(F)==2):
			for i in range(0,3):
				L[i].append(j4)
				L[i].append(len(L[i]))
			return L	
	t=1
	d=[0,0,0]
	j3=[0,0,0]
	A=[j4,j4,j4]# on conserve d'où l'on vient
	for i in range(0,3) :
		j3[i]=F[i][0]
                L[i].append(j4)
		L[i].append(F[i][0])
	while t :
		for i in range(0,3) : #on n'a besoin que de prendre 3 chemins au départ et pas à chaque étape
			g=f.subs(j0=j3[i])
			F=g.roots(ring=k)
			if (len(F)==1) & (F[0][0]==A[i]) :# regarde si on est pas arrivé en bas du cratère
				t=0
				d[i]=1 #on change alors cet indicateur pour dire que l'on a un chemin descendant
				
			else :
				if len(F)>1 : 
					if F[0][0]==A[i] : #regarde si on ne retourne pas sur ses pas
						L[i].append(F[1][0])
						A[i]=j3[i]
						j3[i]=F[1][0]
					else :
						L[i].append(F[0][0])
						A[i]=j3[i]
						j3[i]=F[0][0]
				else :
					d[i]='b' # on indique que l'on est bloqué et que l'on ne peut pas descendre dans ce corps sur cette pente
	for i in range(0,3):
		L[i].append(len(L[i]))
		#print L[i], d[i]
	return L # retourne la liste des 3 "descentes" , et le j-invariant de la courbe de départ (utile si on demande d'en choisir une au hasard) et aussi pour calculer les isogénies correspondantes
def Conversion(t,p,r): # fonction auxilaire qui permet de prendre au hasard un élément de FF(3^32) il faut la compiler après la création du corps
	res=0
	l=[]
	for i in range(0,r):
                s=r-1
		e=s-i
                b=p^e
		l=divmod(t,b)
		if  l[0]!=0:
			res=l[0]*a^e+res
			t=l[1]
		else :
			t=l[1]
	return res
	
