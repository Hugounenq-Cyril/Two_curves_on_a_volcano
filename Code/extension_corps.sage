from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_field


class Tower_two:

	def __init__(self,K,h,name='x',alpha=None,implementation=None):
		"""
        	Initialize a 2-adic extension tower of K of initial height h

		INPUT:
		-K a finite field
		-h an integer
		-name character to design the variable
		-implementation for the method to build the polynomial ring

		OUTPUT:
		A 2-adic tower od height h
	
		EXAMPLES:
		sage: A=FiniteField(97)
		sage: tower=Tower_two(A,12,'y')
		sage: tower._top
		Univariate Quotient Polynomial Ring in y1 over Finite Field of size 97 
		with modulus y1^4096 + 69
		sage: tower._top.cardinality().factor()
		97^4096
		sage: tower._base
		Finite Field of size 97
		sage: tower._levels
		[Finite Field of size 97,
		 Univariate Quotient Polynomial Ring in y1 over Finite Field of size 97
		 with modulus y1^4096 + 69]
		sage: tower._alpha
		28

		
		sage: q=next_prime(2^65)
		sage: towert=Tower_two(FiniteField(q),4,'y')
		sage: q%4
		3
		sage: towert._top
		Univariate Quotient Polynomial Ring in y2 over Finite Field in 
		x of size 36893488147419103363^2 with modulus y2^8 + 2268659324
		5818090076*x + 30949137583800564139
		sage: valuation(towert._top_cardinal,q)==valuation(towert._top.
		cardinality(),q)
		True
		sage: valuation(towert._top_cardinal,q)
		16	
		"""
		if not K.is_field():
           		raise RuntimeError('This must be a field.')	
		if not K.is_finite():
            		raise RuntimeError('The field must be finite.')
		self._levels=[K]
		self._name=name		
		self._P=PolynomialRing(K,self._name + str(len(self._levels)), implementation=implementation)
		self._top=K
		self._top_cardinal=K.cardinality()
		self._base=K
		q=self._top_cardinal
		self._euler_phi=euler_phi(q)		
		if q%4==3:	
			self._base_name='y'
			if (self._name=='y'):
				self._base_name='x'
			K=FiniteField(q**2,self._base_name)
			self._levels.append(K)
			self._top=K
			self._base=K
			self._P = PolynomialRing(K, self._name + str(len(self._levels)), implementation=implementation)
			h=h-1
			self._top_cardinal=self._top_cardinal**2
			self._euler_phi=euler_phi(q^2)
		q=self._top_cardinal #compute a non 2 adic residue
		test=euler_phi(q)/2
		if alpha==None:
			self._alpha=self._top.random_element()
		else:
			self._alpha=alpha
		#self._alpha=35
		while (self._alpha**test==1 or self._alpha==self._top.zero()):
			self._alpha=self._top.random_element()
		#compute a 2 adic tower of heigt h
		self._top=PolynomialQuotientRing_field(self._P,self._P(self._P.gen()**(2**h)-self._alpha) ,self._name + str(len(self._levels)))
		#on met a jour les donnees selon ou l on se trouve sur la tour
		self._levels.append(self._top)
		self._top_cardinal=self._top_cardinal**(2**h)
		
	def frobenius_computation(self,element,power):
		'''
		Compute the Frobenius on an element of a 2-adic tower

		INPUT:
		-self a 2-adic tower
		-element an element who belongs to the 2-adic tower
		-power an integer who is the power of the base-field of the 2-adic 
		tower
		
		OUTPUT:
		Compute the power-th power of element

		EXAMPLE:
		sage: A=FiniteField(66797)
		sage: rand=tower._top.random_element();rand
		47238*x1^15 + 64720*x1^14 + 39952*x1^13 + 35645*x1^12 + 
		48979*x1^11 + 35081*x1^10 + 24782*x1^9 + 262*x1^8 + 37146*x1^7 
		+ 5356*x1^6 + 24170*x1^5 + 51789*x1^4 + 23251*x1^3 + 37276*x1^2
		 + 46087*x1 + 14362
		sage: tower.frobenius_computation
		(rand,66797^13)==rand^(66797^13)
		True
		sage: time(tower.frobenius_computation(rand,66797^15))
		CPU times: user 24 ms, sys: 0 ns, total: 24 ms
		Wall time: 45.9 ms
		30273*x1^15 + 19737*x1^14 + 20814*x1^13 + 36423*x1^12 + 
		31022*x1^11 + 56956*x1^10 + 65478*x1^9 + 66535*x1^8 + 9182*x1^7
		 + 52462*x1^6 + 63167*x1^5 + 42045*x1^4 + 50293*x1^3 + 
		49814*x1^2 + 39630*x1 + 14362
		sage: time(rand^(66797^15))
		CPU times: user 28 ms, sys: 0 ns, total: 28 ms
		Wall time: 34 ms
		30273*x1^15 + 19737*x1^14 + 20814*x1^13 + 36423*x1^12 + 
		31022*x1^11 + 56956*x1^10 + 65478*x1^9 + 66535*x1^8 + 9182*x1^7
		 + 52462*x1^6 + 63167*x1^5 + 42045*x1^4 + 50293*x1^3 + 
		49814*x1^2 + 39630*x1 + 14362
		'''
		K=element.parent()
		n=self.cardinality_field(K)	
		rho= power % (n*self._euler_phi)
		u,v=divmod(rho,n)
		#c=self._alpha**u*(K.gen()**v)
		a=self._alpha**u
		#c1=1
		#retour=K(0)
		L=element.list()		
		#for i in element:
		h=valuation(self.cardinality_field(K),self._base.cardinality())#on calcule la hauteur de l element dans la tour
		if self.test_bi(element):#element est en bivarie
			length=len(L[1].list())
			L1=L[0].list()
			L2=L[1].list()
			for i in range(0,length):#on agit sur chacune des deux entrees separement a cause du cas qui nous interesse
				q,r=divmod(v*(2*i+1),h)
				L2[(r-1)/2]=self._alpha**(q)*a**((2*i+1)%self._euler_phi)*element[1][i]
			if (power!=sqrt(n)):#on s evite de faire un calcul trivial que l on sera ammene a faire
				for i in range(0,length):
					q,r=divmod(v*(2*i),h)
					L1[r/2]=self._alpha**(q)*a**((2*i)%self._euler_phi)*element[0][i]
			L=[L1,L2]
		else :
			for i in range(1,len(L)):
				q,r=divmod(v*i,h)
				L[r]=self._alpha**(q)*a**(i%self._euler_phi)*element[i]		
			#retour+=i*c1
			#c1*=c
		#return retour
		return K(L)
	
	def add_one_level(self):
		'''
		Compute an extra quadratic extension on the 2-adic tower

		INPUT:
		-self a 2-adic tower

		OUPUT:
		-self a 2-adic tower with 2 extras levels one bivariate one 
		univariate of the quadratic extension of the previous top level

		EXAMPLE:
		sage: A=FiniteField(2257117)
		sage: tower=Tower_two(A,5)
		sage: tower._top
		Univariate Quotient Polynomial Ring in x1 over Finite Field of 
		size 2257117 with modulus x1^32 + 1302647
		sage: tower.add_one_level()
		sage: print tower._levels[-1], 
		tower._levels	[-1]==tower._top
		Univariate Quotient Polynomial Ring in x3 over Finite Field of 
		size 2257117 with modulus x3^64 + 1302647 True
		sage: print tower._levels[-2]
		Univariate Quotient Polynomial Ring in x2 over Univariate 
		Quotient Polynomial Ring in x1 over Finite Field of size 
		2257117 with modulus x1^32 + 1302647 with modulus x2^2 + 
		2257116*x1
		sage: print tower._levels[-3]
		Univariate Quotient Polynomial Ring in x1 over Finite Field of 
		size 2257117 with modulus x1^32 + 1302647
		'''
		self._P = PolynomialRing(self._top, self._name + str(len(self._levels)))
		self._top=PolynomialQuotientRing_field(self._P,self._P(self._P.gen()**(2)-self._top.gen()),self._name + str(len(self._levels)))
		self._levels.append(self._top)
		self._top_cardinal=self._top_cardinal**2
		h=valuation(valuation(self._top_cardinal,self._base.cardinality()),2)
		B=PolynomialRing(self._base,self._name + str(len(self._levels)))
		self._top=PolynomialQuotientRing_field(B,B(B.gen()**(2**h)-self._alpha) ,self._name + str(len(self._levels)))
		self._levels.append(self._top)
		
	def floor(self,element,depart=None):
		'''
		Compute the index of the level of the field which the element belongs to in the 2-adic tower

		INPUT:
		-self a 2-adic tower
		-element an element who belongs to the 2-adic tower

		OUTPUT:
		An integer who is the index of the level of the field which the
		element belongs to.

		EXAMPLE:
		sage: A=FiniteField(36740236697)
		sage: tower=Tower_two(A,6)
		sage: rand=tower._top.random_element()
		sage: tower.floor(rand)
		-1
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: tower.floor(rand)
		-11
		sage: rand.parent()==tower._levels[-11]
		True
		sage: tower.floor(rand,-1)
		-11
		'''
		if depart==None:
			if element.parent()==self._levels[0]:
				depart=0 			
			elif self.test_bi(element):
				depart=-2
			else:
				depart=-1
		while (self._levels[depart]!=element.parent()):
			depart=depart-2
		return depart

	def floorfield(self,field,depart=None):
		'''
		Compute the index of the level of the field which the element belongs to in the 2-adic tower

		INPUT:
		-self a 2-adic tower
		-field a field who belongs to the 2-adic tower

		OUTPUT:
		An integer who is the index of the level of the field which the
		element belongs to.

		EXAMPLE:
		
		'''
		if depart==None:
			element=field.random_element()
			if element==self._levels[0]:
				depart=0 			
			elif self.test_bi(element):
				depart=-2
			else:
				depart=-1
		while (self._levels[depart]!=field):
			depart=depart-2
		return depart


	def bivariate_to_univariate(self,element,Test=True):
		'''
		Compute the univariate representation of an element according to it's bivariate represenation

		INPUT:
		-self a 2-adic tower
		-element an element who belongs to the 2 adic tower
		
		OUTPUT:
		-a univariate representation of element
		
		EXAMPLE:
		sage: A=FiniteField(5516801)
		sage: tower=Tower_two(A,3,'t')
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: rand=tower._levels[-2].random_element();rand
		(3368243*t3^15 + 489656*t3^14 + 5217401*t3^13 + 669200*t3^12 + 
		3688246*t3^11 + 5091882*t3^10 + 5420720*t3^9 + 4848371*t3^8 + 
		4602711*t3^7 + 1619338*t3^6 + 2051308*t3^5 + 99385*t3^4 + 
		649798*t3^3 + 3699297*t3^2 + 1979671*t3 + 1797999)*t4 + 
		2152187*t3^15 + 4329412*t3^14 + 1194781*t3^13 + 3551307*t3^12 +
		 3881994*t3^11 + 139190*t3^10 + 3811588*t3^9 + 4707251*t3^8 + 
		438552*t3^7 + 1960523*t3^6 + 3427830*t3^5 + 4859068*t3^4 + 
		820223*t3^3 + 4703466*t3^2 + 2113839*t3 + 584831
		sage: randu=tower._levels[-1].random_element();randu
		4418717*t5^31 + 3954899*t5^30 + 676252*t5^29 + 2657239*t5^28 + 
		2786632*t5^27 + 137005*t5^26 + 881136*t5^25 + 3296396*t5^24 + 
		944204*t5^23 + 4494128*t5^22 + 2307926*t5^21 + 4704121*t5^20 + 
		1253327*t5^19 + 4985720*t5^18 + 4492812*t5^17 + 2139076*t5^16 +
		 1624290*t5^15 + 4838348*t5^14 + 682868*t5^13 + 2005323*t5^12 +
		 3922200*t5^11 + 3207602*t5^10 + 1781287*t5^9 + 4805664*t5^8 + 
		533434*t5^7 + 1147554*t5^6 + 2164213*t5^5 + 486436*t5^4 + 
		2720797*t5^3 + 876219*t5^2 + 676256*t5 + 557455
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: tower.bivariate_to_univariate(rand)
		3368243*t5^31 + 2152187*t5^30 + 489656*t5^29 + 4329412*t5^28 + 
		5217401*t5^27 + 1194781*t5^26 + 669200*t5^25 + 3551307*t5^24 + 
		3688246*t5^23 + 3881994*t5^22 + 5091882*t5^21 + 139190*t5^20 + 
		5420720*t5^19 + 3811588*t5^18 + 4848371*t5^17 + 4707251*t5^16 +
		 4602711*t5^15 + 438552*t5^14 + 1619338*t5^13 + 1960523*t5^12 +
		 2051308*t5^11 + 3427830*t5^10 + 99385*t5^9 + 4859068*t5^8 + 
		649798*t5^7 + 820223*t5^6 + 3699297*t5^5 + 4703466*t5^4 + 
		1979671*t5^3 + 2113839*t5^2 + 1797999*t5 + 584831
		sage: tower.bivariate_to_univariate(randu)==randu
		True
		'''
		if Test==True:
			if self.test_bi(element)==False:
				return element
		j=self.floor(element,-2)
		#image=self._levels[j+1].zero()
		L1=element[0].list()#on recupere les coefficients dans le corps de base de la tour du coefficient de valuation 0 en l element generateur de l extension quadratique
		L2=element[1].list()#on recupere les coefficients du coefficient de valuation 1 en l element generateur de l extension quadratique
		K=self._levels[j+1]
		#g0=self._levels[j+1].gen()
		#g2=g0**2
		#g1=1
		L01=[]		
		for i in range(len(L1)):
			#image+=L1[i]*g1+L2[i]*g1*g0
			#g1*=g2
			L01.append(L1[i])
			L01.append(L2[i])
		#return image
		return K(L01)
		
	def univariate_to_bivariate(self,element,Test=True):
		'''
		Compute the bivariate representation of an element according to it's univariate represenation

		INPUT:
		-self a 2-adic tower
		-element an element who belongs to the 2 adic tower
		
		OUTPUT:
		-a bivariate representation of element

		EXAMPLE:
		sage: A=FiniteField(590295810358705651741)
		sage: tower=Tower_two(A,2,'var')
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: rand=tower._levels[randint(0,len		
		(tower._levels)-1)].random_element()
		sage: tower.test_bi(rand)
		False
		sage: tower.univariate_to_bivariate(rand)
		(438144624500729183827*var5^15 + 120942089090657742635*var5^14 
		+ 156129555714883466203*var5^13 + 426050634633995261833*var5^12
		 + 562135439260484837804*var5^11 + 229147672930488747093*
		var5^10 + 262165143676171352178*var5^9 + 139834227697485604773*
		var5^8 + 143665528252432324217*var5^7 + 582229954009265124834*
		var5^6 + 137247776743245350483*var5^5 + 170707535809781859969*
		var5^4 + 557328708431197174306*var5^3 + 310281612217773889735*
		var5^2 + 378815390545277918934*var5 + 158096202360917180608)*
		var6 + 326545472088197014781*var5^15 + 405314577172505115360*
		var5^14 + 302396889012946540398*var5^13 + 277752148151176116885
		*var5^12 + 61541377272656783930*var5^11 + 45350435637962271005
		*var5^10 + 590281084646677816546*var5^9 + 300203182178169891788
		*var5^8 + 225762999453095092910*var5^7 + 545176520918349701714*
		var5^6 + 246934292990209009038*var5^5 + 497762453521749192182*
		var5^4 + 383287780162870352425*var5^3 + 29379007997975947489*
		var5^2 + 118003439015664062148*var5 + 21880843422525697228
		sage: tower.test_bi(tower.univariate_to_bivariate(rand))
		True
		sage: rand=tower._levels[-1].random_element()
		sage: tower.test_bi(rand)
		False
		sage: tower.test_bi(tower.univariate_to_bivariate(rand))
		True
		sage: rand=tower._levels[-2].random_element()
		sage: tower.test_bi(rand)
		True
		sage: rand==tower.univariate_to_bivariate(rand)
		True
		sage: tower.bivariate_to_univariate		
		(tower.univariate_to_bivariate(rand))==rand
		True
		'''		
		if Test==True:
			if self.test_bi(element):
				return element
		L1=element.list();
		j=self.floor(element,-1)
		#K2=self._levels[j-2]
		K1=self._levels[j-1]		
		#retour1=self._levels[j-1].zero()
		#retour2=retour1
		#g0=self._levels[j-2].gen()
		#g1=1
		L01=[]
		L02=L01
		for i in range(0,len(L1),2):
			#retour1+=L1[i]*g1
			L01.append(L1[i])			
			#retour2+=L1[i+1]*g1
			L02.append(L1[i+1])			
			#g1*=g0
		#retour2=retour2*self._levels[j-1].gen()
		#return retour1+retour2		
		#return K1.gen()*K2(L02)+K2(L01)
		return K1[L01,L02]	
	def test_bi(self,element):
		'''
		Test if the element who belongs to the 2 adic tower self is 
		bivariate.
		
		INPUT:
		-self a 2 adic tower
		-element an element to test if it is bivariate

		OUTPUT:
		A boolean value according to the fact the element is bivariate 
		or not.

		EXAMPLE:
		sage: rand=tower._levels[randint(0,len		
		(tower._levels)-1)].random_element()
		sage: rand.parent()
		Univariate Quotient Polynomial Ring in x1 over Finite Field of 
		size 9223372036854778193 with modulus x1^16 +
		 2968624709650912984
		sage: tower.test_bi(rand)
		False
		sage: rand=tower._levels[randint(0,len
		(tower._levels)-1)].random_element()
		sage: rand.parent()
		Univariate Quotient Polynomial Ring in x2 over Univariate 
		Quotient Polynomial Ring in x1 over Finite Field of size 
		9223372036854778193 with modulus x1^16 + 2968624709650912984 
		with modulus x2^2 + 9223372036854778192*x1
		sage: tower.test_bi(rand)
		True
		'''
		if element.parent()==self._base:
			return False
		else :
			if element[0].parent().gen()==self._base.gen():
				return False
			else:
				return True

	def push1l(self,element):
		'''
		Push the element element in the field below its own in the 
		2-adic tower self

		INPUT:
		-self a 2-adic tower
		-element an element of the 2-adic tower

		OUTPUT:
		The element in a univariate represenation in the field below
		in the 2 adic tower


		EXAMPLE:
		sage: tower=Tower_two(FiniteField		
		(151115727451828646838329),2,'y')
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: rand=tower._levels[2].random_element()
		sage: tower.test_bi(rand)
		True
		sage: tower.push(tower.lift1l(rand))==rand
		False
		sage: tower.push(tower.lift1l		
		(rand))==tower.bivariate_to_univariate(rand)
		True
		sage: tower.push(tower.lift1l(rand))
		86155293705634944065507*y3^7 + 11738399941797272842739*y3^6 + 
		7535855977131842597939*y3^5 + 98757717956189556804059*y3^4 + 
		98164163091034187896974*y3^3 + 10959400026616958963968*y3^2 + 
		18001363466791083306758*y3 + 98760081933253528442968
		sage: rand
		(86155293705634944065507*y1^3 + 7535855977131842597939*y1^2 + 
		98164163091034187896974*y1 + 18001363466791083306758)*y2 + 
		11738399941797272842739*y1^3 + 98757717956189556804059*y1^2 + 
		10959400026616958963968*y1 + 98760081933253528442968
		sage: rand=tower._levels[3].random_element()
		sage: tower.push(tower.lift1l(rand))==rand
		True
		sage: tower.test_bi(rand)
		False
		'''
		if element.parent()==self._levels[0]:
			return False
		elif element.parent()==self._levels[1]:
			return self._levels[0](element[0])	
		elif self.test_bi(element):
			j=self.floor(element,-2)
			#image=self._levels[j+1].zero()
			L1=element[0].list()#on recupere les coefficients dans le corps de base de la tour du coefficient de valuation 0 en l element generateur de l extension quadratique
			#L2=element[1].list()#on recupere les coefficients du coefficient de valuation 1 en l element generateur de l extension quadratique
			K=self._levels[j-2]
			l=len(L1)
			L02=[]
			L01=[]
			for i in range(0,l,2):
				L01.append(L1[i])
				L02.append(L1[i+1])
			#g0=self._levels[j+1].gen()
			#g2=g0**2
			#g1=1
			return K([L01,L02])
		else :
			i=self.floor(element,-1)
			L1=element.list()
			#retour=self._levels[i-2].zero()
			#g0=self._levels[i-2].gen()
			#g1=1
			K=self._levels[i-2]
			L2=[]
			for j in range(0,len(L1),2):
				L2.append(L1[j])
				#retour+=L1[j]*g1
				#g1*=g0
			return K(L2)
		
	def pushD(self,element,nb):
		'''
		Push the element element in the field d levels below in the 
		2-adic tower self

		INPUT:
		-self a 2-adic tower
		-element an element of the 2-adic tower
		-nb the number of level in the tower we want to go down

		OUTPUT:
		The element in a univariate represenation in the field below
		in the 2 adic tower

		EXAMPLE:
		sage: F=FiniteField(101)
		sage: K=Tower_two(F,1)
		sage: K.add_one_level()
		sage: K.add_one_level()
		sage: K.add_one_level()
		sage: K5=K._levels[-7]
		sage: K5
		Univariate Quotient Polynomial Ring in x1 over Finite Field of
		 size 101 with modulus x1^2 + 2
		sage: Test=K5.random_element()
		sage: Test2=Test
		sage: Test2=K.lift1l(Test2); Test2
		78*x3^2 + 49
		sage: Test2=K.lift1l(Test2); Test2
		78*x5^4 + 49
		sage: Test2=K.lift1l(Test2); Test2
		78*x7^8 + 49
		sage: K.pushD(Test2,3)
		78*x1 + 49
		sage: K.pushD(Test2,3)==Test
		True
		'''
		i=0
		Test=True
		while (Test and i<nb ):
			i=i+1 
		#for i in range(nb):	
			element2=self.push1l(element)
			if element==self.lift1l(element2):
				element=element2
			else:
				Test=False
		return element

	def cardinality_field(self,field):
		'''
		Compute the cardinality of a field on the curve
		
		INPUT:
		-field a finite field who belongs to the 2-adic tower self
		-self a 2-adic tower
		
		OUTPUT:
		The cardinality of the field		
		
		EXAMPLE:
		sage: tower=Tower_two(FiniteField		(304835101),4,'x')
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: tower.add_one_level()
		sage: a=tower._levels[randint(0,len(tower._levels)-1)]
		sage: a
		Univariate Quotient Polynomial Ring in x11 over Finite Field of
		 size 304835101 with modulus x11^512 + 66020696
		
		sage: tower.cardinality_field(a)==a.cardinality()
		True
		'''
		if field==self._levels[0]:
			return self._levels[0].cardinality()
		else:			
			i=-1
			while (self._levels[i]!=field):
				i-=1
			i=int((-i-1)/2)
			retour=self._top_cardinal
			while i!=0:
				i-=1
				a,b=retour.factor()[0]
				retour=a**(b/2)
			return retour

	def lift1l(self,element):
		'''
		Lift the element element in the field above its own in the 
		2-adic tower self

		INPUT:
		-self a 2-adic tower
		-element an element of the 2-adic tower

		OUTPUT:
		The element in a univariate represenation in the field above
		in the 2 adic tower


		EXAMPLE:
		sage: tower=Tower_two(FiniteField(274877906957),6,'x')
		sage: rand=tower._top.random_element()
		sage: rand
		81167568574*x1^63 + 5708700412*x1^62 + 141398047299*x1^61 + 
		40017253015*x1^60 + 146498969736*x1^59 + 238819366038*x1^58 + 
		120548415538*x1^57 + 110956335376*x1^56 + 12180378314*x1^55 + 
		203792127200*x1^54 + 129443701333*x1^53 + 224628847406*x1^52 + 
		235425337778*x1^51 + 11153815859*x1^50 + 128390527262*x1^49 + 
		218967627437*x1^48 + 133261544663*x1^47 + 119727481338*x1^46 + 
		90492542981*x1^45 + 34229939685*x1^44 + 169346742904*x1^43 + 
		155056051545*x1^42 + 150103700729*x1^41 + 16755034999*x1^40 + 
		110473761394*x1^39 + 71124840718*x1^38 + 224217049489*x1^37 + 
		247365080360*x1^36 + 263862310336*x1^35 + 16598744618*x1^34 + 
		127560288532*x1^33 + 273176473703*x1^32 + 245480895805*x1^31 + 
		94583737710*x1^30 + 169932652192*x1^29 + 98861728566*x1^28 + 
		56564435349*x1^27 + 95724363943*x1^26 + 231914839762*x1^25 + 
		117748059792*x1^24 + 261794620682*x1^23 + 147255281320*x1^22 + 
		269450523909*x1^21 + 201281363514*x1^20 + 51778156992*x1^19 + 
		130714594642*x1^18 + 32622075949*x1^17 + 206387000135*x1^16 + 
		127084265343*x1^15 + 227530778577*x1^14 + 33114450197*x1^13 + 
		39684548081*x1^12 + 207945037531*x1^11 + 234452740845*x1^10 + 
		229996467388*x1^9 + 236309029388*x1^8 + 92804606257*x1^7 + 
		169902980975*x1^6 + 70276057354*x1^5 + 64893131660*x1^4 + 
		58959369247*x1^3 + 188626126231*x1^2 + 62805321683*x1 + 
		133178898807
		sage: tower.lift1l(rand)
		81167568574*x3^126 + 5708700412*x3^124 + 141398047299*x3^122 + 
		40017253015*x3^120 + 146498969736*x3^118 + 238819366038*x3^116 
		+ 120548415538*x3^114 + 110956335376*x3^112 + 12180378314*
		x3^110 + 203792127200*x3^108 + 129443701333*x3^106 + 
		224628847406*x3^104 + 235425337778*x3^102 + 11153815859*x3^100 
		+ 128390527262*x3^98 + 218967627437*x3^96 + 133261544663*x3^94 
		+ 119727481338*x3^92 + 90492542981*x3^90 + 34229939685*x3^88 + 
		169346742904*x3^86 + 155056051545*x3^84 + 150103700729*x3^82 + 
		16755034999*x3^80 + 110473761394*x3^78 + 71124840718*x3^76 + 
		224217049489*x3^74 + 247365080360*x3^72 + 263862310336*x3^70 + 
		16598744618*x3^68 + 127560288532*x3^66 + 273176473703*x3^64 + 
		245480895805*x3^62 + 94583737710*x3^60 + 169932652192*x3^58 + 
		98861728566*x3^56 + 56564435349*x3^54 + 95724363943*x3^52 + 
		231914839762*x3^50 + 117748059792*x3^48 + 261794620682*x3^46 + 
		147255281320*x3^44 + 269450523909*x3^42 + 201281363514*x3^40 + 
		51778156992*x3^38 + 130714594642*x3^36 + 32622075949*x3^34 + 
		206387000135*x3^32 + 127084265343*x3^30 + 227530778577*x3^28 + 
		33114450197*x3^26 + 39684548081*x3^24 + 207945037531*x3^22 + 
		234452740845*x3^20 + 229996467388*x3^18 + 236309029388*x3^16 + 
		92804606257*x3^14 + 169902980975*x3^12 + 70276057354*x3^10 + 
		64893131660*x3^8 + 58959369247*x3^6 + 188626126231*x3^4 + 
		62805321683*x3^2 + 133178898807
		sage: tower.push(tower.lift1l(rand))==rand
		True
		'''
		if element.parent()==self._levels[0]:
			return self._levels[1](element)
		if self.test_bi(element):
			j=self.floor(element,-2)
			if j==-2:
				self.add_one_level();
				j=j-2
			#image=self._levels[j+1].zero()
			L1=element[0].list()#on recupere les coefficients dans le corps de base de la tour du coefficient de valuation 0 en l element generateur de l extension quadratique
			L2=element[1].list()#on recupere les coefficients du coefficient de valuation 1 en l element generateur de l extension quadratique
			K=self._levels[j+2]
			#g0=self._levels[j+1].gen()
			#g2=g0**2
			#g1=1
			L01=[]
			L02=[]		
			for i in range(len(L1)):
				#image+=L1[i]*g1+L2[i]*g1*g0
				#g1*=g2
				L01.append(L1[i])
				L02.append(0)
				L01.append(L2[i])
				L02.append(0)
			#return image
			return K([L01,L02])
		else :	
			i=self.floor(element,-1)
			if i==-1:
				self.add_one_level();
				i=-3
			K=self._levels[i+2]
			L01=[]
			for j in element:
				L01.append(j)
				L01.append(0)	
			return K(L01)

	def lift(self,element):
		'''
		Lift the element element in the representation of the field above its own in the 
		2-adic tower self

		INPUT:
		-self a 2-adic tower
		-element an element of the 2-adic tower

		OUTPUT:
		The element in the corresponding representation in the representation of the field above its own
		in the 2 adic tower

		If we have a univariate representation of element then we will have an element in the quadratic extension, whereas if we have a bivariate representation of element then we will have it's representation in the univariate form

		EXAMPLES: 
		'''
		if element.parent()==self._levels[0]:
			return self._levels[1](element)
		elif self.test_bi(element):
			j=self.floor(element,-2)
		#image=self._levels[j+1].zero()
			L1=element[0].list()#on recupere les coefficients dans le corps de base de la tour du coefficient de valuation 0 en l element generateur de l extension quadratique
			L2=element[1].list()#on recupere les coefficients du coefficient de valuation 1 en l element generateur de l extension quadratique
			K=self._levels[j+1]
			#g0=self._levels[j+1].gen()
			#g2=g0**2
			#g1=1
			L01=[]		
			for i in range(len(L1)):
				#image+=L1[i]*g1+L2[i]*g1*g0
				#g1*=g2
				L01.append(L1[i])
				L01.append(L2[i])
			#return image
			return K(L01)
		else :
			i=self.floor(element,-1)
			K=self._levels[i+1]
			L1=[]
			L2=[]
			for j in element:
				L1.append(j)
				L2.append(0)
			return K([L1,L2])

	def push(self,element):
		'''
		Push the element element in the representation of the field below its own in the 
		2-adic tower self

		INPUT:
		-self a 2-adic tower
		-element an element of the 2-adic tower

		OUTPUT:
		The element in the corresponding representation in the representation of the field below its own
		in the 2 adic tower

		If we have a bivariate representation of element then we will have an element in the field below in the 2-adic tower in the  univariate representation, whereas if we have a univariate representation of element then we will have it's representation in the bivariate form

		EXAMPLES: 
		'''
		if element.parent()==self._levels[0]:
			return False
		elif element.parent()==self._levels[1]:
			return self._levels[0](element[0])
		elif self.test_bi(element):
			i=self.floor(element,-2)
			L1=element[0].list()
			K=self._levels[i-1]
			return K(L1)
		else :
			#return self.univariate_to_bivariate(element,Test=False)
			L1=element.list();
			j=self.floor(element,-1)
			#K2=self._levels[j-2]
			K1=self._levels[j-1]		
			#retour1=self._levels[j-1].zero()
			#retour2=retour1
			#g0=self._levels[j-2].gen()
			#g1=1
			L01=[]
			L02=[]
			for i in range(0,len(L1),2):
				#retour1+=L1[i]*g1
				L01.append(L1[i])			
				#retour2+=L1[i+1]*g1
				L02.append(L1[i+1])			
				#g1*=g0
			#retour2=retour2*self._levels[j-1].gen()
			#return retour1+retour2		
			#return K1.gen()*K2(L02)+K2(L01)
			return K1([L01,L02])

	def determination_level(self,element):
		'''
		Compute the level which the element belongs to
		'''
		if element.parent()==self._levels[0]:
			return element,0
		i=self.floor(element)
		Test=True
		while(Test==True and element.parent()!=self._base and element.parent()!=self._levels[1]):
			temp=self.push1l(element)		
			if self.lift1l(temp)!=element :
				Test=False
			else :
				element=temp
				i=i-2
		if element.parent()==self._levels[1]:
			q=self._levels[0].cardinality()			
			if element**q==element:
				element=self._levels[0](element[0])
				return element,0
			return element,1
		return element,i
	def init_zeta(self,lambda_,q,m):
		L=[0]
		lambda_=self.frobenius_computation(lambda_,q)
		L.append(lambda_)
		for i in range(2,m-1):
			if i%2==1:
				l2=self.frobenius_computation(L[-1],q)
				L.append(L[1]*l2)
			else :
				l2=self.frobenius_computation(L[i/2],q**(i/2))
				L.append(L[i/2]*l2)
		return L

		

	def computing_em(self,lambda_,m,q,L):
		if m==0:
			print 'probleme m=0'
			return 'probleme m=0'
		elif m==1:
			return L[1]
		else :
			if m%2==1:
				return self.computing_em(lambda_,m-1,q,L)+L[m]
			else :
				return  self.computing_em(lambda_,m/2,q,L)+L[m/2]*self.frobenius_computation(self.computing_em(lambda_,m/2,q,L),q**(m/2))

	def root_computing(self,delta):
		'''
		Input:
		-delta an element of the tower self
		-self the 2-adic tower

		Output:
		The square root of delta without checking if delta is a square.
		If delta is not a square it will gave it's root in the level 
		above in the 2-adic tower.

		Examples:
		'''
		test=False
		c=1
		q=self._base.cardinality()
		delta,i=self.determination_level(delta) #permet de se placer dans le corps le plus petit au sens de l'inclusion
		n=valuation(self.cardinality_field(delta.parent()),q)
		if n==1:
		#we are looking if we are working on a basefield, id est a finite field in which root can be computed directly by sage here
			if self.test_bi(delta):#devrait etre inutile
				delta=delta[0]	
			#delta=self._base(delta[0])
			if delta.is_square():
				return delta.square_root()
			else :
				beta2=delta/self._alpha
				beta=beta2.square_root()
				eta=valuation(self.cardinality_field(self._levels[1]),q)
				return beta*self._levels[1].gen()**(eta/2)
		elif n==2:
			while test==False:
				lambda_=delta**((q-1)/2)
				eta=1+lambda_
				if eta==0:
					c*=delta.parent().random_element()
					delta*=c**2
				else:
					test=True						
		else:
			while test==False:		
				lambda_=delta**((q-1)/2)
				L=self.init_zeta(lambda_,q,n)
				em=self.computing_em(lambda_,n-2,q,L)
				eta=1+lambda_*(1+em)
				if eta==0:
					c*=delta.parent().random_element()
					delta*=c**2
				else:
					test=True
		beta2=delta*(eta**2)
		i=1
		if beta2**q!=beta2: #on teste si beta appartient a Fq
			delta=delta*delta.parent().gen()#on multiplie delta par un residu non quadratique dont on connait la racine dans l extension d apres la construction
			if delta.parent()==self._levels[-1] or delta.parent()==self._levels[-2]:#si on est en haut de la tour on ajoute un etage pour trouver la racine
				self.add_one_level()
			#delta=self.lift1l(delta) on deplace cela ailleurs
			test=False
			while test==False:
				if n==2:
					lambda_=delta**((q-1)/2)
					eta=1+lambda_
					if eta==0:
						c*=delta.parent().random_element()
						delta*=c**2
					else:
						test=True										
				
				else:						
					lambda_=delta**((q-1)/2)
					L=self.init_zeta(lambda_,q,n)
					em=self.computing_em(lambda_,n-2,q,L)
					eta=1+lambda_*(1+em)
					if eta==0:
						c*=delta.parent().random_element()
						delta*=c**2
					else:
						test=True		
			beta2=delta*(eta**2)
			delta=self.lift1l(delta) #on se place dans l extension superieur ou se trouve la racine quadratique
			eta=self.meeting2(eta,delta)
			i=delta.parent().gen()#etape importante on stocke l information que l on s est servi d un residu quadratique		
		if self.test_bi(beta2):
			beta2=beta2[0]
		beta2=self._base(beta2[0])
		if beta2.is_square():
			beta=beta2.square_root()
		else: #ce else doit etre inutile
			beta2=beta2/self._alpha
			beta=beta2.square_root()
			beta2=valuation(self.cardinality_field(self._levels[1]),q)
			beta=beta*self._levels[1].gen()**(beta2/2) #on enleve dans le resultat la partie qui vient de alpha
			beta=self.meeting2(beta,delta)#on releve le resultat dans une extension superieure si necessaire afin que delta et beta soient au meme niveau	
		if c!=1:
			c=self.meeting2(c,eta)
		beta=(delta.parent()(1)/(eta*c*i))*beta #il est important de noter que l on divise par i qui contient (ou pas) la partie que l on a introduite avec le residu quadratique avec lequel on avait multiplie delta plus haut
		#beta2=self.push1l(beta) inutile maintenant que determination a ete change
		#if self.lift1l(beta2)==beta:
			#beta=beta2
		return beta

	def calcul_beta(self,delta): #fonction qui ne sert plus a rien
		test=False
		c=1
		q=self._base.cardinality()
		#K1=self.lift1l(delta).parent()# a supprimer
		deltab=delta*delta.parent().gen()
		delta,i=self.determination_level(delta) #permet de se placer dans le corps le plus adapte pour etre sur de trouver une racine
		n=valuation(self.cardinality_field(delta.parent()),q)
		if n==1:
			if self.test_bi(delta):
				delta=delta[0]	
			delta=self._base(delta[0])
			if delta.is_square():
				return delta.square_root()
			else :
				beta2=delta/self._alpha
				beta=beta2.square_root()
				return beta*self._levels[1].gen()
		else:
			while test==False:		
				lambda_=delta**((q-1)/2)
				L=self.init_zeta(lambda_,q,n)
				em=self.computing_em(lambda_,n-2,q,L)
				eta=1+lambda_*(1+em)
				if eta==0:
					c*=delta.parent().random_element()
					delta*=c
				else:
					test=True
		beta2=delta*(eta**2)
		lambda_=1
		if beta2**q!=beta2:#on avait en fait pas pris un carre dans le corps
			delta,i=self.determination_level(deltab)			
			test=False			
			while test==False:		
				lambda_=delta**((q-1)/2)
				L=self.init_zeta(lambda_,q,n)
				em=self.computing_em(lambda_,n-2,q,L)
				eta=1+lambda_*(1+em)
				if eta==0:
					c*=delta.parent().random_element()
					delta*=c
				else:
					test=True
			lambda_=delta.parent().gen()
			beta2=delta*(eta**2)
		if self.test_bi(beta2):
			beta2=beta2[0]
		if self._base(beta2[0]).is_square():
			beta=self._base(beta2[0]).square_root()
		else:
			beta=self.root_computing(self._levels[0](beta2[0]))
			#il faut mettre beta au meme niveau que eta et c
			print "test racine carre", 'beta2[0]', beta2[0], 'beta2', beta2
			if i!=1:
				if self.test_bi(delta):
					beta=self.lift(beta)
				K=self._levels[i]
				while beta.parent()!=K:
					beta=self.lift1l(beta)
			#print beta, 'beta.parent()==eta.parent()', beta.parent()==eta.parent()
		#if beta2!=beta**2:
			#print 'beta2',beta2,'beta**2',beta**2
		beta=(delta.parent()(1)/(eta*c*lambda_))*beta
		beta2=self.push1l(beta)
		if self.lift1l(beta2)==beta:
			beta=beta2
		return beta

	def meeting(self,element1,element2):
		'''
		Lift une des deux entrees afin que les deux elements soient au meme niveau
		'''
		i1=self.floor(element1)
		i2=self.floor(element2)
		if (i1==0 and i2!=0):
			element1=self._levels[1](element1)
			i1=self.floor(element1)
		if (i2==0 and i1!=0):
			element2=self._levels[1](element2)
			i2=self.floor(element2)
		r=False
		if i1<i2 :
			r=element1
			element1=element2
			element2=r
			r=i1
			i1=i2
			i2=r
			r=True
		if self.test_bi(element1)!=self.test_bi(element2):
			element2=self.lift(element2)
			i2+=1
		while i2<i1:
			element2=self.lift1l(element2)
			i2+=2
		if r:
			return element2,element1
		else:
			return element1,element2

	def meeting2(self,element1,element2):
		'''
		Lift la premiere entree afin que les deux elements soient au meme niveau
		'''
		i1=self.floor(element1)
		i2=self.floor(element2)
		if i1<i2 :
			if self.test_bi(element1)!=self.test_bi(element2):
				element1=self.lift(element1)
				i1+=1
			while i1<i2:
				element1=self.lift1l(element1)
				i1+=2
		return element1
