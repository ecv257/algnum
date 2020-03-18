from algnum import *
from operator import itemgetter
from fractions import Fraction
from math import sqrt

class quadnum(tuple):
	''' representamos un entero cuadrático por una tupla $(b,c,D)= (b+c\sqrt D)/2$ donde $D$ es un discriminante no necesariamente fundamental ie $D\equiv 0,1 (4)$ no es un cuadrado perfecto.
	'''
	
	'''
	Con __slots__ indicamos los atributos quepueden accederse usando ['..'] o ,'...'. Al indicar
	una lista vacñia se indica que el objeto no tiene atributos modificables desde fuera. dicho de otra
	manera es la forma más sencilla de crear clases con elementos inmutables.'''
	__slots__ = []
	
	''' La forma de crear elementos nuevos es bastante sorprendente y no acabo de entenderla.'''
	def __new__(cls,b,c,D):
		return tuple.__new__(cls,(
		b.numerator if type(b)==Fraction and b.denominator == 1 else b,
		c.numerator if type(c)==Fraction and c.denominator == 1 else c,
		D))
	b = property(itemgetter(0))
	c = property(itemgetter(1))
	D = property(itemgetter(2))
	
	@property
	def beta(self):
		return 2*self.b+(self.c if self.D%4 ==1 else 0)
	
	@property
	def tau(self):
		return sign(self.c)
	
	def __add__(self,a):
		if type(a) == quadnum:
			if a.D == self.D: return quadnum(self.b+a.b,self.c+a.c,self.D)
			else: return NotImplemented
		else:
			return quadnum(self.b+a,self.c,self.D)
	def __radd__(self,a):
		return self.__add__(a)
	def __sub__(self,a):
		return self.__add__(-a)
	def __rsub__(self,a):
		return (-self).__add__(a)
	def __pos__(self):
		return self
	def __neg__(self):
		return quadnum(-self.b,-self.c,self.D)
	def __mul__(self,a):
		if type(a) == quadnum:
			if a.D == self.D:
				return quadnum(self.b*a.b+self.c*a.c*(self.D if self.D%4==0 else (self.D-1))//4, self.b*a.c+self.c*a.b+(0 if self.D%4==0 else self.c*a.c),self.D)
			else:
				return NotImplemented
		else:
			return quadnum(self.b*a,self.c*a,self.D)
	def __rmul__(self,a):
		return self.__mul__(a)
	def __truediv__(self,a):
		if type(a) == quadnum:
			if a.D == self.D:return (self*a.__invert__())
			else: return NotImplemented
		else:
			return self*(Fraction(1)/a)
	def __rtruediv__(self,a):
		return a*self.__invert__()
	def __pow__(self,a):
		if type(a)== int :
			return generic_pow(self,a)
		return NotImplemented
	def __invert__(self):
		return self.conj/self.norm
	def __int__(self):
		if self.c == 0 and self.b==int:
			return self.b
		else:
			return NotImplemented
	def __repr__(self):
		return '(' + str(self.b)+','+ str(self.c)+','+str(self.D)+')'
	def __str__(self):
		b = self.b if self.D%4==0 else Fraction(2*self.b+self.c)/2
		c = self.c  if self.D%4==0 else Fraction(self.c)/2
		d = self.D//4 if self.D%4==0 else self.D
		return str(b)+(('+' if c>0 else '')+(str(c) if c !=1 else '')+'v'+str(d) if self.c else '')
	@property
	def norm(self):
		return (self.b*self.b-self.c*self.c*self.D//4 ) if self.D%4==0 else (self.b*self.b+self.b*self.c+self.c*self.c*(1-self.D)//4 )
	@property
	def conj(self):
		return quadnum(self.b+(self.c if self.D%4==1 else 0),-self.c,self.D)
	def __abs__(self):
		#devuelve un float/complex con el valor numérico
		return (sqrt(self.D)*self.c+2*self.b+(self.c if self.D%4==1 else 0))/2
	def __eq__(self,a):
		if type(a)==int or type(a)==Fraction :
			return a==self.b and self.c ==0
		elif type(a)==quadnum:
			if self.D==a.D:
				return self.b==a.b and self.c==a.c
			elif self.b==a.b:
				g = gcd(self.D,a.D)
				if issquare(self.D//g) and issquare(a.D//g):
					return self.c*isqrt(self.D//g)==a.c*isqrt(a.D//g)
		else:
			return super().__eq__(self,a)
	def __lt__(self,a) :
		if self.D<0:
			return NotImplemented
		q=self-a
		b=2*q.b if q.D%4==0 else 2*q.b+q.c
		if b <0 and q.c < 0:
			return True
		if b >0 and q.c > 0:
			return False
		return  b*abs(b)+q.c*abs(q.c)*q.D<0
	def __gt__(self,a) :
		if self.D<0:
			return NotImplemented
		q=self-a
		b=2*q.b if q.D%4==0 else 2*q.b+q.c
		if b <0 and q.c < 0:
			return False
		if b >0 and q.c > 0:
			return True
		return  b*abs(b)+q.c*abs(q.c)*q.D>0
	def __int__(self):
		return (sign(self.c)*isqrt(self.c**2*self.D)+2*self.b+(self.c if self.D%4==1 else 0))//2
		
class mo:
	''' 
	Representamos un modulo orientado por una expresion $[(b \pm \sqrt D)/2, a]$ donde $b\equiv D \pmod 2$ y $a$ es un divisor de la norma $(b^2-D)/4$
	
	Internamente lo representamos como una tupla (b,tau,D,a)
	el método .q recupera el quadnum $(b \pm \sqrt D)/2$ observemos que la representación de q es diferente: q=w
	'''			
	def __init__(self,b,tau,D,a):
		self.b=b
		self.tau=tau
		self.D=D
		self.a=a
		self.d=isqrt(D)
		
	def __init__(self,b,tau,D,a,d):
		self.b=b
		self.tau=tau
		self.D=D
		self.a=a
		self.d=d

	@property
	def q(self):
		return quadnum()

	def check(self):
		return (self.qnorm % self.a ==0)
	
	@property
	def qnorm(self):
		return (self.b**2-self.D)//4
	
	def isreduced(self):
		n=self.qnorm
		if n>=0:
			return False
		ap=n//self.a
		return self.tau*(self.b+2*self.a) > self.d and self.tau*(self.b-2*ap) > self.d

	def next(self):
		k = (self.b+self.tau*self.d)//(2*self.a)
		self.b -= 2*self.a*k
		self.a = self.qnorm//self.a
		self.tau=-self.tau
		return self
	
	def cycle(self)	:
		while not self.isreduced():
			self.next()
		cy = [(self.b,self.a)]
		self.next()
		while (self.b,self.a)!=cy[0]:
			cy.append((self.b,self.a))
			self.next()
		return cy

	def __str__(self):
		b = self.b//2 if self.D%4==0 else self.b
		d = self.D//4 if self.D%4==0 else self.D
		if self.D%4==0:
			return f'[{b}+\u221a{d},{self.a}]'
		return f'[({b}+\u221a{d})/2,{self.a}]'
		 

class classGroup:
	def __init__(self, D):
		self.D=D
		d=isqrt(D)
		self.C=[]
		self.G={}
		self.fD=factor(abs(D))
		z=set()
		for b in range(d-(D-d)%2, 0, -2):
			n= abs(b*b-D)//4
			for a in factor(n).divisors():
				if gcd(b,gcd(a,n//a))==1 and (b,a) not in z:
					m = mo(b,1,D,a,d)
					if m.isreduced():
						c = m.cycle()
						self.C.append(c)
						z.update(c)
	@property
	def h(self):
		return len(self.C)
	def isambiguos(self,k=0,ldiv=[]):
		if ldiv==[]:ldiv=set(self.fD.divisors())
		g=[]
		for b,a in self.C[k]: 
			if abs(a) in ldiv and b%a==0:
				g.append(a)
		return g
	def ambiguos(self):
		self.G = {}
		ldiv = list(factor(self.D).divisors())
		for i,c in enumerate(self.C):
			g = self.isambiguos(i,ldiv)
			if g:
				self.G[i]=g
		return self.G
	def unit(self,k=0):
		b1,a1=self.C[k][0]
		q1=-quadnum(b1,1,self.D)
		a1=-a1
		q=quadnum(1,0,self.D)
		n=1
		for b,a in self.C[k]:
			q2=quadnum(b,sign(a),self.D)
			if q2==q1 and a==a1:
				break
			q=q*q2
			n=n*a
		return q/n
		

test =1
if test:
	m1 =mo(7, 1, 37, 3, 6);
	c1=m1.cycle() 
	m2 =mo(-4, -1, 28, -3, 5);c2=m2.cycle()
	
	for d in range(1,300):
		if  not issquare(d) and d%4==0:
			F=factor(d)
			f=F.sqfreekern()
			if f%4==1 and f!=d:
				H=classGroup(d)
				h2=len(H.ambiguos())
				if h2 != 2**(len(F.P)-1): #len(factor(h2).P)>1:
					print(f'{d} ({factor(d)}):{H.h} {h2} {H.ambiguos()}')

	g=classGroup(148) # no cuadra el númerode genera
	g2=classGroup(2005) #
	g3=classGroup(229)


