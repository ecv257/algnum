# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 11:59:31 2019

@author: ecv25
"""
test = 2
from fractions import Fraction
from decimal import Decimal
import operator

#definimos la función signo 
sign = lambda x: (1, -1)[x<0]

# In[]
def generic_inverse(g):
    print(type(g))
    if type(g) == int:
        return 1/Fraction(g)
    elif type(g) == Fraction:
        return 1/g
    elif type(g) == Decimal:
        return 1/g
    else:
        return g.inverse()
        
def generic_pow(g, n):
    if isinstance(n, int):
        if n == 0 : return type(g)()
        if n < 0 : 
            n = -n 
            g = generic_inverse(g)
        y = g
        f = n.bit_length()-2
        while f >= 0 :
            y = y*y
            if n>>f&1 : 
                y = y*g
            f -= 1
        return y
    else:
        return NotImplemented

#generic_pow(Decimal(1)/3, -7)-(Decimal(1)/3)**-7

# In[] 
def bezout(a,b):
    '''**Algoritmo Extendido de Euclides** o **Algoritmo de Bezout**:
    Dados enteros a,b devuelve enteros x,y tales que ax+by=gcd(a,b). 
    Usamos el algoritmo binario sigue el algoritmo 1.3.8 de Cohen'''
    if type(a)!= int or type(b) != int: 
        return generic_bezout(a,b)
    a, b, sa, sb  = abs(a),abs(b), sign(a), sign(b)
    if b== 0: return sa,0,a
    if a== 0: return 0,sb,b
    if a<b: a,b,f1 = b,a,1;
    else: f1=0
    # haceos un paso del agoritmo de Euclides para eliminar casos aberrantes
    q,r = a//b,a%b
    a,b = b,r
    if b== 0: return ( 0,sb,a ) if f1 == 0 else (sa,0,a)
    k=0
    while not (a&1 | b&1): # mientras a y b sean pares
        k+=1
        a>>=1;b>>=1
    if b&1: f2 = 0
    else: a,b,f2=b,a,1
    u,d,v1,v3 = 1,a,b,b
    if a&1: t1,t3=0,-b
    else: t1,t3 = (1+b)>>1,a>>1
    while t3 != 0:
        while not t3&1:
            t3>>=1
            t1 = (t1+b)>>1 if t1&1 else t1>>1
        if t3 > 0:
            u, d = t1,t3
        else:
            v1, v3 = b-t1,-t3
        t1,t3 = u-v1,d-v3
        if t1<0: t1 = t1+b
    v = (d-a*u)//b 
    d = d<<k
    if f2 == 1: u,v=v,u
    u = u-q*v
    return (sa*u,sb*v,d) if f1 else (sa*v,sb*u,d)

### Para el test de funcionamiento
if test == 1:
    import random
    import math
    import fractions
    M=1000
    for t in range(1,10):
        a=random.randint(-M,M)
        b=random.randint(-M,M)
        c=random.randint(1,M)
        d=random.randint(1,M)
        a = fractions.Fraction(a,c)
        b = fractions.Fraction(b,d)
        #u,v,g2=bezout(a,b)
        g = fractions.gcd(a,b)
        uu,vv,gg2=generic_bezout(a,b)
        print(f'a={a},b={b},g={g}, u={uu}, v={vv}, ua+vb={uu*a+vv*b}, g2={gg2}')



# In[]
class primes:
    def initprimetable(ulim):
        size = (ulim>>1)+(1<<(ulim.bit_length()+1>>1))+3
        dpr = bytearray(size)
        k,q,r,s,fin = 1,0,0,0,(ulim>>1)
        while r <= fin:
            q = dpr.find(0,q+1)
            k = (q<<1)+1
            r = (q*q+q)<<1
            for s in range(r,fin+1,k):
                dpr[s]=1
        r,q =1,1
        while True:
            s=q
            q = dpr.find(0,q+1)
            if q>fin:
                break
            r+=1
            dpr[r] = (q-s) << 1
        dpr[0]=2;dpr[1]=1 #2 and 3
        nprimos = r+1;
        dpr[r+1]=0;
        lastp = ( s << 1) + 1
        return bytes(dpr[:r+2]), nprimos, lastp
    dpr, np, lp = initprimetable(2000000)
    def __init__(self, lim = 0, start=0):
        self.lim = min(lim,self.lp) if lim >0 else self.lp
        self.n = 1; self.p = self.dpr[0]
        while self.p < start:
            self.p += self.dpr[self.n]
            self.n += 1

    def __iter__(self):
        while self.p < self.lim:
            yield self.p
            self.p += self.dpr[self.n]
            self.n += 1

    def __getitem__(self,n):
        if type(n) == int:
            if self.n > n:
                self.n = 1
                self.p = 2
            while self.n < n:
                self.p += self.dpr[self.n]
                self.n = self.n+1
            return self.p
        elif type(n) == slice:
            start = n.start if n.start != None else 1
            stop = n.stop+1 if n.stop != None else self.np+1
            step = n.step if n.step != None else 1
            return list(self.__getitem__(m) for m in iter(range(start, stop, step)))
# In[]
        
from heapq import heappush,heappop

def isqrt(n):
    ''' calcula la raiz cuadrada entera de n usando el método de Newton '''
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x

class factor:
    def __init__(self, n):
        self.P = []
        if type(n) == Fraction:
            fd = n.denominator
            n = n.numerator
        else : 
            fd = 1
        if type(n) == int:
            # trial division hasta el limite de primos
            if n < 0 :
                self.P.append((-1,1))
                n = -n
            for p in primes():                
                k = 0
                while n%p == 0:
                    k += 1
                    n //=p
                if k > 0:
                    self.P.append((p,k))
                if p*p > n :
                    break
            if n > 1:
                self.P.append((n,1))
            if fd > 1:
                n = fd
                for p in primes():        
                    k = 0
                    while n%p == 0:
                        k += 1
                        n //=p
                    if k > 0:
                        self.P.append((p,-k))
                    if p*p > n :
                        break
                if n > 1:
                    self.P.append((n,-1))
                self.P.sort()

    def sumdiv(self):
        a = 1
        for b in self.P:
            a *= (b[1]**(b[2]+1)-1)//(b[1]-1)
        return a
    
    def numdiv(self):
        a = 1
        for b in self.P:
            a *= (b[2]+1)
        return a

    def phi(self):
        a = 1
        for b in self.P:
            a *= (b[1]-1)*b[1]**(b[2]-1)
        return a

    def mu(self):
        a = 1
        for b in self.P:
            if b[2]>1: 
                a = 0
                break
            else:
                a = -a
        return a

    def issquare(self):
        pass
    
    def issquarefree(self):
        pass
    
    def ispowerfree(self,m):
        pass
    
    def ispower(self):
        pass


    def divisors(self):
        '''Devuelve los divisores de n en orden creciente
        Para ello mantenemos un heap con ternas (d,p,k) d es el divisor que estamos
        analizando, p^k es la mayor potencia de p que lo divide, y agregamos a la cola
        dq para cada q primo >= p con q==p sii p^{k+1}\\n '''
        if len(self.P):
            self.h = [(1,0,0)]
            while len(self.h):
                yield self.h[0][0]
                m, b, k  = heappop(self.h)
                if self.P[b][1] > k:
                    heappush(self.h,(m*self.P[b][0],b,k+1))
                while b+1 < len(self.P):
                    b+=1
                    heappush(self.h,(m*self.P[b][0],b,1))
        else:
            yield 1

list(factor(120).divisors())

# In[]
from operator import itemgetter
class mod(tuple):
    
    
    ''' Con __slots__ indicamos los atributos qu epueden accederse usando ['..'] o ,'...'. Al indicar 
    una lista vacñia se indica que el objeto no tiene atributos modificables desde fuera. dicho de otra 
    manera es la forma más sencilla de crear clases con elementos inmutables.'''
    __slots__ = []
    
    ''' La forma de crear elementos nuevos es bastante sorprendente y no acabo de entenderla.'''
    def __new__(cls,a,m):
        return tuple.__new__(cls,(a%m,m))
    a = property(itemgetter(0))
    m = property(itemgetter(1))
    
    def __add__(self,a):
        if type(a) == mod and a.m == self.m:
            return mod(self.a+a.a,self.m)
        else:
            return mod(self.a+a,self.m)
    def __radd__(self,a):
        return mod(self.a+a,self.m)
    def __sub__(self,a):
        if type(a) == mod and a.m == self.m:
            return mod(self.a-a.a,self.m)
        else:
            return mod(self.a-a,self.m)
    def __rsub__(self,a):
        return mod(a-self.a,self.m)
    def __mul__(self,a):
        if type(a) == mod and a.m == self.m:
            return mod(self.a*a.a,self.m)
        else:
            return mod(self.a*a,self.m)
    def __rmul__(self,a):
        return mod(self.a*a,self.m)
    def __truediv__(self,a):
        if type(a) == mod and a.m == self.m:
            return self*a.__invert__()
        else:
            return mod(self.a*(~mod(a,self.m)).a,self.m)
    def __rtruediv__(self,a):
        return mod(a*self.__invert__().a,self.m)
    def __pow__(self,a):
        if type(a)== int :
            return generic_pow(self,a)
        return NotImplemented
    def __invert__(self):
        u,v,d = bezout(self.a,self.m)
        if abs(d) == type(self.a)(1):
            return mod(d*u,self.m)
    def __int__(self):
        if a.type == int:
            return self.a
        return NotImplemented
    def __repr__(self):
        return 'mod(' + str(self.a)+','+str(self.m)+')' 

# In[]

def chinese(*args):
    ''' los argumentos son de la forma (x_i,m_i), Calcula la única clase  ... '''
    x,m = 1,1
    for arg in args :
        u,v,d = bezout(m,arg.m)
        x = u*m*arg.a + v*arg.m*x
        m = m*arg.m
        x = x%m
    return mod(x,m)

# In[]
def Gauss(a,b):
    ''' a e b son vectores (ndarrays) en un espacio euclídeo determina la combinación lineal entera de X e Y 
    con norma mínima'''
    A,B = np.dot(a,a),np.dot(b,b)
    if A < B: A,B,a,b = B,A,b,a
    while True:
        n = np.dot(a,b)
        r = int(round(n/B))
        T = A-2*r*n+r*r*B
        if T >= B: break
        t,a = a-r*b,b
        b = t;A = B;B = T
    return b

# In[]
class poly(tuple):
    ''' Definimos la clase de polinomios de nuevo como objetos inmutables. Basados en un tupla, y una 
    'variable' que puede ser cualquier objeto. (Podría ser a su vez otro polinomio) Entenderemos que 
    dos polinomios son iguales si las tuplas son iguales y el objeto var de cada polinomio es el mismo.
    
    Las operaciones por el momento están definidas cuando el objeto es el mismo en los dos operandos o 
    cuando uno de los operandos es un escalar. ''' 
    
    def __new__(cls,a,v):
        return  tuple.__new__(cls,a,v)
    deg = property(operator.meth   odcaller('__len__'))
    var = property(operator.itemgetter(-1))
    l = property(operator.itemgetter(0))
    
    def __add__(self,a):
        if type(a) == mod and a.m == self.m:
            return mod(self.a+a.a,self.m)
        else:
            return mod(self.a+a,self.m)
    def __radd__(self,a):
        return mod(self.a+a,self.m)
    def __sub__(self,a):
        if type(a) == mod and a.m == self.m:
            return mod(self.a-a.a,self.m)
        else:
            return mod(self.a-a,self.m)
    def __rsub__(self,a):
        return mod(a-self.a,self.m)
    def __mul__(self,a):
        if type(a) == mod and a.m == self.m:
            return mod(self.a*a.a,self.m)
        else:
            return mod(self.a*a,self.m)
    def __rmul__(self,a):
        return mod(self.a*a,self.m)
    def __truediv__(self,a):
        if type(a) == mod and a.m == self.m:
            return self*a.__invert__()
        else:
            return mod(self.a*(~mod(a,self.m)).a,self.m)
    def __rtruediv__(self,a):
        return mod(a*self.__invert__().a,self.m)
    def __pow__(self,a):
        if type(a)== int :
            return generic_pow(self,a)
        return NotImplemented
    def __invert__(self):
        u,v,d = bezout(self.a,self.m)
        if abs(d) == type(self.a)(1):
            return mod(d*u,self.m)
    def __int__(self):
        if a.type == int:
            return self.a
        return NotImplemented
    def __repr__(self):
        return 'mod(' + str(self.a)+','+str(self.m)+')' 
