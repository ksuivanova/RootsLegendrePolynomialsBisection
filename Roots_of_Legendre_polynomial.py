from pylab import *
import numpy as np
import math

def bisec(a,b,f,tol): #Bisection method

    fa = f(a)
    m=(a+b)/2
    while abs(a-b)>=tol:
        m=(a+b)/2
        fm = f(m)
        if fm == 0:
            return m
        if sign(fm*fa)<0:
            b = m
        else:
            a=m
            fa = fm #this condition is actually useless, we can delete it, bc in this particular case sign(f(a))=sign(f(m)) 
    return m

def horner(poly, x):
    
        n = len(poly)    
        result = poly[0]
        for i in range(1,n):
            result = result*x+poly[i]
        return result 

   
def sign_count(L):  #Counts number of sign changes in a list L
    
    c=0
    s=sign(L[0])
    for i in L[1:]:
        if s*i<0:
        	c=c+1
        	s=sign(i)
    return c

def legendre_polynomial(n,x):  #Definition of n'th legendre polynomial  as recursive function, here it is stored in a list

    for k in range(n):
        L= [1,x]
        for j in range(1,n):
        	t=x*L[j]-(j)**2/float(4*(j)**2-1)*L[j-1]
        	L.append(t)
    return L

def LnHorner(n):  #Definition of n'th legendre polynomial  using Horner's method
    oldzeros=[]
    zeros=[]
    a = -1
    b = 1
    tol = 10**(-7)
    if (n==0):
        poly = [1]
    elif(n==1):
        poly = [1,0]
        zeros=[0]
        
    else:   
        poly = [1, 0]
        zeros=[0]      
        polyminus1 = [1]
        polyminus2 = []
       
        for i in range(2,n+1):
            polyminus2 = polyminus1
            polyminus1 = poly
            polyminus1.extend([0])
            polyminus2[:0]=[0,0]  
            poly= [x-( (i-1)**2/(4*(i-1)**2-1) )*y for x, y in zip(polyminus1, polyminus2)]
            zeros.append(b)
            oldzeros.extend(zeros)
            oldzeros.insert(0,a)
            Leg = lambda x: horner(poly, x)
            for j in range(i):
                zeros[j] = bisec(oldzeros[j],oldzeros[j+1],Leg,tol)
    return poly, zeros

   
def leg_root1(i,n):    #i'th  root of n'th legendre polynomial using the sign changes
    a = -1
    b = 1
    tol = 10**(-7)
    while abs(b-a)>  tol:
        x = (a+b)/2.
        L=legendre_polynomial(n,x)
        
        if sign_count(L) > i-1:
            a = x
        else:
            b = x
    return x

#main
if __name__ == '__main__':
    
    print('using the sign changes',leg_root1(8,10))
    poly, zeros = LnHorner(10)
    print('using Bisection',zeros[8])
   