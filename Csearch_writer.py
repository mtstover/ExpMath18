#Execute from command line as
#python 'Csearch_writer.py' n mu1 d
#in which n is replaced by the desired number of terms
#and mu1,d are replaced by the eigenvalue and the degree

#returns lower bound for C(mu1,d)
#and spherical coordinates of Chebyshev coefficients for bound polynomial

import sys
args=sys.argv
n=int(args[1])
mu1=args[2]
d=args[3]

f=open('searchprogram.py','w')
f.write('from numpy import *\n')
f.write('import copy\n')

#Maxfunc
f.write('def maxfunc(g,a,b,n):\n'+'    if (g(a)>0 or g(b)>0):\n'+'        return 0\n'+'    N=int(5*(b-a))\n'+'    M=4\n'+'    iternum=3\n'+'    deltax=(b-a)/N\n'+'    Dnext=g(a+deltax/1000)-g(a)\n'+'    index=[]\n'+'    for i in range(1,N+1):\n'+'        Dprev=Dnext\n'+'        Dnext=g(a+deltax*(i+1.0/1000.0))-g(a+deltax*(i-1.0/1000.0))\n'+'        val=a+(i-1)*deltax\n'+'        if (g(val)>0):\n'+'            return 0\n'+'        if (Dprev >= 0 and Dnext < 0):\n'+'            index.append(val)\n'+'    count=0\n'+'    while (count < iternum):\n'+'        deltax=deltax/M\n'+'        altindex=[]\n'+'        for j in index:\n'+'            Dnext=g(j+deltax/1000)-g(j)\n'+'            for i in range(1,M+1):\n'+'                Dprev=Dnext\n'+'                Dnext=g(j+deltax*(i+1.0/1000.0))-g(j+deltax*(i-1.0/1000.0))\n'+'                if (Dnext < 0):\n'+'                    altindex.append(j+(i-1)*deltax)\n'+'                    break\n'+'        index=altindex\n'+'        count +=1\n'+'    value=max(g(a),g(b))\n'+'    for i in index:\n'+'        if (g(i)> value):\n'+'            value=g(i)\n'+'    return value\n')
    
#Basis Polynomials
f.write('def V0(x):\n')
f.write('   return(1)\n')
f.write('def V1(x):\n')
f.write('   return(x)\n')

for i in range(n-1):
    f.write('def V'+str(i+2)+'(x):\n')
    f.write('   return(x*V'+str(i+1)+'(x)-V'+str(i)+'(x))\n')

#Coordinate TF
f.write('def polyfloat(smaxarray):\n')
for i in range(n-1):
    f.write('   theta'+str(i)+'=smaxarray['+str(i)+']\n')
f.write('   core0=sin(theta0)\n')
for i in range(n-3):
    f.write('   core'+str(i+1)+'=core'+str(i)+'*sin(theta'+str(i+1)+')\n')
f.write('   c1=float(cos(theta0))\n')
for i in range(n-2):
    f.write('   c'+str(i+2)+'=float(cos(theta'+str(i+1)+')*core'+str(i)+')\n')
f.write('   c'+str(n)+'=float(sin(theta'+str(n-2)+')*core'+str(n-3)+')\n')
f.write('   return([c1')
for i in range(n-1):
    f.write(',c'+str(i+2))
f.write('])\n')

#Calculate C
f.write('def C(theta0')
for i in range(n-2):
    f.write(',theta'+str(i+1))
f.write(',L,epsilon):\n')
f.write('   coeff=polyfloat([theta0')
for i in range(n-2):
    f.write(',theta'+str(i+1))
f.write('])\n')
for i in range(n):
    f.write('   c'+str(i+1)+'=coeff['+str(i)+']\n')
f.write('   def F(x):\n')
f.write('      return(c1*V1(x)')
for i in range(n-1):
    f.write('+c'+str(i+2)+'*V'+str(i+2)+'(x)')
f.write(')\n')
f.write('   if(F(2-epsilon)<0 and F(-L)<F(2-epsilon) and F(0)<0):\n')
f.write('       M=maxfunc(F,-L,2-epsilon,'+str(n)+')\n')
f.write('   else:\n')
f.write('       M=0\n')
f.write('   Cbd=0\n')
f.write('   if(M-F(L)!=0):\n')
f.write('       Cbd=M/(M-F(L))\n')
f.write('   return(Cbd)\n')

#Smax function
f.write('def smax(theta01,theta02')
for i in range(n-2):
    f.write(',theta'+str(i+1)+'1')
    f.write(',theta'+str(i+1)+'2')
f.write(',grain,L,epsilon):\n')
f.write('   cm=min(0,(2-epsilon)/(2-epsilon-L))\n')
f.write('   anglearraym=[0')
for i in range(n-2):
    f.write(',0')
f.write(']\n')
f.write('   anglearray=copy.deepcopy(anglearraym)\n')
def indentgen(indent):
    h=''
    for i in range(indent):
        h=h+'   '
    return(h)
for i in range(n-1):
    f.write(indentgen(i)+'   for j'+str(i)+' in range(grain+1):\n')
    f.write(indentgen(i+1)+'   anglearray['+str(i)+']=theta'+str(i)+'1+(theta'+str(i)+'2-theta'+str(i)+'1)*j'+str(i)+'/grain\n')
f.write(indentgen(n)+'c=C(anglearray[0]')
for i in range(n-2):
    f.write(',anglearray['+str(i+1)+']')
f.write(',L,epsilon)\n')
f.write(indentgen(n)+'if(cm<c):\n')
f.write(indentgen(n+1)+'cm=c\n')
f.write(indentgen(n+1)+'anglearraym=copy.deepcopy(anglearray)\n')
f.write('   anglearraym.append(cm)\n')
f.write('   return(anglearraym)\n')

#Bound search function

#This data defines search parameters for each number of terms.
#Parameters have been optimized for fastest convergence,
#but the n=3,4,5 cases have been modified to increase precision for
#the calculations used in our paper.
paramslist=[[0.4604064686093708, 20, 34], [0.3062091953808147, 8, 10], [0.18, 16, 20], [0.18781219173304786, 15, 18], [0.16335254427484058, 5, 7], [0.13610501554746368, 5, 12], [0.11970051327385134, 5, 13], [0.10684329854240615, 5, 15], [0.09649044434230802, 5, 17], [0.08797279621411469, 5, 18], [0.08084078384200599, 5, 20], [0.0747809484060451, 5, 22], [0.06956794402800048, 5, 23], [0.0650355523277687, 5, 25], [0.06105845203373528, 5, 27], [0.05754035183957318, 5, 28], [0.0544060342012993, 5, 30]]

f.write('def msearch(mu1,d):\n')
f.write('   L=float(d/sqrt(d-1))\n')
f.write('   epsilon=float(2-mu1/sqrt(d-1))\n')
f.write('   ps='+str(paramslist[n-3])+'\n')
f.write('   k=ps[0]\n   grain=ps[1]\n   loops=int(ps[2])\n   print(str(loops)+\' loops\')\n   print(\'grain \'+str(grain))\n')
f.write('   pti=[pi/4')
for i in range(n-2):
    f.write(',pi/4')
f.write(']\n')
f.write('   baser=pi/4\n   fforward=False\n')
f.write('   searchresult=[0')
for i in range(n-1):
    f.write(',0')
f.write(']\n')
f.write('   for i in range(loops):\n')
f.write('       if(fforward==False):\n')
f.write('           print(\'running loop\'+str(i+1))\n')
f.write('           windr=baser*exp(-k*i)\n')
for i in range(n-1):
    f.write('           t'+str(i)+'1=max(pti['+str(i)+']-windr,0)\n')
    f.write('           t'+str(i)+'2=max(pti['+str(i)+']+windr,pi/2)\n')
f.write('           check=smax(t01,t02')
for i in range(n-2):
    f.write(',t'+str(i+1)+'1,t'+str(i+1)+'2')
f.write(',grain,L,epsilon)\n')
f.write('           if(check['+str(n-1)+']>searchresult['+str(n-1)+']):\n')
f.write('               searchresult=copy.deepcopy(check)\n')
f.write('               del check[-1]\n')
f.write('               pti=copy.deepcopy(check)\n')
f.write('           print(\'loop \'+str(i+1)+\'; arc error +/- \'+str(windr))\n'+'           print(\'current maximum bound:\')\n'+'           print(searchresult)\n'+'           print(\'___________\')\n')
f.write('   print(\'search complete:\')\n')
f.write('   print(\'coefficients:\')\n')
f.write('   print(polyfloat(pti))\n')
f.write('   print(\'searchresult:\')\n')
#f.write('   print(searchresult)\n')
f.write('   return searchresult\n')

f.write('res=msearch('+mu1+','+d+')\n')
f.write('print(res)\n')
f.write('print(\' C lower bound: \')\n')
f.write('print res['+str(n-1)+']\n')

f.close()
execfile('searchprogram.py')
