#!/usr/bin/env python
# coding: utf-8

# optimOD.py 
# see README.txt file for details.
#----------------------------------
import numpy as np
from scipy.optimize import minimize

#--------------------------------------------
# Necesary functions
#-----------------------------
def readOD(dir):
# Read all the files required for the problem.
# All files must be located in directory: dir
# Files to read:
# namelistOD.txt, file-cf1.txt, file-ref.txt,file-cf2.txt,file-min.txt
# See README.txt for details.
#-------------------------------------------------
# Read file namelistOD.txt.
# Lines starting with # are comments.
# In the first  line: Ns,Nw,rs2w,Ssc,Swc,Sscmin,Sscmax,Swcmin,Swcmax
# In the second line: SscminR (vector of length Ns)
# In the third  line: SscmaxR (vector of length Ns)
# In the fourth line: SwcminR (vector of length Nw)
# In the fifth  line: SwcmaxR (vector of length Nw)
    global Ns,Nw,rs2w,Ssc,Swc,Sscmin,Sscmax,Swcmin,Swcmax
    global SscminR,SscmaxR,SwcminR,SwcmaxR
    file_param=dir+'/namelistOD.txt'
    with open(file_param,'r') as file:
        param={}
        icount=-1
        for line in file.readlines():
            t=line.split()
            icount=icount+1
            #print(icount,len(t))
            if(len(t)==0 or t[0].strip()[0]=='#'):
                icount=icount-1 
            else:
                if(icount==0):
                    Ns=int(t[0].strip())
                    Nw=int(t[1].strip())
                    rs2w=float((t[2].strip()))
                    Ssc=float((t[3].strip()))
                    Swc=float((t[4].strip()))
                    Sscmin=float((t[5].strip()))
                    Sscmax=float((t[6].strip()))
                    Swcmin=float((t[7].strip()))
                    Swcmax=float((t[8].strip()))
                if(icount==1):
                    SscminR=[]
                    for i in range(len(t)):
                        SscminR=np.append(float(t[i].strip()),SscminR)
                if(icount==2):
                    SscmaxR=[]
                    for i in range(len(t)):
                        SscmaxR=np.append(float(t[i].strip()),SscmaxR)
                if(icount==3):
                    SwcminR=[]
                    for i in range(len(t)):
                        SwcminR=np.append(float(t[i].strip()),SwcminR)
                if(icount==4):
                    SwcmaxR=[]
                    for i in range(len(t)):
                        SwcmaxR=np.append(float(t[i].strip()),SwcmaxR)
    if(Ssc <0 or Swc <0 or Ssc>1 or Swc>1 or abs(Ssc+Swc-1)> 0.01):
        print('Check Ssc or Swc values')
    if(Sscmin>1 or Sscmax<=0 or Swcmin>1 or Swcmax<=0):
        print('Check Sscmin, Sscmax, Swcmin, Swcmax')
# put Sscmin/max and Swcmin/max into a single vector
    global SminR,SmaxR
    SminR=np.append(SscminR,SwcminR)
    SmaxR=np.append(SscmaxR,SwcmaxR)
# guarantee that Ssc and Swc add upto 1. 
    tot=Ssc+Swc
    Ssc=Ssc/tot  
    Swc=Swc/tot
# read the file: file-cf1.txt
# AA : matrix of shape NTTx(Ns+Nw)
# first NS columns correspond to the matrix As and last Nw columns corrspond to Aw.
# See README.txt for details.
    global NTT,AA
    AA=np.loadtxt(dir+'/file-cf1.txt')
    NTT=AA.shape[0] # number of row
    l1=AA.shape[1]  # number of columns
    if(l1 != Ns+Nw):
        print('wrong dimensions in file-cf1.txt')
# read the file: file-ref.txt. One column and NTT rows.
# BB : vector of length NTT. 
# See README.txt for details.
    global BB
    BB=np.loadtxt(dir+'/file-ref.txt')
    if(BB.shape[0] != NTT):
        print('wrong dimensions in file-ref.txt')
# read the file: file-cf2.txt
# CC : matrix of shape 12x(Ns+Nw)
# first NS columns correspond to the matrix As and last Nw columns corrspond to Aw.
# each row correspods to one month
# See README.txt for details.
    global NMP,CC
    CC=np.loadtxt(dir+'/file-cf2.txt')
    NMP=CC.shape[0] # number of row
    l1=CC.shape[1]  #number of columns
    if(l1 != Ns+Nw):
        print('wrong dimensions in file-cf2.txt')
# read the file: file-min.txt. One column and NMP rows
# MM : vector of length NMP. 
# See README.txt for details.
    global MM
    MM=np.loadtxt(dir+'/file-min.txt')
    if(MM.shape[0] != NMP):
        print('wrong dimensions in file-min.txt')
    for j in range(NMP):
        ccmax=max(CC[j,:])
        if(ccmax < MM[j]):
            print('wrong construction of file-min.txt')
    return

def P_anom(S):  # Total production anomalies.
    p=np.dot(AA,S)-BB
#    print('p:',p.shape,AA.shape,S.shape)
    return p
   
def fun1(S): # Defines the minimizing function. Eq. (1) in README.txt.
    p=P_anom(S)
    # f=np.var(p) not the same if the mean is not zero.
    f=sum(p*p/(len(p)))
    return f 

def fconstr1(S):
# defines the condition of minimun capacity of production. Eq. (2) in README.txt.
# See section 1.4 of READM.txt for explanation.
    Prod=np.dot(CC,S) # produccion
    fcons=0
    for i in range(NMP):
        if(Prod[i] < MM[i]):
            fcons=fcons+(Prod[i]-MM[i])**2
    return fcons

def fconstr2(S):
# defines el constrain of minimun/maximun threshold per sub-region. Eq. (3) in README.txt.
    fcons=0
    for i in range(len(SminR)):
        if(S[i] < SminR[i]):
            fcons=fcons+(S[i]-SminR[i])**2
        if(S[i] > SmaxR[i]):
            fcons=fcons+(S[i]-SmaxR[i])**2
    return fcons

def fconstr3(S):
# defines el constrain total shares of each technology. Eq. (4) in README.txt.
    fcons=(sum(S[0:Ns])-Ssc)**2+(sum(S[Ns:Ns+Nw])-Swc)**2
    return fcons

def fconstr4(S):    
#define el constrain de ratio solar to wind technology: rs2w, eq. (5).
    Ss=sum(S[:Ns])
    Sw=sum(S[Ns:])
    fcons=0
    if(rs2w>1 and Ss<Sw):
        fcons=fcons+(Ss-Sw)**2
    elif(0<rs2w and rs2w<1 and Ss>Sw):
        fcons=fcons+(Ss-Sw)**2
    return fcons

def fconstr5(S):    
#define el constrain de ratio solar to wind technology: rs2w, eq. (6).
    Ss=sum(S[:Ns])
    Sw=sum(S[Ns:])
    fcons=0
    if(Ss<Sscmin):
        fcons=fcons+(Ss-Sscmin)**2
    elif(Ss>Sscmax):
        fcons=fcons+(Ss-Sscmax)**2
    if(Sw<Swcmin):
        fcons=fcons+(Sw-Swcmin)**2
    elif(Sw>Swcmax):
        fcons=fcons+(Sw-Swcmax)**2
    return fcons

def funOD(X,C1,C2,C3):
# defines function to minimize for the OD problem.
    Ctot=np.dot(X,X)
    S=X*X/Ctot
    f0=fun1(S)
    fcons1=fconstr1(S)
    fcons2=fconstr2(S)
    fcons3=fconstr3(S)
    f=f0+C1*fcons1+C2*fcons2+C3*fcons3
#    print('f:',f0,fcons1,fcons2,fcons3)
    return f

def funODS(X,C1,C2,C4,C5):
# defines function to minimize for the OD problem.
    Ctot=np.dot(X,X)
    S=X*X/Ctot
    f0=fun1(S)
    fcons1=fconstr1(S)
    fcons2=fconstr2(S)
    fcons4=fconstr4(S)
    fcons5=fconstr5(S)
    f=f0+C1*fcons1+C2*fcons2+C4*fcons4+C5*fcons5
#    print('f:',f0,fcons1,fcons2,fcons4)
    return f

def funtryOD(X):
# used to try to verify possible solutions within the constrains
    Ctot=np.dot(X,X)
    S=X*X/Ctot
    fcons1=fconstr1(S)
    fcons2=fconstr2(S)
    fcons3=fconstr3(S)
    f=fcons1+fcons2+fcons3
#    print(f0,fcons1,fcons2,fcons3)
    return f
def funtryODS(X):
# used to try to verify possible solutions within the constrains
    Ctot=np.dot(X,X)
    S=X*X/Ctot
    fcons1=fconstr1(S)
    fcons2=fconstr2(S)
    fcons4=fconstr4(S)
    fcons5=fconstr5(S)
    f=fcons1+fcons2+fcons4+fcons5
#    print(f0,fcons1,fcons2,fcons4,fcons5)
    return f

# Read all data of the OD or ODS problem.
dir='.'
readOD(dir)
#------------------------
print('Ns=',Ns)
print('Nw=',Nw)
print('rs2w=',rs2w)
print('Ssc=',Ssc)
print('Swc=',Swc)
print('Sscmin=',Sscmin)
print('Sscmax=',Sscmax)
print('Swcmin=',Swcmin)
print('Swcmax=',Swcmax)
print('SscminR:',SscminR)
print('SscmaxR:',SscmaxR)
print('SwcminR:',SwcminR)
print('SwcmaxR:',SwcmaxR)
print('NTT=',NTT)
print('NMP=',NMP)
#------------------------------------------------------------------------------------
# Check constrains
# inicializa X
X=np.random.rand(Ns+Nw)   #se podria leer de un fichero
# search for a possible solution of the constrains with tolerance=tol
# If a solution is found it is used as initial guess for the full minimization problem
tol=1.e-6   
for i in range(Ns+Nw):
    X[i]=np.random.uniform(np.sqrt(max(SminR[i],0)),np.sqrt(min(SmaxR[i],1)))
res=minimize(funtryOD,X,tol=tol)
X=res.x
ff=funtryOD(X)
if(ff > 1.e-11):
    print('It seems that constrains are not compatible')
    print('____________________________________________')
    print('Fisrt Constrain:')
    Ctot=np.dot(X,X)
    S=X*X/Ctot
    Prod=np.dot(CC,S)
    for i in range(NMP):
        if(Prod[i]<MM[i]):
            print(i,Prod[i],'<',MM[i])
    print()
    print('Second Constrain:')
    for i in range(len(SminR)):
        if(i<Ns):
            Ty='Solar'
            ii=i+1
        else:
            Ty='Wind'
            ii=i+1-Ns
        if(S[i] < SminR[i]):
            print(Ty+' region:',ii,'   S[i]=%6.3f < Smin[i]=%6.3f' %(S[i],SminR[i]))
        if(S[i] > SmaxR[i]):
            print(Ty+' region:',ii,'   S[i]=%6.3f > Smax[i]=%6.3f' %(S[i],SmaxR[i]))           
    print()
    print('Third Constrain:')
    print('Sum(Ss[i])=%8.3f  Ssc=%8.3f' %(sum(S[:Ns]),Ssc))
    print('Sum(Sw[i])=%8.3f  Swc=%8.3f' %(sum(S[Ns:]),Swc))

#Optimize OD problem
tol=1.e-7
C1=10.0
C2=100.0
C3=100.0
fc1=100
fc2=100
fc3=100
while(fc1>tol or fc2>tol or fc3>tol):
#    print('C:',C1,C2,C3)
#    print('f:',fc1,fc2,fc3)
    if(fc1>tol):
        C1=2.0*C1
    if(fc2>tol):
        C2=2.0*C2
    if(fc3>tol):
        C3=2.0*C3
    res=minimize(funOD,X,args=(C1,C2,C3))
    X=res.x
    Ctot=np.dot(X,X)
    S=X*X/Ctot
    fc1=fconstr1(S)
    fc2=fconstr2(S)
    fc3=fconstr3(S)
f0=fun1(S)
print('fmin=%5.3f  fcons1=%9.3e  fcons1=%9.3e fcons3=%9.3e' % (f0,fc1,fc2,fc3))
Ss=sum(S[:Ns])
Sw=sum(S[Ns:])
print('Ssc=%5.3f   Sum(S_s)=%5.3f' %(Ssc,Ss))
print('Swc=%5.3f   Sum(S_w)=%5.3f' %(Swc,Sw))
print('_____________')
for i in range(Ns):  
    print(' Ss(%2d)= %5.3f' %(i+1,S[i]))
print('_____________')
for i in range(Nw):   
    print(' Sw(%2d)= %5.3f' %(i+1,S[i+Ns])) 
print('_____________')
Prod=np.dot(CC,S) # produccion
for i in range(NMP):
    print(' Prod[%2d]=%8.3f, MM[%2d]=%8.3f' %(i+1,Prod[i],i+1,MM[i]))





