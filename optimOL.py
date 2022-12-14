#!/usr/bin/env python
# coding: utf-8

# In[1]:


# optimOL.py 
# see README.txt file for details.
#----------------------------------
import numpy as np
from scipy.optimize import minimize


# In[2]:


#--------------------------------------------
# Necesary functions
#-----------------------------
def readOL(dir):
# Read all the files required for the problem.
# All files must be located in directory: dir
# Files to read:
# namelistOL.txt, file-cf1.txt, file-ref.txt,file-cf2.txt,file-min.txt
# See README.txt for details.
#-------------------------------------------------
# Read file namelistOL.txt.
# Lines starting with # are comments.
# In the first  line: Ns,Nw,rs2w,Isc,Iwc,Iscmin,Iscmax,Iwcmin,Iwcmax
# In the second line: IscminR (vector of length Ns)
# In the third  line: IscmaxR (vector of length Ns)
# In the fourth line: IwcminR (vector of length Nw)
# In the fifth  line: IwcmaxR (vector of length Nw)
    global Ns,Nw,rs2w,Isc,Iwc,Iscmin,Iscmax,Iwcmin,Iwcmax
    global IscminR,IscmaxR,IwcminR,IwcmaxR
    file_param=dir+'/namelistOL.txt'
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
                    Isc=float((t[3].strip()))
                    Iwc=float((t[4].strip()))
                    Iscmin=float((t[5].strip()))
                    Iscmax=float((t[6].strip()))
                    Iwcmin=float((t[7].strip()))
                    Iwcmax=float((t[8].strip()))
                if(icount==1):
                    IscminR=[]
                    for i in range(len(t)):
                        IscminR=np.append(IscminR,float(t[i].strip()))
                if(icount==2):
                    IscmaxR=[]
                    for i in range(len(t)):
                        IscmaxR=np.append(IscmaxR,float(t[i].strip()))
                if(icount==3):
                    IwcminR=[]
                    for i in range(len(t)):
                        IwcminR=np.append(IwcminR,float(t[i].strip()))
                if(icount==4):
                    IwcmaxR=[]
                    for i in range(len(t)):
                        IwcmaxR=np.append(IwcmaxR,float(t[i].strip()))
# put Iscmin/max and Iwcmin/max into a single vector
    global IminR,ImaxR
    IminR=np.append(IscminR,IwcminR)
    ImaxR=np.append(IscmaxR,IwcmaxR)
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
# CC : matrix of shape NMPx(Ns+Nw)
# first NS columns correspond to the matrix As and last Nw columns corrspond to Aw.
# each row correspods to one time interval
# See README.txt for details.
    global NMP,CC
    CC=np.loadtxt(dir+'/file-cf2.txt')
    NMP=CC.shape[0] # number of row
    l1=CC.shape[1]  # number of columns
    if(l1 != Ns+Nw):
        print('wrong dimensions in file-cf2.txt')
# read the file: file-min.txt. One column and NMP rows
# MM : vector of length NMP. 
# See README.txt for details.
    global MM
    MM=np.loadtxt(dir+'/file-min.txt')
    if(MM.shape[0] != NMP):
        print('wrong dimensions in file-min.txt')
    return

def P_anom(I):  # Total production anomalies.
    p=np.dot(AA,I)-BB
    return p
   
def fun1(I): # Defines the minimizing function for OL and OLS problems, Eq. (7) in README.txt.
    p=P_anom(I)
    # f=np.var(p) not the same if the mean is not zero.
    f=sum(p*p/(len(p)))
    return f 

def fconstr1(I):
# defines the constrain of minimun capacity of production. Eq. (8) in README.txt.
    Prod=np.dot(CC,I) # produccion
    fcons=0
    for i in range(NMP):
        if(Prod[i] < MM[i]):
            fcons=fcons+(Prod[i]-MM[i])**2
    return fcons

def fconstr2(I):
# defines the constrain of minimun/maximun threshold per sub-region. Eq. (9) in README.txt.
    fcons=0
    for i in range(len(IminR)):
        if(I[i] < IminR[i]):
            fcons=fcons+(I[i]-IminR[i])**2
        if(I[i] > ImaxR[i]):
            fcons=fcons+(I[i]-ImaxR[i])**2
    return fcons

def fconstr3(I):
# defines the constrain of total shares of each technology. Eq. (10) in README.txt.
    fcons=(sum(I[0:Ns])-Isc)**2+(sum(I[Ns:Ns+Nw])-Iwc)**2
    return fcons

def fconstr4(I):    
# defines the  constrain de ratio solar to wind technology: rs2w, eq. (11).
    Is=sum(I[:Ns])
    Iw=sum(I[Ns:])
    fcons=0
    if(rs2w>1 and Is<Iw):
        fcons=fcons+(Is-Iw)**2
    elif(0<rs2w and rs2w<1 and Is>Iw):
        fcons=fcons+(Is-Iw)**2
    return fcons

def fconstr5(I):    
#defines the constrain of total ammount of installed capacity, eq. (12).
    Is=sum(I[:Ns])
    Iw=sum(I[Ns:])
    fcons=0
    if(Is<Iscmin):
        fcons=fcons+(Is-Iscmin)**2
    elif(Is>Iscmax):
        fcons=fcons+(Is-Iscmax)**2
    if(Iw<Iwcmin):
        fcons=fcons+(Iw-Iwcmin)**2
    elif(Iw>Iwcmax):
        fcons=fcons+(Iw-Iwcmax)**2
    return fcons

def funOL(X,C1,C2,C3):
# defines function to minimize for the OD problem.
    I=X*X
    f0=fun1(I)
    fcons1=fconstr1(I)
    fcons2=fconstr2(I)
    fcons3=fconstr3(I)
    f=f0+C1*fcons1+C2*fcons2+C3*fcons3
#    print('f:',f0,fcons1,fcons2,fcons3)
    return f

def funODL(X,C1,C2,C4,C5):
# defines function to minimize for the OD problem.
    I=X*X
    f0=fun1(I)
    fcons1=fconstr1(I)
    fcons2=fconstr2(I)
    fcons4=fconstr4(I)
    fcons5=fconstr5(I)
    f=f0+C1*fcons1+C2*fcons2+C4*fcons4+C5*fcons5
#    print('f:',f0,fcons1,fcons2,fcons4,foncs5)
    return f

def funtryOL(X):
# used to try to verify possible solutions within the constrains
    I=X*X
    fcons1=fconstr1(I)
    fcons2=fconstr2(I)
    fcons3=fconstr3(I)
    f=fcons1+fcons2+fcons3
#    print(f0,fcons1,fcons2,fcons3)
    return f
def funtryOLS(X):
# used to try to verify possible solutions within the constrains
    I=X*X
    fcons1=fconstr1(I)
    fcons2=fconstr2(I)
    fcons4=fconstr4(I)
    fcons5=fconstr5(I)
    f=fcons1+fcons2+fcons4+fcons5
#    print(f0,fcons1,fcons2,fcons4,fcons5)
    return f


# In[3]:


# Read all data of the OL or OLS problem.
dir='.'
readOL(dir)
print('Ns=',Ns)
print('Nw=',Nw)
print('rs2w=',rs2w)
print('Isc=',Isc)
print('Iwc=',Iwc)
print('Iscmin=',Iscmin)
print('Iscmax=',Iscmax)
print('Iwcmin=',Iwcmin)
print('Iwcmax=',Iwcmax)
print('IscminR:',IscminR)
print('IscmaxR:',IscmaxR)
print('IwcminR:',IwcminR)
print('IwcmaxR:',IwcmaxR)
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
    X[i]=np.random.uniform(np.sqrt(max(IminR[i],0)),np.sqrt(min(ImaxR[i],1)))
res=minimize(funtryOL,X,tol=tol)
X=res.x
ff=funtryOL(X)
if(ff > 1.e-11):
    print('It seems that constrains are not compatible')
    print('____________________________________________')
    print('Fisrt Constrain:')
    I=X*X
    Prod=np.dot(CC,I)
    for i in range(NMP):
        if(Prod[i]<MM[i]):
            print(i,Prod[i],'<',MM[i])
    print()
    print('Second Constrain:')
    for i in range(len(IminR)):
        if(i<Ns):
            Ty='Solar'
            ii=i+1
        else:
            Ty='Wind'
            ii=i+1-Ns
        if(I[i] < IminR[i]):
            print(Ty+' region:',ii,'   I[i]=%6.3f < Imin[i]=%6.3f' %(I[i],IminR[i]))
        if(I[i] > ImaxR[i]):
            print(Ty+' region:',ii,'   I[i]=%6.3f > Imax[i]=%6.3f' %(I[i],ImaxR[i])) 
    print()
    print('Third Constrain:')
    print('Sum(Is[i])=%8.3f  Isc=%8.3f' %(sum(I[:Ns]),Isc))
    print('Sum(Iw[i])=%8.3f  Iwc=%8.3f' %(sum(I[Ns:]),Iwc))


# In[4]:


#Optimize OL problem
tol=1.e-6
C1=10.0
C2=100.0
C3=100.0
fc1=100
fc2=100
fc3=5000
while(fc1>tol or fc2>tol or fc3>tol):
#    print('C:',C1,C2,C3)
#    print('f:',fc1,fc2,fc3)
    if(fc1>tol):
        C1=2.0*C1
    if(fc2>tol):
        C2=2.0*C2
    if(fc3>tol):
        C3=2.0*C3
    res=minimize(funOL,X,args=(C1,C2,C3))
    X=res.x
    I=X*X
    fc1=fconstr1(I)
    fc2=fconstr2(I)
    fc3=fconstr3(I)
f0=fun1(I)
print('fmin=%5.3f  fcons1=%9.3e  fcons1=%9.3e fcons3=%9.3e' % (f0,fc1,fc2,fc3))
Is=sum(I[:Ns])
Iw=sum(I[Ns:])
print('Isc=%5.3f   Sum(Is)=%5.3f' %(Isc,Is))
print('Iwc=%5.3f   Sum(Iw)=%5.3f' %(Iwc,Iw))
print('_____________')
for i in range(Ns):  
    print(' Is(%2d)= %6.3f' %(i+1,I[i]))
print('_____________')
for i in range(Nw):   
    print(' Iw(%2d)= %6.3f' %(i+1,I[i+Ns])) 
print('_____________')
Prod=np.dot(CC,I) # produccion
for i in range(NMP):
    print(' Prod[%2d]=%9.3f, MM[%2d]=%9.3f' %(i+1,Prod[i],i+1,MM[i]))


# In[51]:




