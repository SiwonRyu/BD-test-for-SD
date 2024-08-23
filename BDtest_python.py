import operator as op
import numpy as np
import pandas as pd
import math
import time
import xlrd

# Importing data part
data = pd.read_excel ('C:\data_temp\ex.xlsx') #for an earlier version of Excel, you may need to use the file extension of 'xls'
print(data.shape)
x1 = np.array(data)[:,0]
x2 = np.array(data)[:,1]
x3 = np.array(data)[:,2]
x4 = np.array(data)[:,3]

# Parameters
[alpha,n1,n2,ngrid,B,R,gen] = [0.05,100,100,40,100,1,2]

# Generate : when gen = 1, use following desing and when gen <> 1, real imported data(x1 and x2) will be used
if gen == 1 :
    # Design
    [mu1, sigma1, mu2, sigma2] = [0, 1, 2, 1]
    sample1 = mu1 + sigma1 * np.random.randn(n1,1,1,R)
    sample2 = mu2 + sigma2 * np.random.randn(n2,1,1,R)

    minx = min(sample1.min(3).min(0), sample2.min(3).min(0))[0]
    maxx = max(sample1.max(3).max(0), sample2.max(3).max(0))[0]
else :
    sample1 = x1[:,None,None,None] #This dimension adjustement is for fix dimension in following functions
    sample2 = x2[:,None,None,None]

    minx = min(sample1.max(0), sample2.max(0))[0][0]
    maxx = max(sample1.max(0), sample2.max(0))[0][0]
grid = np.linspace(minx,maxx, num=ngrid)
s = 1

# Test part function
def bdtest(sample1,sample2,grid,s):
    start_time = time.time()
    if gen > 1 :
        repeat = 1
    else :
        repeat = R

    def operator(x,grid,s):
        return (np.apply_along_axis(op.lt, 1, x, grid)*np.apply_along_axis(op.sub, 1, x, grid)**(s-1)/math.factorial(s-1))[:,:,0]

    def ecdf(x,grid,s):
        return operator(x,grid,s).mean(0)[None]

    def stat(sample1,sample2,grid,s):
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        return (n1*n2/(n1+n2))**(0.5) * (ecdf(sample1,grid,s) - ecdf(sample2,grid,s)).max(1)[:,None]

    def op_J(x,grid,s):
        n = x.shape[0]
        u = np.random.randn(n,1,B,repeat)
        y = (operator(x,grid,s)-ecdf(x,grid,s))
        return n**(0.5) * (y*u).mean(0)[None]

    def multiplier1(sample1,sample2,grid,s):
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        lbd = n2/(n1+n2)
        return (lbd**0.5 * op_J(sample1,grid,s) - (1-lbd)**0.5 * op_J(sample2,grid,s)).max(1)[:,None]

    def multiplier2(sample1,sample2,grid,s):
        return op_J(sample2,grid,s).max(1)[:,None]

    def bsample(x,nb):
        n = x.shape[0]
        index = np.random.choice(n,[n,1,nb])
        return x[index][:,:,:,0,0]

    def boot1(sample1,sample2,grid,s):
        n2 = sample2.shape[0]
        bsample2 = bsample(sample2,B)
        return (n2)**0.5*(ecdf(bsample2,grid,s)-ecdf(sample2,grid,s)).max(1)[:,None]

    def boot2(sample1,sample2,grid,s):
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        bpsample = bsample(np.concatenate([sample1, sample2], axis=0), B)
        bpsample1 = bpsample[:n1,:,:]
        bpsample2 = bpsample[n1+1:,:,:]
        return (n1*n2/(n1+n2))**0.5*(ecdf(bpsample1,grid,s) - ecdf(bpsample2,grid,s)).max(1)[:,None]

    def boot3(sample1,sample2,grid,s):
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        bsample1 = bsample(sample1, B)
        bsample2 = bsample(sample2, B)
        return (n1 * n2 / (n1 + n2)) ** 0.5 * ((ecdf(bsample1,grid,s) - ecdf(sample1,grid,s))
            - (ecdf(bsample2,grid,s) - ecdf(sample2,grid,s))).max(1)[:,None]


    bd = stat(sample1,sample2,grid,s)
    pval_m1 = (multiplier1(sample1,sample2,grid,s) > bd).mean(2)[:,:,None]
    pval_m2 = (multiplier2(sample1,sample2,grid,s) > bd).mean(2)[:,:,None]
    pval_b1 = (boot1(sample1,sample2,grid,s) > bd).mean(2)[:,:,None]
    pval_b2 = (boot2(sample1,sample2,grid,s) > bd).mean(2)[:,:,None]
    pval_b3 = (boot3(sample1,sample2,grid,s) > bd).mean(2)[:,:,None]

    if s == 1 :
        Torder = 'FSD'
    elif s == 2 :
        Torder = 'SSD'
    elif s == 3 :
        Torder = 'TSD'
    else :
        Torder = str(s) + 'th order SD'

    if R > 1 & gen == 1 :
        print('Barrett-Donald Test for Stochastic Dominance'
              '\n* H0 : sample1', Torder, 'sample2')
        print('* Design : \n\tsample1 ~ N(%3.1f' % mu1, ',%3.1f)' % sigma1,
                  '\n\tsample2 ~ N(%3.1f' % mu2, ',%3.1f)\n' % sigma2)
        print('\n* SD order \t\t\t = %6d' % s,
              '\n* # of bootstrap \t = %6d' % B,
              '\n* # of grid points \t = %6d' % ngrid,
              '\n* # of repetition \t = %6d\n' % R)
        print('* Rejection probabilities :')
        print('\t- Multiplier 1 \t = %5.4f' % (pval_m1 < alpha).mean(3))
        print('\t- Multiplier 2 \t = %5.4f' % (pval_m2 < alpha).mean(3))
        print('\t- Bootstrap 1  \t = %5.4f' % (pval_b1 < alpha).mean(3))
        print('\t- Bootstrap 2  \t = %5.4f' % (pval_b2 < alpha).mean(3))
        print('\t- Bootstrap 3  \t = %5.4f' % (pval_b3 < alpha).mean(3))

    else :
        print('Barrett-Donald Test for Stochastic Dominance'
              '\n* H0 : sample1', Torder, 'sample2')
        if gen == 1:
            print('* Design : \n\tsample1 ~ N(%3.1f' % mu1, ',%3.1f)' % sigma1,
                  '\n\tsample2 ~ N(%3.1f' % mu2, ',%3.1f)' % sigma2)
        print('\n* SD order \t\t\t = %6d' % s,
              '\n* # of bootstrap \t = %6d' % B,
              '\n* # of grid points \t = %6d\n' % ngrid)
        print('* BD statistic \t\t = %5.4f' % bd)
        print('* p-values :')
        print('\t- Multiplier 1 \t = %5.4f' % pval_m1)
        print('\t- Multiplier 2 \t = %5.4f' % pval_m2)
        print('\t- Bootstrap 1  \t = %5.4f' % pval_b1)
        print('\t- Bootstrap 2  \t = %5.4f' % pval_b2)
        print('\t- Bootstrap 3  \t = %5.4f' % pval_b3)
    et = time.time() - start_time
    print('\n* Time elapsed : %5.2f Sec' % et)

bdtest(sample1,sample2,grid,1)