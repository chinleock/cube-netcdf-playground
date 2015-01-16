import numpy as NUM
from numpy import random as RAND
from pymc import DiscreteUniform as DU
from pymc import Exponential as EXP
from pymc import Poisson as POISSON
from pymc import MCMC


def getData():
#     data = RAND.random(300) * 100.
    data = NUM.array([2.,5.,67.,34.,8.,234.,45.,26.,23.,-9999.,1.,6.,1.,65.,4.,1.,3.,-9999.,34.,6.,1.,54.,26.,126.,23.,8.,21.,-9999.,21.,373.])
#     naIndex = RAND.randint(0,30, 3)
#     data[naIndex] = None
#     data = data.reshape(10,5,6)
    x = NUM.array([ 4, 5, 4, 0, 1, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6,
                    3, 3, 5, 4, 5, 3, 1, 4, 4, 1, 5, 5, 3, 4, 2, 5,
                    2, 2, 3, 4, 2, 1, 3, None, 2, 1, 1, 1, 1, 3, 0, 0,
                    1, 0, 1, 1, 0, 0, 3, 1, 0, 3, 2, 2, 0, 1, 1, 1,
                    0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2,
                    3, 3, 1, None, 2, 1, 1, 1, 1, 2, 4, 2, 0, 0, 1, 4,
                    0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1])
    return x

# def getOverallMinMax(data):
#     d = data
#     return NUM.nanmin(d), NUM.nanmax(d)

def getEarlyLateMean():
    ##totalTime = data.shape[0]
    ##changeTimePoint = RAND.randint(1,totalTime)
    ##return NUM.nanmean(data[0:switchpoint].ravel()), NUM.nanmean(data[switchpoint:-1].ravel())
    earlyMean = EXP('early_mean', beta = 1., value = 1)
    lateMean = EXP('late-mean', beta = 1., value = 1)
    return earlyMean, lateMean
     

# def getRate(switchpoint, earlyMean, lateMean, length):
#     rate = NUM.empty(length)
#     rate[:switchpoint] = earlyMean
#     rate[switchpoint:] = lateMean
#     return rate
x = NUM.array([ 4, -9999., 4, 0, 1, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6,
                3, 3, 5, 4, 5, 3, 1, 4, 4, 1, 5, 5, 3, 4, 2, 5,
                2, 2, 3, 4, 2, 1, 3, -9999., 2, 1, 1, 1, 1, 3, 0, 0,
                1, 0, 1, 1, 0, 0, 3, 1, 0, 3, 2, 2, 0, 1, 1, 1,
                0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2,
                3, 3, 1, -9999., 2, 1, 1, 1, 1, 2, 4, 2, 0, 0, 1, 4,
                0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0])
shape = (11,2,5)
x = x.reshape(shape)

def missingDataModel():
    switchpoint = DU('switchpoint', lower=0, upper=110)
    early_mean = EXP('early_mean', beta=1., value=1)
    late_mean = EXP('late_mean', beta=1., value=1)
    
    def rate(s=switchpoint, e=early_mean, l=late_mean):
        ''' Concatenate Poisson means '''
        out = NUM.empty(shape)
        out[:s] = e
        out[s:] = l
        return out
    maskVal = NUM.ma.masked_values(x, value = -9999.)
    maskVal.set_fill_value(-0.0000001)


    model = POISSON('data_model', mu=rate(), value=maskVal, observed=True)
    return locals()

if __name__ == '__main__':
#     data = getData()
#     ##overallMin, overallMax = getOverallMinMax(data)
#     switchpoint = DU('switchpoint', lower = 0, upper = 110)
#     ##switchpoint = DU('switchpoint', lower = overallMin, upper = overallMax)
#     earlyMean, lateMean = getEarlyLateMean()
#     rate = getRate(switchpoint, earlyMean, lateMean, len(data))
#     maskVal = NUM.ma.masked_values(data, value = -9999.)
#     model = POISSON('data_model', mu = rate, value = maskVal, observed = True)
    missingMCMC = MCMC(missingDataModel())
    missingMCMC.sample(100)