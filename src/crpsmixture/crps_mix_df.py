import numpy as np
from _crpsmixture import crpsmixGw
#import cpp_int_lims
from scipy import stats
from scipy import integrate




#def crps_tdis_lims(preds, y, h, df):
#    if hasattr(y,  "__len__") == False:
#        y = np.array([y])
#    if len(preds.predictions) != len(y):
#        raise ValueError("preds same length as y")
#    if h < 0:
#        raise ValueError("h must be positive")
#    if df < 0:
#        raise ValueError("df must be positive")
#    #if hasattr(y,  "__len__"):
#    len_preds = [len(x.points) for x in preds.predictions]
#    len_cumsum = np.cumsum(len_preds)
#    len_cumsum = np.insert(len_cumsum, 0, 0)
#    mean = np.concatenate([x.points for x in preds.predictions])
#    weights = np.concatenate([np.diff(np.insert(x.ecdf, 0, 0)) for x in preds.predictions])
#    crps = cpp_int_lims.cpp_int_lims(y,mean,weights, len_cumsum, h, df, -1*float('Inf'), float('Inf'))
#    return crps

def crps_int_t(y, w, m, h, df, rel_tol):
    #ind = (y > upper) - (y < lower)
    F1 = lambda x: np.sum(w*stats.t.cdf(x, df = df, loc = m, scale = h))**2
    F2 = lambda x: (1-np.sum(w*stats.t.cdf(x, df = df, loc = m, scale = h)))**2
    s1, err1 = integrate.quad(F1, -np.inf, y, epsabs = rel_tol)
    s2, err2 = integrate.quad(F2, y, np.inf, epsabs = rel_tol)
    #s2 = 0
    return s1 + s2

def crps_t_int(preds, y, h, df, rel_tol = 1e-6):
    n = len(y)
    crps_val = np.zeros(n)
    def prep_int(predictions, y):
        w = np.diff(np.insert(predictions.ecdf, 0, 0)) 
        m = predictions.points
        return crps_int_t(y, w, m, h, df, rel_tol)  
    A = list(map(prep_int, preds.predictions, y))
    return np.mean(A)

def smooth_crps(preds, y, h, df=None):
    if df == None:
        n = len(y)
        crps_val = np.zeros(n)
        for i in range(n):
            m = preds.predictions[i].points
            w = np.diff(np.insert(preds.predictions[i].ecdf, 0, 0))
            crps_val[i] = crpsmixGw(m, w, y[i], h)
        return np.mean(crps_val)
    else:
        #return crps_tdis_lims(preds, y, h, df)
        return crps_t_int(preds, y, h, df)
