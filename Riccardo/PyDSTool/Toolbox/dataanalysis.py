"""
Data analysis utilities

Robert Clewley, August 2006
"""

from scipy.linalg import norm
from scipy import extract, mean, std
from PyDSTool import *
from numarray import random_array, array, dot, zeros, transpose, shape
from mdp.nodes import PCANode
from matplotlib.pylab import figure, plot, hold, polyfit, \
     bar, subplot, axes, title, xlabel, ylabel, draw, semilogy
from copy import copy


# -------------------------------------------------------------------------

# Principal Component Analysis utilities

def doPCA(in_pts, D, tol=.8):
    pts=copy(in_pts)
    pts.shape = (len(in_pts),D)
    p=PCANode(output_dim=tol)
    p.train(pts)
    p.stop_training()
    return p


def get_evec(p, dim):
    # extra call to array() to make it contiguous
    v = array(transpose(p.v)[dim-1])
    v.shape = (totdim, 3)
    return v


def pca_dim(pca, covering, data, refpt, tol=0.8):
    if type(refpt) != int:
        refpt = refpt[0]
    pca[refpt] = doPCA(nhd(data, covering, refpt), tol)
    print "PCA eigenvalues for nhd #%i:"%refpt
    print pca[refpt].d
    return len(pca[refpt].d)


def get_residual(pcanode, pts):
    pvals = pcanode.execute(pts, pcanode.get_output_dim())
    pproj = mat(pvals)*mat(transpose(pcanode.v))
    res = pts-array(pproj) # orthogonal complement of projection onto PCs
    resval_part = sqrt(sum(res**2))
    resval = sqrt(sum(resval_part**2))  # ?!
    return resval


def plot_PCA_residuals(data, D=None, newfig=True, marker='o'):
    if D is None:
        D = shape(data)[1]
    p=doPCA(data, D, D)
    spec=zeros((D,1),Float)
    for i in range(1,D):
        spec[i] = norm(p.d[:i])
    res = transpose(sqrt(spec[-1]**2-spec**2)/spec[-1])[0]
    print "2-norm of PCA spectrum =", spec[-1]
    if newfig:
        figure()
        style='k'+marker
    else:
        style=marker
    semilogy(res,style)
##    title('PCA residuals')
    xlabel(r'$\rm{Dimension}$',fontsize=20)
    ylabel(r'$\rm{PCA \ residual}$',fontsize=20)
    return p, res


def plot_PCA_spectrum(data, D):
    p=doPCA(data, D, D)
    spec=zeros((D,1),Float)
    for i in range(1,D):
        spec[i] = norm(p.d[:i])
    figure();plot(spec,'ko')
    title('2-norm of spectrum of singular values')
    xlabel('D')
    ylabel('|s|')
    return p


# -------------------------------------------------------------------------


def get_linear_regression_residual(pfit, x, y, weight='', w=0):
    res = 0
    if weight == 'lo':
        for i in range(len(x)):
            res += ((y[i]-(pfit[0]*x[i]+pfit[1]))/(i+1)/w)**2
    elif weight == 'hi':
        l = len(x)
        for i in range(l):
            res += ((y[i]-(pfit[0]*x[i]+pfit[1]))/(l-i+1)/w)**2
    else:
        for i in range(len(x)):
            res += (y[i]-(pfit[0]*x[i]+pfit[1]))**2
    return sqrt(res)


def fitline(x, y, lo=0, hi=1, doplot=1, quiet=1, linewidth=2):
    """Fitline takes the position of low and high fit limits
    and returns the slope. Integer indices in the x data set can
    be specified for the low and high limits, otherwise use a
    fraction between 0 and 1.
    
    In application to fractal dimension estimation, if
    x = log d (radius, a.k.a. inter-point distance) and
    y = log v (index, a.k.a. estimate of volume at a given radius)
    then the slope is the dimension estimate of the data set.
    """
    lendat = len(x)
    assert hi > lo
    if lo >= 0 and hi <= 1:
        loix = int(lendat*lo)
        hiix = int(lendat*hi)
    elif lo >=0 and hi > 0:
        loix = int(lo)
        hiix = int(hi)
    else:
        raise ValueError, "Invalid low and high fit limits"
    if quiet==0:
        print "lo ix %d, hi ix %d, max ix %d"%(loix, hiix, lendat-1)
    pfit = polyfit(x[loix:hiix],y[loix:hiix],1)
#    print pfit
    if doplot==1:
        plot([x[loix],x[hiix]],
             [x[loix]*pfit[0]+pfit[1],x[hiix]*pfit[0]+pfit[1]],
             linewidth=linewidth)
    return pfit[0]




def mean_bins(d):
    a = d.indepvararray*d.coordarray
    s = sum(a.tolist()[0])
    return s / sum(d.coordarray.tolist()[0])
    
def std_bins(d):
    di=d.indepvararray.tolist()
    dx=d.coordarray[0]
    ds = []
    for i in range(len(di)):
	ds.extend([di[i] for j in range(int(dx[i]))])
    return std(ds)

def closest_val(x, alist, eps=1e-7):
    for a in alist:
        if abs(x-a)<eps:
            return a
    print alist
    raise ValueError, "No nearby value found for %f in list at this tolerance!"%x

def bin_dict(d, interval=1, maxD=None):
    if maxD is None:
        maxD = int(around(max(d.tolist())[0]))+1
##    if interval < 1:
##        dbins_keys = interval*arange(0, round((maxD+1+interval)/interval), typecode='f')
##        dbins_keys = dbins_keys.tolist()
##        dbins = {}.fromkeys(dbins_keys, 0)
##        for x in d:
##            dbin = closest_val(x[0], dbins_keys)
##            dbins[dbin] = dbins[dbin]+1
##    else:
    dbins_keys = arange(0, maxD+1+interval, interval, typecode='f').tolist()
    dbins = {}.fromkeys(dbins_keys, 0)
    for x in d:
        dbin = closest_val(x[0], dbins_keys, interval)
        #dbin = int(around(x[0]/interval))*interval
        dbins[dbin] = dbins[dbin]+1
    return dbins

def bin_array(d, interval=1, maxD=None):
    dbins_dict = bin_dict(d, interval, maxD)
    x, y = sortedDictLists(dbins_dict)
    return (array(x), array(y))


def whiten(data):
    wdata=zeros(shape(data),Float)
    for d in range(shape(data)[1]):
        x = data[:,d]
        wdata[:,d] = (x-mean(x))/std(x)
    return wdata


def second_diff(data, i):
    return data[i+1]+data[i-1]-2*data[i]


def find_knees(data, tol=1., inlog=False, verbose=False):
    """
    inlog option finds knees in the logs of the data entries.
    tol=1. works well if inlog=False
    tol=0.3 works well if inlog=True
    """
    knee_ixs = []
    for i in range(1,len(data)-1):
        if inlog:
            frac = log(data[i+1])+log(data[i-1])-2*log(data[i])
        else:
            d2 = data[i+1]+data[i-1]-2*data[i]
            try:
                frac = d2/data[i]
            except ZeroDivisionError:
                frac = Inf
        if verbose:
            print i, data[i], frac
        if frac > tol and frac < Inf:
            knee_ixs.append((i,frac))
    if verbose:
        print "High second derivatives at: ", knee_ixs, "\n"
    knees = []
    found = False  # in a contiguous segment of high second derivatives
    curr_kixs = []
    old_kix = None
    for i in range(len(knee_ixs)):
        process = False
        kix, frac = knee_ixs[i]
        if verbose:
            print "Processing knee at index", kix
        if kix-1 == old_kix:
            if i == len(knee_ixs)-1:
                # this is the last index so have to process
                process = True
            if found:
                curr_kixs.append(i)
            else:
                curr_kixs = [i-1, i]
                found = True
        else:
            process = old_kix is not None
        if verbose:
            print old_kix, found, curr_kixs
        if process:
            if found:
                found = False
                if verbose:
                    print "Current knee indices:", curr_kixs,
                    print [knee_ixs[k] for k in curr_kixs]
                all_d2 = array([knee_ixs[k][1] for k in curr_kixs])
                ixs_sort = argsort(all_d2)
                max_ix = ixs_sort[-1]
                knees.append(knee_ixs[curr_kixs[max_ix]][0])
                curr_kixs = []
            else:
                if verbose:
                    print "Appending knee index", old_kix
                knees.append(old_kix)
        old_kix = kix
    # add final singleton index of high derivative, if any
    if not found and old_kix not in knees:
        if verbose:
            print "Appending knee index", old_kix
        knees.append(old_kix)
    return knees


def colormap(mag, cmin, cmax):
    """
    Return a tuple of floats between 0 and 1 for the red, green and
    blue amplitudes.
    Function originally named floatRgb and written by Alexander Pletzer,
    from the Python Cookbook.
    """

    try:
          # normalize to [0,1]
          x = float(mag-cmin)/float(cmax-cmin)
    except:
          # cmax = cmin
          x = 0.5
    blue = min((max((4*(0.75-x), 0.)), 1.))
    red  = min((max((4*(x-0.25), 0.)), 1.))
    green= min((max((4*math.fabs(x-0.5)-1., 0.)), 1.))
    return (red, green, blue)
