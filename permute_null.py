"""
Permute correlation null model.

.. author:: Albert Zhou

"""


import os, string
import numpy as np
from kagami.comm import *
from kagami.dtypes import Table


def _conv_fname(fname):
    vchrs = "-_.() %s%s" % (string.ascii_letters, string.digits)
    fname = ''.join(c for c in fname if c in vchrs)
    fname = fname.replace(' ','_')
    return fname

def _corr(x,y):
    vx, vy = (x != 0), (y != 0)
    if (np.sum(vx) < 3 or np.sum(vy) < 3) or \
       (np.sum(np.logical_and(~vx, ~vy))/np.sum(np.logical_or(vx, vy)) > 0.5) or \
       (np.isclose(np.var(x),0) or np.isclose(np.var(y),0)): return np.nan
    return np.corrcoef(x,y)[0,1]

def _perm_null(chm, tax, nperms):
    m = chm.shape[0]
    assert tax.shape[0] == m
    def _pm(n):
        idx = np.array([(
            np.random.choice(m, size = n, replace = False),
            np.random.choice(m, size = n, replace = False)
        ) for _ in range(nperms)])
        return np.array(smap(idx, unpack(lambda xi,xj: _corr(chm[xi],tax[xj]))))
    pmtx = np.array(smap(np.arange(3,m+1), _pm))
    return pmtx

def _perm_null_mp(params):
    cx, tx, nperms, ofname = params
    pm = _perm_null(cx, tx, nperms)
    np.savez(ofname, permatx = pm)

def _shoot(chmtab, taxtab, opath, chemid = None, taxnid = None, nperms = 10000, nprocs = None):
    chmtab = Table.loadhdf(chmtab)
    taxtab = Table.loadhdf(taxtab)
    
    rids = np.array(
        # sorted(np.intersect1d(chmtab.rows_, taxtab.rows_), key = lambda x: int(x.split('_',1)[-1]))
        sorted(np.intersect1d(chmtab.rows_, taxtab.rows_), key = lambda x: int(x))
    )
    chmtab, taxtab = chmtab[rids], taxtab[rids]
    
    chemid = [chemid] if available(chemid) else chmtab.cols_
    taxnid = [taxnid] if available(taxnid) else taxtab.cols_
    
    checkOutputDir(opath)
    
    pms = [(
        chmtab[:,ci].X_[:,0], taxtab[:,ti].X_[:,0], nperms,
        os.path.join(opath, _conv_fname(f'permutation_null_{ci}_{ti}_nperms_{nperms}.npz')),
    ) for ci in chemid for ti in taxnid]
    pmap(pms, _perm_null_mp, nprocs = nprocs)
    
    
# main
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('chmtab',  help = 'chemical table')
    parser.add_argument('taxtab',  help = '(summarised) taxon table')
    parser.add_argument('opath',  help = 'results output path')
    parser.add_argument('-c', '--chem', help = 'chemical name', default = None)
    parser.add_argument('-t', '--tax', help = 'taxon name', default = None)
    parser.add_argument('-p', '--npms', help = 'number of permutations', type = int, default = 10000)
    parser.add_argument('-n', '--nprocs', help = 'number of processors', type = int, default = None)    
    args = parser.parse_args()

    _shoot(args.chmtab, args.taxtab, args.opath, chemid = args.chem, taxnid = args.tax, nperms = args.npms, nprocs = args.nprocs)
    