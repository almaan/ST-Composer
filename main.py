#!/usr/bin/env python3


import os.path as osp
import os
import datetime
import re
import argparse as arp
import logging

import pandas as pd
import numpy as np

import PIL.Image as Image
from sklearn.cluster import KMeans

    


def st_array_from_image(rawimg,
                        ncolors,
                        bg = 128):

    # Synthetic data specs
    nspots = 23
    side = 625
    ncolors = ncolors + 1
    ch = 3

    # check image size and crop to 625x625 if necessary
    w,h = rawimg.size
    if w != side or h != side:
        oriW = w/2.0
        oriH = h/2.0
        dr =   side/2.0
        box = (oriW-dr,oriH-dr,oriW+dr,oriH + dr)
        rawimg = rawimg.crop(box)

    w, h = rawimg.size
    rawimg = np.array(rawimg)

    # extract different "morphological" regions based on color
    km = KMeans(ncolors)
    idx = km.fit_predict(rawimg.reshape(w*h,ch))

    # get current label of background and set to -1
    bglab = np.round(np.mean(idx[np.all(rawimg.reshape(w*h,ch) > (bg-1), axis = 1) * \
               np.all(rawimg.reshape(w*h,ch) < (bg+1), axis = 1)])).astype(int)

    idx[ idx == bglab] = -1

    # map region id to coordinate id
    newimg = idx.reshape(w,h)

    # generate array coordinates
    arrX = np.arange(nspots + 3)
    arrY = np.arange(nspots + 3)
    xx,yy = np.meshgrid(arrX[1:-1],arrY[1:-1])
    arrcrd = np.hstack((xx.reshape(xx.shape[0]*xx.shape[0],1),
                        yy.reshape(yy.shape[0]*yy.shape[0],1))).astype(int)

    # generate pixel coordinates
    pixX = np.ceil(arrX * w/(nspots + 2) )
    pixY = np.ceil(arrY * h/(nspots + 2) )
    XX,YY = np.meshgrid(pixX[1:-1],pixY[1:-1])
    pxcrd = np.hstack((XX.reshape(XX.shape[0]*XX.shape[0],1),
                       YY.reshape(YY.shape[0]*YY.shape[0],1))).astype(int)

    # get region id for each spot
    sptidx = np.array([newimg[pxcrd[x,1],pxcrd[x,0]] for x in range(pxcrd.shape[0])])

    # remove spots not covered by the "tissue"
    arrcrd = arrcrd[sptidx > -1,:]
    pxcrd = pxcrd[sptidx > -1,:]
    sptidx = sptidx[sptidx > -1]

    return arrcrd,sptidx


def MakeSCSets(n_cells,
               mixture,):
    
    """Make Single Cell dataset"""
    
    
    
    n_types = mixture.Z
    n_genes = mixture.G
    pvals = np.ones(n_types) / n_types
    # set number of samples from each cell type
    nmem = np.random.multinomial(n_cells,pvals)
    
    yobs = np.zeros((n_cells,n_genes))
    meta = np.zeros((n_cells,1),dtype = int)
    start_idx = 0
    
    # get observations from all cell types
    for z in range(n_types):
        if not mixture.ct[z].DE:
            mixture.ct[z].setDEGenes()
            
        yobs[start_idx:start_idx+nmem[z],:] = mixture.ct[z].getCounts(nmem[z],bias = False)
        meta[start_idx:start_idx+nmem[z],:] = z
        
        start_idx += nmem[z]
    
    
    return {'counts':yobs, 'meta':meta}
    
        
            

def MakeSTSets(labels,
               mixture,
               conc,
               upper_cell_bound):
    
    """Make Spatial Transcriptomics dataset"""
    
    n_spots = labels.shape[0]
    n_genes = mixture.G
    n_types = mixture.Z
    
    sidx = np.argsort(labels)
    xobs = np.zeros((n_spots,n_genes))
    wobs = np.zeros((n_spots,n_types))
    
    regions,mems = np.unique(labels,return_counts = True)

    start_idx = 0
    for region in range(regions.shape[0]):
        cps = np.ceil(np.random.random()*upper_cell_bound).astype(int)
        logger.info(f'region {region} will have {cps} cells per spot')
        mixed_data = mixture.getCounts(cps,
                                       mems[region],
                                       conc,
                                       partial = True)
        
        
        xobs[start_idx:start_idx + mems[region],:] = mixed_data['counts']
        wobs[start_idx:start_idx + mems[region],:] = mixed_data['proportions']
        start_idx += mems[region]
    
    xobs = xobs[sidx,:]
    wobs = wobs[sidx,:]
    
    
    return {'counts':xobs, 'proportions':wobs}

class Mixture:
    """Mixture of different cell types
    
    A class which represents a mixture of cell types.
    Will generate observations of joint expression profiles.
    These observations can represented as partial (scaled) or
    raw (non-scaled).
    
    """
    
    def __init__(self,
                 n_types,
                 base_ct,
                 ):
        
        self.Z = n_types
        self.G = base_ct.G
        self.ct = [base_ct]
        self.ct[0].setDEGenes()
        
        # generate multiple cell types based on base
        for z in range(1,self.Z):
            tmp = CellType(self.G,base_ct.alpha,base_ct.beta)
            # share non-DE genes properties
            tmp.__dict__.update(self.ct[0].__dict__)
            # set DE genes for each cell type
            tmp.setDEGenes()
            # store new cell type
            self.ct.append(tmp)
            
    def __call__(self,z):
        """Get z cell types of mixture"""
        return self.ct[z]
    
    def getCounts(self,
                  n_cells,
                  n_obs,
                  conc,
                  partial = False,
                  ):
        """Get observations from mixture"""
        
        # get probability of picking each cell type
        probs = np.random.dirichlet(conc)
        # draw cells from each cell type
        nmem = np.random.multinomial(n_cells,
                                     probs,
                                     size = n_obs,
                                     )
        
        # get proportions of cell type
        prop = nmem / nmem.sum(axis=1).reshape(-1,1)
        # generate observations
        xobs = np.zeros((n_obs,self.G))
        # iterate over observations
        for obs in range(n_obs):
            # iterate over cell types
            for z in range(self.Z):
                xobs[obs,:] += self.ct[z].getCounts(nmem[obs,z]).sum(axis = 0)
        
        # scale data to mimic partial observations
        if partial:
            alpha = (np.random.random(n_obs)).reshape(-1,1)
            xobs = np.floor(xobs / n_cells * alpha)
        
        return {'counts' : xobs, 'proportions':prop}


class CellType:
    """Cell type class"""
    
    def __init__(self,
                 n_genes,
                 alpha = 1.0,
                 beta = 1.0,
                 p_dge = 0.05,
                 o_prob  = 0.01,
                 o_mu = 0,
                 o_sigma = 1,
                 ):
        
        self.G = n_genes
        self.alpha = alpha
        self.beta = beta
        self.p_dge = p_dge
        # original mean values
        lmbda = np.random.gamma(alpha,beta,self.G)
        # outlier indicator
        omega = np.random.binomial(1,o_prob, self.G)
        # outlier factor
        l_psi = self.scalingFactors(selected = omega,
                                    neg_prob = 0.5,
                                    mu = o_mu,
                                    sigma = o_sigma
                                    )
        # add outliers
        self.olmbda = np.median(lmbda)*omega*l_psi
        self.olmbda += lmbda*(1-omega).round(0)
        
        # set poisson parameter
        self.lmbda = self.olmbda
        self.DE = False
   
    def setDEGenes(self,
                   z_mu = 0,
                   z_sigma = 1,
                   neg_prob = 0.5,
                   ):
        
        self.DE = True
        # specifiy which genes should be DE
        self.zidx = np.round(np.random.binomial(1,self.p_dge, self.G))
        # get scaling factors for DE genes
        z_psi = self.scalingFactors(selected = self.zidx,
                                    neg_prob = neg_prob,
                                    mu = z_mu,
                                    sigma = z_sigma )
        
        # scale DE-genes
        zlmbda = np.mean(self.olmbda)*self.zidx*z_psi
        zlmbda += self.olmbda*(1-self.zidx).round(0)
        
        # update poission parameter to include DE-gene configuration
        self.lmbda = zlmbda
    
    def scalingFactors(self,
                       selected,
                       neg_prob,
                       mu  = 0,
                       sigma = 1):
        
        """Generate Scaling Factors
        
        Will return a vector with scaling factors at the
        position of the selected elements and zeros for
        those not selected.
        
        Parameters
        ----------
        selected : np.ndarray
            one-hot encoded vector indicating elements to be scaled
        neg_prob : float
            probability of downscaling
        mu : float
            mu-value in LogNormal distribution
        sigma : float
            sigma-value in LogNormal distribution
        
        """
        
        # get total number of features
        n_feats = selected.shape[0]
        # get number of features
        n_sel = selected.sum()
        # set direction of scaling
        dire = np.power(-1, np.random.binomial(1,neg_prob,n_sel))
        # draw scaling factors from log-normal distribution
        psi = np.random.lognormal(mu,sigma,n_sel)
        # maintain original direction probability structure
        dire[psi < 1] = -1*dire[psi < 1]
        sf = np.zeros(n_feats)
        # fill size factor vector
        sf[selected.astype(bool)] = np.power(psi, dire)
        
        return sf
    
    def dropGenes(self,
                  n_cells,
                  d_prob = 0.01):
        
        """Generate mask for genes to be dropped"""
        
        mask = np.random.binomial(1,
                                  d_prob,
                                  size = (self.G,n_cells))
        mask = (1 - mask).round(0)
        
        return mask
    
    def getCounts(self,
                 n_cells,
                 l_mu = 0.0,
                 l_sigma = 1.0,
                 bias = False,
                 d_prob= 0.01,
                 ):
        """Get Observations"""
        
        # get library scaling factors
        l_psi = self.scalingFactors(selected = np.ones(n_cells).astype(int),
                                    neg_prob = 0.0,
                                    mu = l_mu,
                                    sigma = l_sigma,
                                    )
            
        # create GxC matrix of poission parameter
        lmbda_mat =  np.repeat(self.lmbda.reshape(-1,1),n_cells, axis = 1)
        lmbda_mat = lmbda_mat * l_psi.reshape(1,-1)
        
        # introduce gene bias if specified
        if bias:
            beta = np.random.random(self.G).reshape(self.G,1)
            lmbda_mat *= beta
        
        # draw observations from poisson distribution
        xobs = np.random.poisson(lmbda_mat)
        # drop genes 
        yobs = xobs * self.dropGenes(n_cells,d_prob)
        # return CxG matrix
        yobs = yobs.T
        
        return yobs
    
def main(n_genes,
         rawimg,
         n_regions,
         n_types,
         n_sc_cells,
         p_dge,
         upper_cell_bound,
         gamma_alpha,
         gamma_beta,
         conc,
         ):

    # generate st-array from image
    array,labels = st_array_from_image(rawimg,
                                       n_regions)
    
    logger.info('converted image to st-array format')
    
    # generate coordinate based spot-indices 
    crd = ['x'.join([str(array[x,0]),str(array[x,1])]) for x\
          in range(array.shape[0]) ]
    crd = pd.Index(crd)
    
    # generate base cell type to use in mixture
    base_ct = CellType(n_genes, gamma_alpha, gamma_beta,p_dge)
    logger.info(f'generated cell base type')
    # initiate cell mixture
    mixture = Mixture(n_types, base_ct)
    logger.info(f'assembled mixture of sc-types')
    
    # container for results
    data = dict()
    # generate ST-data based on image segmentation
    st_data = MakeSTSets(labels,mixture,conc, upper_cell_bound)
    # format ST-data
    genes = pd.Index(['gene_' + str(g) for g in range(n_genes)])
    
    data.update({'st_cnt':pd.DataFrame(st_data['counts'],
                                       index = crd,
                                       columns = genes)})
    
    celltype = ['Celltype_' + str(x) for x in range(n_types)]
    celltype = pd.Index(celltype)
    
    data.update({'st_proportions' : pd.DataFrame(st_data['proportions'],
                                                 index = crd,
                                                 columns = celltype)})
    logger.info('generated st-data')
    
    # generate SC-data from same mixture as ST-data
    sc_data = MakeSCSets(n_sc_cells, mixture)
    # format SC-data
    cellnames = pd.Index(['C' + str(c) for c in range(n_sc_cells)])
    data.update({'sc_cnt' : pd.DataFrame(sc_data['counts'],
                                         index = cellnames,
                                         columns = genes)})
    
    bio_celltype = np.array(['Celltype_' + str(int(z[0])) for z \
                             in sc_data['meta']]).reshape(-1,1)
    
    meta = np.hstack((sc_data['meta'].reshape(-1,1),bio_celltype))
    metacols = pd.Index(['celltype','bio_celltype'])
    
    data.update({'sc_meta' : pd.DataFrame(meta,
                                          index = cellnames,
                                          columns = metacols)})
    
    logger.info('generated sc-data')
    
    return data

if __name__ == '__main__':
    
    prs = arp.ArgumentParser()
    
    prs.add_argument('-i','--image',
                     required = True,
                     type = str,
                     default  = None,
                     help = ''
                     )
    
    
    prs.add_argument('-o','--odir',
                     required = True,
                     type = str,
                     default  = None,
                     help = ''
                     )
    
    
    prs.add_argument('-scc','--n_sc_cells',
                     type = int,
                     required = True,
                     default  = None,
                     help = ''
                     )
    
    prs.add_argument('-nr','--n_regions',
                     type = int,
                     required = True,
                     default  = None,
                     help = ''
                     )
    
    
    prs.add_argument('-nt','--n_types',
                     type = int,
                     required = True,
                     default  = None,
                     help = ''
                     )
    
    
    prs.add_argument('-ng','--n_genes',
                     type = int,
                     required = True,
                     default  = None,
                     help = '',
                     )
    
    
    prs.add_argument('-pde','--p_dge',
                     type = float,
                     required = False,
                     default  = 0.05,
                     help = '',
                     )
    
    
    prs.add_argument('-dc',
                     '--dirichlet_concentration',
                     type = float,
                     default = 1,
                     help = '',
                     )
    
    prs.add_argument('-t','--tag',
                     type = str,
                     default  = None,
                     help = '',
                     )
    
    
    prs.add_argument('-ga','--gamma_alpha',
                     type = float,
                     default  = 1.0,
                     help = '',
                     )
    
    
    prs.add_argument('-gb','--gamma_beta',
                     type = float,
                     default  = 1.0,
                     help = '',
                     )
    
    
    prs.add_argument('-ucb','--upper_cell_bound',
                     type = int,
                     default  = 30,
                     help = '',
                     )
    
    args = prs.parse_args()
    
    # setup logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    c_handler = logging.StreamHandler()
    c_handler.setLevel(logging.DEBUG)
    c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    c_handler.setFormatter(c_format)
    logger.addHandler(c_handler)

    logger.info(f'Using arguments : {args}')

    # set identifier for synthetic set
    if not args.tag:
        tag = str(datetime.datetime.now())
        tag = re.sub('-|:| |\.|','',tag)
    else:
        tag = args.tag
    
    if not osp.exists(args.odir):
        os.mkdir(args.odir)
        logger.info(f'creating output directory {args.odir}')
    
    logger.info(f'using image {args.image}')
    # load specified image
    rawimg = Image.open(args.image).convert('RGB')
    conc = np.ones(args.n_types) * args.dirichlet_concentration
    
    main_args = dict(n_genes = args.n_genes,
                     rawimg = rawimg,
                     n_regions = args.n_regions,
                     n_types = args.n_types,
                     n_sc_cells = args.n_sc_cells,
                     p_dge = args.p_dge,
                     upper_cell_bound = args.upper_cell_bound,
                     gamma_alpha = args.gamma_alpha,
                     gamma_beta = args.gamma_beta,
                     conc = conc,
                     )
    
    data = main(**main_args)
    # save data
    for k,v in data.items():
        opth = osp.join(args.odir,'.'.join([tag,k,'tsv']))
        v.to_csv(opth, sep = '\t',header = True, index = True)
        logger.info(f'saved {k} file to >> {opth} ')
    
    logger.info('completed synthetic data generation')