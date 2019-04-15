#!/usr/bin/env python3

"""

Generates paired synthetic Spatial Transcriptomics and Single Cell Data 

The ST-data uses an image as the "histology" defining regions with
different cell compositions of K-specified celltypes. 

The SC-data will be in the form of a balanced set of single cells
coming from either of the K-celltypes.

The synthetic data includes
- rate              - for each celltype and gene
- logodds           - for each gene
- true proportions  - for each spot and celltype
- st_count matrix   - for each spot and gene
- sc_count matrix   - for each single cell and gene
- sc_meta_data      - for each single cell
- visualization     - of synthetic data


"""

from PIL import Image
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import torch as t

from torch.distributions import NegativeBinomial as nb
from scipy.stats import dirichlet
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

import pandas as pd
import os.path as osp
import os
import re



plt.style.use('dark_background')

# Functions ----------------------------

def generate_st_array(rawimg,
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
    
    return arrcrd, pxcrd, sptidx, Image.fromarray(rawimg)

def make_fake_data(sptidx,
                   arrcrd,
                   alpha,
                   n_sc ,
                   n_celltypes,
                   n_genes,
                   upper_cell_limit):
    
    out = {}
    
    # generate parameter values for genes and celltyps
    n_spots = sptidx.shape[0]
    probs = np.random.uniform(0,1,(n_genes,1)) # probability in nb-distribution
    logits = np.log(probs / (1.0-probs)) # logodds in nb-distribution
    # TODO: Think about this
    rates = np.random.gamma(3,0.2, size = (n_genes,n_celltypes)) # rates in nb-distribution
    
    # region labels and cardinality 
    idx, n_members = np.unique(sptidx, return_counts = True) 
    
    # propotions and count matrix 
    props = np.zeros((n_spots,n_celltypes))
    cmat =  np.zeros((n_spots,n_genes))

    sigma = 0.1 # scaling factor
    probs = t.tensor((probs.astype(np.float32)))

    # generate st-data
    for k,(lab,mem) in enumerate(zip(idx,n_members)):
        # get spots within morphological region
        inset = np.where(sptidx == lab)[0]
        # probability of observig celltypes within spot of region 
        pvals = dirichlet.rvs(np.repeat(alpha,n_celltypes))[0]
        # number of cells at each spot
        n_cells = np.random.randint(1,upper_cell_limit, size = mem)
        
        # number of cells from each celltype in all spots of region 
        celltypes  = np.vstack([np.random.multinomial(n = n_cells[x], pvals = pvals).reshape(1,-1) \
                                           for x in range(mem)])
        
        # true proportions of celltypes within spots of region
        props[inset,:] = celltypes / celltypes.sum(axis=1).reshape(-1,1)    
        
        # iterate over each celltype
        for z in range(n_celltypes):
            # use rate for specific cell type
            s_total_count = t.tensor((sigma * rates[:,z]).astype(np.float32).reshape(-1,1))
            # iterate over each spot in the region 
            for spt, j in enumerate(celltypes[:,z]):
                    # generate counts from j cells from cellype z to spot spt 
                    samples = nb(total_count = s_total_count, probs = probs).sample(t.tensor([j]))
                    # add counts to total count matrix
                    cmat[inset[spt],:] += samples.sum(dim = 0).numpy().reshape(-1,)
     
    # make DataFrame for ST-count matrix
    index = pd.Index(list(map( lambda x: ''.join(['X',str(x[0]),'x',str(x[1])]), arrcrd)))
    colnames_cmat = [''.join(['Gene_',str(x)]) for x in range(n_genes)]
    df_cmat = pd.DataFrame(cmat, index = index , columns = colnames_cmat)
    
    # make DataFrame for ST-celltype proportions
    colnames_proportions = [''.join(['Celltype_',str(x)]) for x in range(n_celltypes)]
    df_proportions = pd.DataFrame(props, index = index, columns = colnames_proportions)
    
    # add DataFrames to output
    out.update({'st_cnt':df_cmat,
                'st_prop':df_proportions})
    
    # if number of single cell observations should be generated
    if n_sc > 0:
        # get number of cells from each celltype. Use uniform probability
        sc_mem = np.random.multinomial(n_sc, np.ones(n_celltypes) / n_celltypes )
        sc_mat = []
        
        for z in range(n_celltypes):
            # use 10x higher rates for sigle cells
            z_total_counts = 2*t.tensor(rates[:,z].astype(np.float32).reshape(-1,1)) 
            # sample cells from celltype z
            sc_mat.append(nb(total_count= z_total_counts, 
                                             probs = probs).sample(t.tensor([sc_mem[z]])).numpy())
            
        
        # format single cell matrix
        sc_mat = np.vstack(sc_mat)
        sc_mat = sc_mat.reshape(n_sc,n_genes)

        # shuffle single cells
        ridx = np.arange(n_sc)
        np.random.shuffle(ridx)
        
        # make DataFrame for SC-count matrix
        index_sc = pd.Index([''.join(['C',str(x)]) for x in range(n_sc)])
        df_scmat = pd.DataFrame(data = sc_mat[ridx,:],
                                index = index_sc,
                                columns = colnames_cmat)

        # make DataFrame for SC-meta data
        idx_ct = np.vstack([np.repeat([k],n).reshape(-1,1) for (k,n) in enumerate(sc_mem)]).reshape(-1,1)[ridx]
        bio_ct = np.array(['_'.join(['Celltype',str(k)]) for k in idx_ct]).reshape(-1,1)[ridx]
        df_meta = pd.DataFrame(np.hstack((idx_ct,bio_ct)),
                               index = index_sc,
                               columns = ['celltype','bio_celltype'] )
        
        # make DataFrame distribution Parameters
        df_rates = pd.DataFrame(data = rates*2.0, 
                                index = colnames_cmat,
                                columns = colnames_proportions)
        
        df_logits = pd.DataFrame(data = logits,
                                 index = colnames_cmat,
                                 columns = ['logits'] )
        
        # add DataFrames to output
        out.update({'sc_meta':df_meta,
                    'sc_cnt':df_scmat,
                    'R':df_rates,
                    'logits':df_logits})
    
    return out

def hex_2_rgb(hexval):
    
    value = hexval.lstrip('#')
    nchr = len(value)
    rgb = tuple(int(value[i:i + nchr // 3], 16) for i in range(0, nchr, nchr // 3))
    
    return rgb

if __name__ == '__main__':
    
    import datetime
    import argparse as arp
    
    parser = arp.ArgumentParser()
    
    parser.add_argument('-i','--image',
                        required = True,
                        help = ('Link to image',
                                ' file to be used.',
                                ' if larger than 625x625',
                                ' pixels center-cropping',
                                ' is performed.'
                                ))
    
    parser.add_argument('-o','--odir',
                        required = False,
                        default = '',
                        help = ('output directory',
                                ' if none is specified',
                                ' a directory will be',
                                ' created in the same folder'
                                ' as the input directory',
                                ))
    
    parser.add_argument('-g','--n_genes',
                        required = False,
                        type = int,
                        default = 500,
                        help = ('number of genes to',
                                ' use in the synthetic',
                                ' data'
                                ))
    
    
    parser.add_argument('-s','--n_single_cells',
                        required = False,
                        type = int,
                        default = 2500,
                        help = ('number of single',
                                ' cells to use in'
                                ' synthetic data',
                                ' set to negative',
                                ' in order to only generate'
                                ' st-data.'
                                ))
    
    parser.add_argument('-z','--n_celltypes',
                        required = False,
                        type = int,
                        default = 5,
                        help = ('number of celltypes',
                                ' to use within tissue',
                                ' default is 5.'
                                ))
    
    
    parser.add_argument('-r','--regions',
                        required = False,
                        type = int,
                        default = 3,
                        help = ('number of regions within',
                                ' provided image. Do not count'
                                ' background.'
                                ))
    
    parser.add_argument('-u','--upper_cell_limit',
                        required = False,
                        type = int,
                        default = 30,
                        help = ('maximum number of',
                                ' cells allowed to be',
                                ' observed at each spot.'
                                ))

    parser.add_argument('-t','--tag',
                    required = False,
                    default = '',
                    help = ('tag to be attached',
                            ' to results. If none',
                            ' is specified then a',
                            ' date and time based',
                            ' tag is generated.'
                            ))
    
    
    parser.add_argument('-a','--alpha',
                    required = False,
                    type = float,
                    default =  0.5,
                    help = ('alpha value to use',
                            ' in dirichlet prior',
                            ' for celltypes within'
                            ' each morphological region',
                            ' Default is 0.5',
                            ))
    
    parser.add_argument('-b','--background',
                        required =False,
                        default = '#808080',
                        help = ('HEX-value of',
                                ' background color',
                                ' default is #808080',
                                ))

    args = parser.parse_args()

# Generate Data ---------------------------    

    
    imgpth = args.image
    
    # add tag to results for uniqueness
    if args.tag:
        tag = args.tag
    else:
        tag = re.sub('|\.|:','',str(datetime.datetime.today())).replace(' ','-')
    
    # set output directory
    if not args.odir:
        odir = osp.join(osp.dirname(args.image),'.'.join(['fake_tissue',tag]))
    else:
        odir = args.outdir
    
    # create output diretory if it does not exists
    if not osp.isdir(odir):
        os.mkdir(odir)
        
    # read from input 
    n_regions = args.regions
    rawimg = Image.open(imgpth).convert('RGB')
    
    # generate spot coordinates based on image
    arrcrd, pxcrd, sptidx, newimg = generate_st_array(rawimg, n_regions)
    
    # save image used in morphology generation
    newimg.save(osp.join(odir,'.'.join(['image',tag,'png'])))
    
    # make DataFrame for array-to-pixel conversion including region class
    spot_data = dict(x = arrcrd[:,0],
                y = arrcrd[:,1],
                pixel_x = pxcrd[:,0],
                pixel_y = pxcrd[:,1],
                region = sptidx)
    
    df_spot = pd.DataFrame(spot_data, columns = spot_data.keys())
    
    # generate fake data
    data = make_fake_data(sptidx = sptidx,
                          arrcrd = arrcrd,
                          alpha = args.alpha,
                          n_celltypes=args.n_celltypes,
                          n_genes=args.n_genes,
                          n_sc = args.n_single_cells,
                          upper_cell_limit = args.upper_cell_limit)
    
    
# Visualize ----------------------------
    
    # dimensionality reduction procedures
    dr3 = PCA(n_components = 3)
    dr2 = TSNE(n_components=2)
    
    rel_freq = lambda x: x.values / x.values.sum(axis=1).reshape(-1,1)
    
    # Get colors for st_data
    tr = dr3.fit_transform(rel_freq(data['st_cnt']))
    mx = np.max(tr,axis = 1).reshape(-1,1)
    mn = np.min(tr,axis = 1).reshape(-1,1)
    rgb_st = (tr -mn) / (mx - mn) 
    
    # get 2d t-SNE representation of SC-data
    sc_tsne = dr2.fit_transform(data['sc_cnt'].values)
    # color t-SNE by expression similarity
    rgb_sc = dr3.fit_transform(rel_freq(data['sc_cnt']))
    mx = np.max(rgb_sc,axis = 1).reshape(-1,1)
    mn = np.min(rgb_sc,axis = 1).reshape(-1,1)
    rgb_sc = (rgb_sc -mn) / (mx - mn) 
    
    # generate plots
    fig, ax = plt.subplots(1,4,figsize = (24,4))
    ax[0].scatter(pxcrd[:,0], pxcrd[:,1], c = rgb_st, s = 120, edgecolor = 'k', alpha = 0.8)
    for ii in np.unique(sptidx): 
        ax[1].scatter(pxcrd[sptidx == ii,0], pxcrd[sptidx == ii,1], s = 120, edgecolor = 'k', alpha = 0.8)
    
    ax[0].set_title('Dimensionality reduction of relative counts')
    ax[0].invert_yaxis()
    ax[1].set_title('Regions')
    ax[1].invert_yaxis()
    ax[2].set_title('Original Image')
    ax[2].imshow(rawimg)
    ax[3].scatter(sc_tsne[:,0],sc_tsne[:,1], c = rgb_sc)
    ax[3].set_title('t-SNE of Single Cell Data')
    for aa in ax: aa.axes.set_aspect('equal')
    # save figure
    fig.savefig(osp.join(odir,'.'.join(['viz',tag,'png'])))
