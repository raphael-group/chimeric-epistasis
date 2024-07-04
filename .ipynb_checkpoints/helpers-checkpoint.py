import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from itertools import combinations,chain,combinations_with_replacement,product,chain
from collections import Counter

from scipy.special import factorial


import seaborn as sns

import glob
import os

import pickle
from pprint import pprint
import re

#################################
# HELPERS
#################################

# Create powerset of input
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

# Create list of all partitions of input
def partition(collection):
    if len(collection) == 1:
        yield [ collection ]
        return

    first = collection[0]
    for smaller in partition(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
        # put `first` in its own subset 
        yield [ [ first ] ] + smaller


#################################
# EPISTASIS HELPERS
#################################
# Compute the multiplicative/chimeric measures for any subset of loci

# INPUT
# f_vals: a dictionary where key=tuple of mutated indices and value=fitness
# Example:
# f_vals = {(): 1.0, (2,): 1.37294933339, (1,): 2.40124349574, (1, 2): 1.75153255202, (0,): 0.0619096557142, (0, 2): 0.120596667317, (0, 1): 0.274083334811, (0, 1, 2): 0.210456846373}
# order: largest epistasis order to compute
# log: if True, return log(multiplicative measure)

# OUTPUT:
# dictionary where key=tuple of indices and value=multiplicative epistasis measure for those indices
# e.g. (0,1) -> \epsilon^M_{0,1}, or the pairwise multiplicative measure for mutations in indices 0,1
def mult_epistasis_subset(f_vals, order, log=False):
    mult_eps={}
    for o in range(order+1):
        for subset in combinations(max(f_vals.keys(), key=len),o):
            if len(subset) <= 1:
                mult_eps[subset]=f_vals[subset]
            else:
                numerator=f_vals[subset]
                denominator=1 # product of all smaller epistasis coefficients
                for r in range(len(subset)):
                    for ssb in combinations(subset,r):
                        denominator = denominator * mult_eps[ssb]
                mult_eps[subset] = np.round( numerator / denominator, 10)
    if log:
        for t in mult_eps:
            mult_eps[t]=np.log(mult_eps[t])
    return mult_eps


# Same as above, but for chimeric epistasis measure
# Uses formula for cumulant in terms of partitions https://en.wikipedia.org/wiki/Cumulant#Joint_cumulants
def chimeric_epistasis_subset(f_vals, order):
    # L=len([t for t in f_vals.keys() if len(t)==1])
    chim_eps={}
    chim_eps[()]=0
    for o in range(1,order+1):
        for subset in combinations(max(f_vals.keys(), key=len),o):
            cum_val=0
            # compute cumulant formula for subset
            for pi in partition(list(subset)):
                to_add = factorial((len(pi)-1)) * (-1)**(len(pi)-1)
                for B in pi:
                    to_add *= f_vals[tuple(B)]
                cum_val+= to_add
            chim_eps[subset]=cum_val
    return chim_eps

################################################
# Multiplicative and chimeric epistasis measures
################################################
# INPUT
# f_vals: a dictionary where key=tuple of loci and mutation at each loci, value=fitness
# Example:
# f_vals = {(0,0,0): 1.0, (0,0,1): 1.37294933339, (0,1,0): 2.40124349574, (0,1,1): 1.75153255202, (1,0,0): 0.0619096557142, (1,0,1): 0.120596667317, (1,1,0): 0.274083334811, (1,1,1): 0.210456846373}
# NOTE: you can also have different mutations at a locus, e.g. [1,2,0] means mutation 1 at locus 0 and mutation 2 at locus 1

# order: number of mutations in epistasis coefficient
# log: if True, return log(multiplicative measure)

# OUTPUT:
# (1) dictionary where key=list of mutations and value=multiplicative epistasis measure for those mutations
# e.g. [1,1,0] -> \epsilon^M_{0,1}, or the pairwise multiplicative measure for mutation 0 at locus 0 and mutation 1 at locus 1

# (2) dictionary w/ same format as above but with chimeric measures
# (3) dictionary w/ same format as above but with fitness standard deviations
def mult_and_chim_epistasis(f_vals, order, log=True, fitness_threshold=0.01):
    L=len(list(f_vals.keys())[0])
    mult_eps={}
    chim_eps={}
    f_stdev={}

    for mut_list in f_vals: # eg [0,5,6,0]
        if np.count_nonzero(mut_list) == order:
            # create dictionary to pass into mult_epistasis_subset
            nonzero_inds=np.nonzero(mut_list)[0] # eg [1,2]
            inds_to_mut={i:mut_list[i] for i in nonzero_inds} # eg {1:5,2:6}
            fv_to_pass={} # to pass into mult_epistasis_subset
            for ind_subset in powerset(nonzero_inds):
                genotype=np.zeros(L).astype(int)
                for i in ind_subset:
                    genotype[i]=inds_to_mut[i]
                if tuple(genotype) in f_vals:
                    fv_to_pass[tuple(ind_subset)]=f_vals[tuple(genotype)]

            if np.min(list(fv_to_pass.values()))>fitness_threshold and len(fv_to_pass.values())==2**order:
                mult_eps[tuple(genotype)]=list(mult_epistasis_subset(fv_to_pass, order, log=log).values())[-1]
                chim_eps[tuple(genotype)]=list(chimeric_epistasis_subset(fv_to_pass, order).values())[-1]
                f_stdev[tuple(genotype)]=np.std(list(fv_to_pass.values()))
    return mult_eps,chim_eps,f_stdev

def compute_std_vs_corr_sdf(results_dict,order_list=[3,4,5]):
    stdev_dict={o:{} for o in order_list}
    corr_dict={o:{} for o in order_list}
    sdfrac_dict={o:{} for o in order_list}

    for order in order_list:
        rs=results_dict[order]
        protein_list=list(rs.keys())
        for protein in protein_list:
            mult_eps,chim_eps,f_stdev=rs[protein]
            mult_vals=list(mult_eps.values())
            chim_vals=list(chim_eps.values())
            N=len(chim_vals)

            if N>10: # don't include proteins with too few measured interactions    
                stdev_dict[order][protein]=np.mean(list(f_stdev.values()))
                corr_dict[order][protein]=np.corrcoef(mult_vals,chim_vals)[0,1]
                sdfrac_dict[order][protein]=np.sum( np.sign(mult_vals) != np.sign(chim_vals) ) / len(chim_vals)
    return stdev_dict,corr_dict,sdfrac_dict

################################################
# Plotting functions
################################################

def plot_mult_vs_chim(results_dict, key, order_list=[3,4,5]):
    fig,axs=plt.subplots(1,len(order_list),figsize=(6*len(order_list),6))
    
    for ind,epi_order in enumerate(order_list):
        mult_dict,chimeric_dict,_=results_dict[epi_order][key]
        mult_vals=list(mult_dict.values())
        chimeric_vals=list(chimeric_dict.values())
        
        axs[ind].scatter(chimeric_vals,mult_vals,s=30,edgecolor='black')
        axs[ind].axvline(0,c='black',ls='--')
        axs[ind].axhline(0,c='black',ls='--')
        
        axs[ind].tick_params(axis='x', labelsize=20)
        axs[ind].tick_params(axis='y', labelsize=20)
        
        axs[ind].set_xlabel(f'{epi_order}-way log-multiplicative measure',fontsize=20)
        axs[ind].set_ylabel(f'{epi_order}-way chimeric measure',fontsize=20)
    
        corr=np.round( np.corrcoef(mult_vals,chimeric_vals)[0,1], 4)
        axs[ind].set_title(f'{epi_order}-way epistasis \n Pearson correlation: {corr}',fontsize=15)
        axs[ind].spines[['right', 'top']].set_visible(False)
    
    plt.tight_layout()

def plot_stdev_vs_corr_sdf(stdev_dict,corr_dict,sdf_dict):
    order_list=list(stdev_dict.keys())
    O=len(order_list)
    fig,axs=plt.subplots(2,O,figsize=(5*O,10))

    for order_ind,order in enumerate(order_list):
        xpts=list(stdev_dict[order].values())
        ypts_list=[list(corr_dict[order].values()), list(sdf_dict[order].values())]
    
        xlabel=f'Standard deviation of \n fitness {order}-tuples'
        ylabel_list=['Correlation', 'Sign disagreement fraction']
    
        for i, ypts in enumerate(ypts_list):
            ax=axs[i,order_ind]
            ax.scatter(xpts,ypts,edgecolor='black')
    
            ax.set_xlabel(xlabel, fontsize=18)
            ax.set_ylabel(ylabel_list[i], fontsize=18)
            
            ax.set_title(f'{order_ind+3}-way interactions',fontsize=20)
            sns.regplot(x=xpts,y=ypts,ax=ax,scatter_kws={'edgecolors': 'black'}, line_kws={'color': 'red','ls': '--'},
                        fit_reg=True,truncate=True,ci=67)
    
            ax.tick_params(axis='x', labelsize=16)  # For x-axis tick labels
            ax.tick_params(axis='y', labelsize=16)  # For y-axis tick labels
            ax.spines[['right', 'top']].set_visible(False)
    
    plt.tight_layout()