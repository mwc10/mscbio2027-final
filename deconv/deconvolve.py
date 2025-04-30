import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import anndata as ad
import polars as pl
from anndata import AnnData

from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.sparse import csc_array

import multiprocessing
from pathlib import Path
import argparse

# https://github.com/gokceneraslan/fit_nbinom/blob/master/fit_nbinom.py
# ??

MAX_OPT_ITER = 25

def negbinom_loglik(x, mu, size) :
    '''
    Negative binomial log-likelihood
    '''
    # display(f'x: {x.shape} \n mu: {mu.shape} \n size: {size.shape}')

    return -np.sum(
        gammaln(x + size)
        - gammaln(size)
        - gammaln(x + 1)
        + size * np.log(size / (size + mu))
        + x * np.log(mu / (size + mu))
    )

def fit_gene_weights(x_col, f, sc, sigma, init):
    '''
    Optimize weights for a single gene (col of X) 
    '''

    def objective(p_log, x_col, f, sc):
        p = np.exp(p_log)
        sigma = p[0] + 0.2
        w = p[1:]
        mu = sc[:, None] * (f @ w)

        # print('mu pre-flat', mu.shape)
        return negbinom_loglik(x_col, mu, np.minimum(sigma, 1e3))
    
    init_log = np.log(np.concatenate([[max(sigma - 0.2, 1)], init]))
    # match R's optim default
    opts = {'maxiter': MAX_OPT_ITER}
    res = minimize(objective, init_log, method='BFGS', args=(x_col, f, sc), options=opts)
    # res = minimize_scalar(objective, init_log, method='BFGS', args=(x_col, f, sc), options=opts)
    p_opt = np.exp(res.x)

    return res.fun, p_opt[0] + 0.2, p_opt[1:]

def wrap_fgw(args: tuple):
    i, x_col, f, sc, sigma, init = args

    # print(f'starting {i}')

    return (i,) + fit_gene_weights(x_col, f, sc, sigma, init)

def expression_in_annotations(x, f, sigma=None, max_iter=5, verbose=True, cores=None):
    n_spots, n_genes = x.shape
    n_factors = f.shape[1]

    if sigma is None:
        sigma = np.ones(n_genes)
    sigma = np.clip(sigma, 0.5, 100)

    # Initialize W and scaling factors
    # W = np.ones((n_factors, n_genes))
    # W = np.ones((n_factors, n_genes))
    # they do!
    #   W = matrix(colSums(x), ncol=ncol(x), nrow=ncol(f), byrow=TRUE)/sum(f); # Init W (components, col=genes)
    W = np.tile(x.sum(axis=0), (n_factors, 1)) / f.sum()
    print(f'W init shape -> {W.shape} as f{(n_factors, n_genes)}')

    #W_init = W.copy()
    sc = np.clip(np.mean(x, axis=1) / np.mean(f @ W, axis=1), 0.1, 10)
    print(f'initial sc: {sc.shape}')

    for iter in range(max_iter):
        if verbose:
            print(f"Iteration {iter + 1}")

        # Step 1: Fit weights W for each gene (columns of x)
        W_new = np.zeros_like(W)
        new_sigma = np.zeros(n_genes)
        total_ll = 0

        # parallizeable
        # itr = ((i, x[:, i], f, sc, sigma[i], W[:, i]) for i in range(n_genes))
        # res = map(wrap_fgw, itr)
        # completed = 0
        # for i, ll, s, w in res:
        #     W_new[:, i] = w
        #     new_sigma[i] = s
        #     total_ll += ll

        #     completed += 1
        #     print(f'res {i}')

        #     if verbose and (completed-1 % 50) == 0:
        #         print(f'     finished {completed} genes of {n_genes} -> {(completed / n_genes):0.2f}')

        with multiprocessing.Pool() as pool:
            itr = ((i, x[:, i], f, sc, sigma[i], W[:, i]) for i in range(n_genes))
            res = pool.imap_unordered(wrap_fgw, itr, chunksize=2)
            completed = 0
            for i, ll, s, w in res:
                W_new[:, i] = w
                new_sigma[i] = s
                total_ll += ll

                completed += 1
                # print(f'res {i}')

                if verbose and (completed-1 % 250) == 0:
                    print(f'     finished {completed} genes of {n_genes} -> {(completed / n_genes):0.2f}')


        # for i in range(n_genes):
        #     ll, s, w = fit_gene_weights(x[:, i], f, sc, sigma[i], init=W[:, i])
        #     W_new[:, i] = w
        #     new_sigma[i] = s
        #     total_ll += ll

            

        sigma = np.clip(new_sigma, 0.5, 1e3 - 10)
        W = W_new
        if verbose:
            print(f"  W log-likelihood: {-total_ll:.2f}")

        # Step 2: Fit scaling factors for each spot
        def fit_scaling_factor(i):
            def objective(log_sc):
                sc_i = np.exp(log_sc)
                mu = sc_i * (f[i, :] @ W)
                return negbinom_loglik(x[i, :], mu, sigma)
            opts = {'maxiter': MAX_OPT_ITER}
            res = minimize(objective, np.log(sc[i]), method='BFGS', options=opts)
            return np.exp(res.x[0])

        # probably also parallizeable
        sc = np.array([fit_scaling_factor(i) for i in range(n_spots)])

        w

    return {
        'W': W.T,  # Return shape: [genes × components]
        'sigma': sigma,
        'sc': sc,
    }


## Data stuff
def filter_components(df: pl.DataFrame):
    '''
    Match the normalization and selection of the categories done in the code
    '''
    keep = ['Fat tissue','in situ','Lactiferous duct','Lymphoid nodule','Necrosis',
    "Tumor", "Stroma", "Lymphocyte", 'Vessels']

    frac = [
        (pl.col(c) / pl.sum_horizontal(keep)).alias(c)
        for c in keep
    ]

    return df.select('tnbc_id', 'slide_id', *keep).with_columns(frac)

def sort_categories(df: pl.DataFrame, adata: AnnData) -> pl.DataFrame:
    '''
    Once `df` is filtered to one sample, sort the rows to match
    the counts matrix
    '''
    order = pl.DataFrame({
        'slide_id': adata.obs_names,
        'sort_order': list(range(len(cnt.obs_names)))
    })

    return (
        df.join(order, on='slide_id')
            .sort('sort_order')
            .select(df.columns)
    )


def cli():
    p = argparse.ArgumentParser(description='Deconvolve ST counts into annotation category factors')
    p.add_argument('')

if __name__ == '__main__':
    df_classif = pl.read_csv('../data/predicted_classifications.csv')
    cnt = ad.read_h5ad('../data/anndata/counts/cnt-1.h5ad')

    # multiprocessing.set_start_method('fork')

    test = csc_array(cnt.X)
    test_cat = (
        filter_components(df_classif)
            .filter(pl.col('tnbc_id') == 1)
            .pipe(sort_categories, cnt)
            .drop('tnbc_id', 'slide_id')
            .to_numpy()
    )

    X_in = test[:, test.sum(axis=0) > 5]

    # # test with only 250 spots
    # X_in = X_in[:250, :]
    # test_cat = test_cat[:250, :]

    res = expression_in_annotations(X_in.toarray(), test_cat)
