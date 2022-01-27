from utils import mtx_from_param_to_np
from utils import arr_from_param_to_np
from utils import np_mtx_to_param
from utils import np_arr_to_param
from scipy.integrate import odeint
from lmfit import Parameters, Parameter
import numpy as np

b = 1
tau = 10.0
g_data_tf = []


def f(x, t, ps):
    global g_data_tf
    alphas_p = []
    for i, x_i in enumerate(x):
        alphas_p_i = arr_from_param_to_np(f'alpha__{i}', ps)
        alphas_p.append(alphas_p_i)
    
    alphas_0_p = arr_from_param_to_np("alphaO", ps)
    miu_p = arr_from_param_to_np("miu", ps)
    Ms_p = arr_from_param_to_np("Ms", ps)
    dx = np.ones(x.shape)
    for i, x_i in enumerate(x):
        tf_idx = g_data_tf[i]
        if tf_idx.size != 0:
            tfs_val = x[g_data_tf[i]]
        else:
            tfs_val = np.ones(alphas_p[i].shape)
        y_i = alphas_0_p[i] + np.dot(alphas_p[i], tfs_val) 
        dx[i] = -miu_p[i]* x_i + Ms_p[i]/(1+ np.exp(-y_i))
    dx = np.nan_to_num(dx)
    return dx

def g(t, x0, ps):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x0 = arr_from_param_to_np('x0', ps)
    x = odeint(f, x0, t, args=(ps,))
    return x

def residual(ps, ts, data):
    x0 = arr_from_param_to_np('x0', ps)
    model = g(ts, x0, ps)
    sq = model - np.transpose(data)
    sq = np.multiply(sq, sq)
    return sq.ravel()
    
def init_params(gene_names, data_Gene, data_TF):
    params = Parameters()
    num_genes = len(gene_names)
    
    np_arr_to_param('alphaO', np.random.random(num_genes), params, -10000.0, 10000.0)
    np_arr_to_param('miu', np.random.random(num_genes), params, 0.0, 1.0)
    np_arr_to_param('x0', data_Gene[:,0], params, 0.0, 1.0, vary=False)
    np_arr_to_param('Ms', np.random.random(num_genes), params,  0.0, 10000)
    global g_data_tf
    g_data_tf = data_TF
    for idx, gene in enumerate(gene_names):
        tflen = len(data_TF[idx])
        np_arr_to_param(f'alpha__{idx}', 
                        np.random.random(tflen), 
                        params, -100000.0, 100000.0)
   
    return params