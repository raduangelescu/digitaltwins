from utils import mtx_from_param_to_np
from utils import arr_from_param_to_np
from utils import np_mtx_to_param
from utils import np_arr_to_param
from scipy.integrate import odeint
from lmfit import Parameters, Parameter
import numpy as np

g_data_tf = []
x0 = []

def f(x, t, ps):
    global g_data_tf #, drug_gene_inf
    global alphas_p, alphas_0_p, miu_p, Ms_p, x0
    drug_gene_inf = arr_from_param_to_np("drug_gene_inf", ps)
    dx = np.ones(x.shape)
    b = ps['b'].value 
    tau = ps['tau'].value
    drug_time = b * np.exp(-t/tau)
    for i, x_i in enumerate(x):
        tf_idx = g_data_tf[i]
        if tf_idx.size != 0:
            tfs_val = x[g_data_tf[i]]
        else:
            tfs_val = np.ones(alphas_p[i].shape)
        y_i = alphas_0_p[i] + np.dot(alphas_p[i], tfs_val) 
        
        #drug_influence = drug_time * drug_gene_inf[i] * x[i] 
        dx[i] = -miu_p[i]* x_i + Ms_p[i]/(1+ np.exp(-y_i)) + drug_time * drug_gene_inf[i]
    dx = np.nan_to_num(dx)
    return dx

def g(t, x0, ps):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = odeint(f, x0, t, args=(ps,))
    return x

def residual(ps, ts, data):
    x0 = arr_from_param_to_np('x0', ps)
    model = g(ts, x0, ps)
    sq = model - np.transpose(data)
    sq = np.multiply(sq, sq)
    return sq.ravel()
    
def init_params_sim(gene_names, data_Gene, data_TF, p, dg_inf):
    global alphas_p, alphas_0_p, miu_p, Ms_p, x0, drug_gene_inf
    global g_data_tf
    g_data_tf = data_TF

    num_genes = len(gene_names)
    alphas_p = []
    for i, gname in enumerate(gene_names):
        alphas_p_i = arr_from_param_to_np(f'alpha__{i}', p)
        alphas_p.append(alphas_p_i)

    alphas_0_p = arr_from_param_to_np("alphaO", p)
    miu_p = arr_from_param_to_np("miu", p)
    Ms_p = arr_from_param_to_np("Ms", p)
    x0 = data_Gene[:,0]
    
    #drug_gene_inf = dg_inf
    params = Parameters()
    num_genes = len(gene_names)
    
    np_arr_to_param('x0', data_Gene[:,0], params, 0.0, 1.0, vary=False)
    
    np_arr_to_param('drug_gene_inf', np.random.random(num_genes), params, -10000.0, 10000.0)
    params.add('b', value = 1.0, min = 0.0, max = 1000.0)
    params.add('tau', value = 10.0, min=0.0, max = 1000.0)

    return params