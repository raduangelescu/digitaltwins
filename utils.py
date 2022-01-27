import logging
import numpy as np
import networkx as nx
from lmfit import Parameters, Parameter


def get_logger(name):
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(name)
    fh = logging.FileHandler(f"{name}.log")
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger

def np_mtx_to_param(prefix_param, mtx, ps, mn, mx, vary = True):
    for i in range(0, mtx.shape[0]):
        for j in range(0, mtx.shape[1]):
            ps.add(f'{prefix_param}_{i}_{j}', value= float(mtx[i][j]), min=mn, max=mx, vary = vary)
    ps[f'{prefix_param}_shape'] = Parameter(f'{prefix_param}_shape', value=float(0.0), user_data = mtx.shape, vary=False)

def np_arr_to_param(prefix_param, arr, ps, mn, mx, vary = True):
    for i in range(0, arr.shape[0]):
        ps.add(f'{prefix_param}_{i}', value= float(arr[i]), min=mn, max=mx, vary = vary)
    ps[f'{prefix_param}_shape'] = Parameter(f'{prefix_param}_shape', value=float(0.0), user_data = arr.shape, vary=False)

def mtx_from_param_to_np(prefix_param, ps):
    shape = ps[f'{prefix_param}_shape'].user_data
    param = np.zeros(shape)
    for i in range(0, shape[0]):
        for j in range(0, shape[1]):
            param[i][j] = ps[f'{prefix_param}_{i}_{j}']
    return param

def arr_from_param_to_np(prefix_param, ps):
    shape = ps[f'{prefix_param}_shape'].user_data
    param = np.zeros(shape[0])
    for i in range(0, shape[0]):
        param[i] = ps[f'{prefix_param}_{i}']
    return param

def generate_graph_from_data(alphas_best_fit, tfs, filter_genes):
    graph = nx.DiGraph(data=alphas_best_fit)
    num_cyto = alphas_best_fit.shape[0]
    nodes = []
    for t in range(0, num_cyto):
        nd = filter_genes[t]
        nodes.append(nd)
        graph.add_node(nd)
    up = []
    down = []
    edges = []
    edges_width_up = []
    edges_width_down = []
    for i in range(0, num_cyto):
        tf_idx = tfs[i]
        for j in range(0, len(tf_idx)):
            gene_idx = tf_idx[j]
            influence = alphas_best_fit[i][j] 
            if influence > 0.0:
                up.append([nodes[i],nodes[gene_idx]])
                edges_width_up.append(influence)
            if influence < 0.0:
                down.append([nodes[i],nodes[gene_idx]])
                edges_width_down.append(influence)
            
            edges.append((nodes[i], nodes[gene_idx]))


    graph.add_edges_from(edges)

    return graph, nodes, up, down