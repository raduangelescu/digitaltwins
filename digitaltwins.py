import numpy as np
import json
import os
import pickle
import matplotlib.pyplot as plt
import scipy
from utils import get_logger
from lincsregulatedb import LincsRegulateDb
from htftarget import hTFTargetDb
from mirdb import MirDB
from dataimport import DataImport
import model
import model_drug
from lmfit import minimize, Parameters, Parameter, report_fit
from pathlib import Path
from utils import generate_graph_from_data, arr_from_param_to_np
import networkx as nx
import itertools

class DigitalTwinsConfig:
    def __init__(self,
                 data_folder):
        self.data_folder = data_folder

def min_progress(params, iter, resid, *args, **kws):
    sum = np.sum(resid)
    logger = get_logger('DigitalTwins')
    logger.info(f"iter: {iter}, totalerror: {sum}")
    if iter % 100 == 0:
        logger.info(f"params {params}")

class DigitalTwins:
    def __init__(self, config_filename):
        with open(config_filename, ) as f:
            config_json = json.load(f)
        self.logger = get_logger('DigitalTwins')
        self.config = DigitalTwinsConfig(**config_json)
        lfile = os.path.join(self.config.data_folder, "LINCS_L1000")
        self.tfdb =  hTFTargetDb("data")
        self.mirdb = MirDB("data")
        self.lincs = LincsRegulateDb(lfile)
        self.data_source = DataImport()
    
    def get_first_k_variance_genes(self, data, gene_names, k):
        gene_variances = np.var(data, axis=0)
        list_idx = np.argsort(gene_variances).tolist()[-k:]
        variances = gene_variances[list_idx]
        self.logger.info(f"sorted gene variances {variances}")
        filtered_gene_data = data[:, list_idx]
        filtered_gene_names = np.array(gene_names)[list_idx].tolist()
        return filtered_gene_names, filtered_gene_data

    def generate_variance_histogram(self, data):
        gene_variances = np.var(data, axis=0)
        plt.hist(gene_variances, bins = [0,0.11,0.12,0.13,0.14,0.15,0.152, 0.153, 0.155,0.156,0.157,0.158,0.16,0.17,0.18,0.19, 0.2]) 
        plt.title("Gene expression histogram") 
        plt.savefig("hist.png")
    
    def compute_tfs(self, gene_names):
        the_tfs = []
        inner_idx = {}
        for idx, gene in enumerate(gene_names):
            inner_idx[gene] = idx
        for gene in gene_names:
            tfs = self.tfdb.get_gene_tfs(gene)
            np_tfs = []
            for tf in tfs:
                if tf in inner_idx:
                    np_tfs.append(inner_idx[tf])
            np_tfs = np.array(np_tfs)
            the_tfs.append(np_tfs)
        return the_tfs
    
    def expand_tf(self, fsens_genes, all_data, max_num_tf):
        check_fsens_genes = set(fsens_genes)
        for gene in fsens_genes:
            gene_data = all_data[:, self.gene_index[gene]]
            tfs = self.tfdb.get_gene_tfs(gene)
            tf_rank = []
            tf_name = []
            for tf in tfs:
                if tf in check_fsens_genes:
                    continue
                if tf not in self.gene_index:
                    continue
                tf_data = all_data[:, self.gene_index[tf]]
                tf_rank.append(np.abs(np.corrcoef(gene_data, tf_data)[1,0]))
                tf_name.append(tf)
            tf_rank = np.array(tf_rank)
            tf_name = np.array(tf_name)
            tfs_idx = np.argsort(tf_rank)[-max_num_tf:]
            tf_name = tf_name[tfs_idx]
            for tf in tf_name:
                check_fsens_genes.add(tf)
        names =[]
        data = []
        for gene in check_fsens_genes:
            names.append(gene)
            data.append(all_data[:, self.gene_index[gene]])
        return names, np.transpose(np.array(data)) 

    def normalize_data(self, data):
        data = np.array(data)
        max_raw_data = data.max()
        return  data/max_raw_data
    
    def populate_gene_idx_map(self, gene_names):
        self.all_genes = gene_names
        self.gene_index = {}
        for index, gene in enumerate(self.all_genes):
            self.gene_index[gene] = index

    def interpolate_data(self, data, interpolation_type):
        real_times = np.array([0, 6, 12, 24, 48])
        t_sim = np.linspace(real_times[0], real_times[-1], 100)
        interpolated_data = []
        if interpolation_type == "hermite":
            # cubic hermite
            coeff = np.polynomial.hermite.hermfit(real_times,data,3)
            interpolated_data = np.polynomial.hermite.hermval(t_sim, coeff)
            return interpolated_data, t_sim
        for i in range(data.shape[1]):
            interp_data = data[:,i]
            if interpolation_type == "linear":
                f_i = scipy.interpolate.interp1d(real_times, interp_data)
            elif interpolation_type == "pchip":
                f_i = scipy.interpolate.PchipInterpolator(real_times, interp_data)
            elif interpolation_type == "akima":
                f_i = scipy.interpolate.Akima1DInterpolator(real_times, interp_data)
            elif interpolation_type == "cubic":
                f_i = scipy.interpolate.CubicSpline(real_times, interp_data)
            else:
                self.logger.error("CANNOT UNDERSTAND INTERPOLATION")
            interpolated_data.append(f_i(t_sim))
        interpolated_data = np.array(interpolated_data)
        interpolated_data = np.nan_to_num(interpolated_data)
        return interpolated_data, t_sim
    
    def get_model_genes(self, gene_names, norm_data, model_params):
        genes, data = self.get_first_k_variance_genes(norm_data, 
                                                    gene_names, 
                                                    model_params['num_main_genes'])
        genes, data = self.expand_tf(genes, norm_data, model_params['max_number_tf_factors'])
        sort_genes_idx = np.argsort(np.array(genes))
        genes = np.array(genes)[sort_genes_idx].tolist()
        data = data[:, sort_genes_idx]
        return genes, data
    
    def get_pkl_filename(self, model_params):
        params_suffix = f"i{model_params['interpolation_type']}_tf{model_params['max_number_tf_factors']}_mg{model_params['num_main_genes']}"
        output = os.path.join("results",model_params['name'], f"best_fit_result_{params_suffix}.pkl")
        return output
        
    def get_startplot_filename(self, model_params):
        params_suffix = f"i{model_params['interpolation_type']}_tf{model_params['max_number_tf_factors']}_mg{model_params['num_main_genes']}"
        start_plot_filename = os.path.join("results",model_params['name'], f"start_{params_suffix}.png")
        return start_plot_filename
    
    def get_fitplot_filename(self, model_params):
        params_suffix = f"i{model_params['interpolation_type']}_tf{model_params['max_number_tf_factors']}_mg{model_params['num_main_genes']}"
        output_plot = os.path.join("results",model_params['name'], f"best_fit_result_{params_suffix}.png")
        return output_plot
    
    def load_model(self, filename):
        best_fit = None
        with open(filename, 'rb') as fp:
            best_fit = pickle.load(fp)
        return best_fit
    
    def compute_drug_rank(self, effect):
        return np.count_nonzero(effect)

    def apply_drug(self, dgs_sig, drug):
        return dgs_sig + drug
    
    def search_mirna(self, genes, DGS_sig):
        search_genes = []
        for i in range(0, DGS_sig.shape[0]):
            if DGS_sig[i] == -1:
                search_genes.append(genes[i])
        # score mirnas by how many genes they down regulate
        mirdb_ranks = {}
        mir = self.mirdb.get_mirnas_for_all_genes(search_genes)
        for idx, gene in enumerate(search_genes):
            if gene not in mir:
                continue
            mirnas = mir[gene]
            for mirna in mirnas:
                if mirna not in mirdb_ranks:
                    mirdb_ranks[mirna] = 1
                else:
                    mirdb_ranks[mirna] = mirdb_ranks[mirna] + 1
        # sort mirnas by downregulation
        best_mirna = {'score':-1, 'name':''}
        for mirna in mirdb_ranks:
            if mirdb_ranks[mirna] > best_mirna['score']:
                best_mirna['score'] = mirdb_ranks[mirna]
                best_mirna['name'] = mirna
        return best_mirna 

    def search_drug(self, drugs, genes, DGS_sig, drug_test_num):
        drug_reps = []
        best_drug_combination_effect = []
        #create drug representations
        for idx, drug in enumerate(drugs):
            gene_inf_dict = self.lincs.get_genes_drug(drug)
            dg_inf = list(map(lambda g: gene_inf_dict[g] if g in gene_inf_dict else 0, genes))
            dg_inf = np.array(dg_inf)
            drug_reps.append([dg_inf, idx, drug])
        
        best_drug_combination = {'rank': 99999, 'details':[]} 

        all_tests = list(itertools.combinations(drug_reps, drug_test_num))
        #start search
        drug_combination_effect = DGS_sig
        total_num_of_tests = len(all_tests)
        for idx, drug_combination in enumerate(all_tests):
            if idx % 100000 == 0:
                self.logger.info(f"percentage done: {100 * (idx/total_num_of_tests)}")
            for drug in drug_combination:
                drug_combination_effect = self.apply_drug(drug[0], DGS_sig)
            #normalize effect
            drug_combination_effect = np.clip(drug_combination_effect, -1, 1)
            #compute rank 
            drug_combination_rank = float(self.compute_drug_rank(drug_combination_effect))
            if best_drug_combination['rank'] > drug_combination_rank:
                best_drug_combination['rank'] = drug_combination_rank
                best_drug_combination['details'] = drug_combination
                best_drug_combination['idx'] = idx
                best_drug_combination['drug'] = drug
                best_drug_combination['genes'] = genes
                best_drug_combination_effect = drug_combination_effect

        self.logger.info(f"best drug combination {best_drug_combination}")
        return best_drug_combination, best_drug_combination_effect

    def diff_model(self, 
                data_sens, 
                data_res, 
                genes, 
                diff_params, 
                search_2_drug_combinations=False):
     
        res_model_filename = diff_params['res_from_sens_file']
        sens_model_filename = diff_params['sens_from_sens_file']
        #interpolate data 
        interp_data_sens, interp_time = self.interpolate_data(data_sens, diff_params['interpolation_type'])
        interp_data_res, interp_time = self.interpolate_data(data_res, diff_params['interpolation_type'])
        tfs = self.compute_tfs(genes)
        
        res_fit = self.load_model(res_model_filename)
        sens_fit = self.load_model(sens_model_filename)
        
        self.logger.info(f"Number of function evaluations: {res_fit.nfev}")
        self.logger.info(f"Number of variables in fit: {res_fit.nvarys}")
        self.logger.info(f"Number of data points: {res_fit.ndata}")
        self.logger.info(f"Number of degrees of freedom in fit: {res_fit.nfree}")
        self.logger.info(f"ChiSqr: {res_fit.chisqr}")
        self.logger.info(f"Reduced ChiSqr: {res_fit.redchi}")
        self.logger.info(f"Akaike Information Criterion statistic: {res_fit.aic}")
        self.logger.info(f"Bayesian Information Criterion statistic: {res_fit.bic}")

        #recompute sens data
        model_sens_data = interp_data_sens + sens_fit.residual.reshape(interp_data_sens.shape)
        model_res_data = interp_data_res + res_fit.residual.reshape(interp_data_res.shape)
        alphas_best_fit = []
        for i, x_i in enumerate(genes):
            alphas_p_i = arr_from_param_to_np(f'alpha__{i}', res_fit.params)
            alphas_best_fit.append(alphas_p_i)
        alphas_best_fit = np.array(alphas_best_fit)
        graph, nodes, up, down = generate_graph_from_data(alphas_best_fit, tfs, genes)
        pos = nx.circular_layout(graph)
        plot_colors = [plt.cm.tab20(i) for i in range(0,len(genes))]
        nx.draw_networkx_nodes(graph, pos, node_color = plot_colors)
        nx.draw_networkx_labels(graph, pos)
        nx.draw_networkx_edges(graph, pos,edgelist=up, arrowstyle='-|>', edge_color='green', alpha=0.5,
                                connectionstyle='arc3,rad=0.1', min_source_margin=20, min_target_margin=20)#, width=edges_width_up)
        nx.draw_networkx_edges(graph, pos, edgelist=down, arrowstyle='-[', arrowsize=5, edge_color='red', alpha=0.5,
        connectionstyle='arc3,rad=0.1', min_source_margin=20, min_target_margin=20)#, width=edges_width_down)

        plt.savefig("network.png")

        rest_best_fit_params = res_fit.params
        
        drugs = self.lincs.get_drugs_with_genes(genes)
       
        sq = model_sens_data - model_res_data
        sq = np.multiply(sq, sq)
        base_line = np.sum(sq)
        num_tests = 50
        for i in range(0, num_tests):
            params_suffix = f"ilinear_tf3_mg5_idx{i}"
            output = os.path.join("results", f"cure_result_{params_suffix}.pkl")
            if os.path.isfile(output):
                continue
            dg_inf = []
            params = model_drug.init_params_sim(genes, interp_data_res, tfs, rest_best_fit_params, dg_inf)
            
            #actually fit data
            result = minimize(model_drug.residual, params, args=(interp_time, interp_data_sens), method='leastsq', iter_cb= min_progress)
          
            with open(output, 'wb') as fp:
                pickle.dump(result, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
        #evaluate fits
        rankings =  []
        self.logger.info("__________________________________________________")
        self.logger.info(f"no drug rank (sens - ref distance): {base_line}")
        for i in range(0, num_tests):
            params_suffix = f"ilinear_tf3_mg5_idx{i}"
            fit_cure_name = os.path.join("results", f"cure_result_{params_suffix}.pkl")
            
            fit_cure = self.load_model(fit_cure_name) 
            sq = fit_cure.residual.reshape(interp_data_sens.shape)
            sq = np.multiply(sq, sq)
            rank_fit_cure = np.sum(sq)
            rankings.append({'idx':i,'rank':rank_fit_cure})
        sorted_rankings = sorted(rankings, key=lambda x: x['rank'], reverse=False)
        self.logger.info(f"sorted rankings: {sorted_rankings}")
        top_model_idx = sorted_rankings[0]['idx']
        params_suffix = f"ilinear_tf3_mg5_idx{top_model_idx}"
        fit_cure_name = os.path.join("results", f"cure_result_{params_suffix}.pkl")
        fit_cure = self.load_model(fit_cure_name) 
        Z = arr_from_param_to_np("drug_gene_inf", fit_cure.params)
        b = fit_cure.params['b'].value
        tau = fit_cure.params['tau'].value
        self.logger.info(f"drug gene influence: {Z}")
        self.logger.info(f"tau: {tau}")
        self.logger.info(f"b: {b}")
        drug_evolution = b * np.exp(-interp_time/tau)
        plt.clf()
        plt.plot(interp_time, drug_evolution)
        plt.ylabel("Drug Concentration")
        plt.xlabel("Time")
        plt.savefig("drug_evolution.png")
        
        binarizer = lambda t: -1 if t<0 else 1
        DGS_sig_func = np.vectorize(binarizer)
        DGS_sig = DGS_sig_func(Z)
        self.logger.info("searching single drug")
        best_drug, residual = self.search_drug(drugs, genes, DGS_sig, 1)
        self.logger.info(f"best single drug {best_drug}")
        self.logger.info(f"residual {residual}")
        
        if search_2_drug_combinations:
            self.logger.info("searching 2 drug combinations")
            best_drug = self.search_drug(drugs, genes, DGS_sig, 2)
            self.logger.info(f"best two drug combination {best_drug}")
        
        best_mirna = self.search_mirna(genes, residual)
        self.logger.info(f"best mirna {best_mirna}")
        
        
        self.logger.info("done")

    def fit_model(self, data, genes, model_params):
        #create folders if not present
        Path(os.path.join("results",model_params['name'])).mkdir(parents=True, exist_ok=True)
        #interpolate data 
        interp_data, interp_time = self.interpolate_data(data, model_params['interpolation_type'])
        
        #plot interpolated data
        num_genes = len(genes)
        plot_colors = [plt.cm.tab20(i) for i in range(0,num_genes)]
        figure = plt.gcf()
        figure.set_size_inches(10, 10)
        for i in range(0, num_genes):
            plt.plot(interp_time, interp_data[i, :], '-', color=plot_colors[i], label=genes[i])
        plt.legend(shadow=True, fancybox=True)
        start_plot_filename = self.get_startplot_filename(model_params)
        plt.savefig(start_plot_filename, dpi=300)
        
        #compute transcription factors links from within the gene pool
        tfs = self.compute_tfs(genes)
        
        params = model.init_params(genes,
                            interp_data,
                            tfs)
        
        #actually fit data
        result = minimize(model.residual, params, args=(interp_time, interp_data), method='leastsq', iter_cb= min_progress)
        
        #compute our model result
        final = interp_data + result.residual.reshape(interp_data.shape)
        self.logger.info("done fitting")
        #save results
        output = self.get_pkl_filename(model_params)
        with open(output, 'wb') as fp:
            pickle.dump(result, fp, protocol=pickle.HIGHEST_PROTOCOL)
        
        self.logger.info("wrote results to file")
        plt.clf()
        output_plot = self.get_fitplot_filename(model_params)
        for i in range(0, len(genes)):
            plt.plot(interp_time, interp_data[i, :], 'o', color=plot_colors[i])
            plt.plot(interp_time, final[i, :],'-', label=genes[i], color=plot_colors[i], linewidth=2)
        plt.legend(shadow=True, fancybox=True)
        plt.savefig(output_plot, dpi=300)

    def get_data_from_gene_names(self, genes, data):
        ret_data = []
        for gene in genes:
            gene_data = data[:, self.gene_index[gene]]
            ret_data.append(gene_data)
        return np.transpose(np.array(ret_data))

    def run(self, experiment_name, model_params, force_fit = False):
        pkl_fit_filename = self.get_pkl_filename(model_params)
        gene_names, sensitive, resistant, unknown = self.data_source.get_data_and_filter(experiment_name)
        self.populate_gene_idx_map(gene_names)
        norm_data_sens = self.normalize_data(sensitive)
        norm_data_res = self.normalize_data(resistant)
        sens_genes, sens_data = self.get_model_genes(gene_names, norm_data_sens, model_params)
        
        if os.path.isfile(pkl_fit_filename) == False or force_fit:
            self.logger.info("could not find result of fit, doing fit")
            
            if(model_params['name'] == "sens_from_sens"): 
                self.fit_model(sens_data, sens_genes, model_params)
            else:
                self.get_model_genes(gene_names, norm_data_res, model_params)
                res_from_sens = self.get_data_from_gene_names(sens_genes, norm_data_res)
                self.fit_model(res_from_sens, sens_genes, model_params)
            self.logger.info("finished fitting model")
        
        self.logger.info("testing hypothesis")
        sens_genes, sens_data = self.get_model_genes(gene_names, norm_data_sens, model_params)
        res_from_sens = self.get_data_from_gene_names(sens_genes, norm_data_res)
    
        diff_params = {
            'interpolation_type': model_params['interpolation_type'],
            'res_from_sens_file': os.path.join("results","res_from_sens","best_fit_result_ilinear_tf3_mg5.pkl"),
            'sens_from_sens_file': os.path.join("results","sens_from_sens","best_fit_result_ilinear_tf3_mg5.pkl")
        }

        self.diff_model( sens_data, res_from_sens, sens_genes, diff_params)        
        
model_params = {
            'interpolation_type': "linear",
            'num_main_genes': 5,
            "max_number_tf_factors": 3,
            "name": "res_from_sens"
        }

twins = DigitalTwins('data/config.json')
twins.run('GSE128722', model_params, force_fit = False)
