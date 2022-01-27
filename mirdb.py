import csv
import os
import requests
import pandas as pd 

class MirDB:
    def __init__(self, folder_name):
        self.folder_name = folder_name
        self.gene_to_mirna = {}
    

    def get_gene_mirnas(self, gene):
        if gene not in self.gene_to_mirna:
            return [] 
        return self.gene_to_mirna[gene]

    def get_mirnas_for_gene(self, gene_symbol, dict_gene):
        headers = {'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9',
        'Accept-Encoding':'gzip, deflate',
        'Content-Type':'application/x-www-form-urlencoded',
        'Host':'mirdb.org',
        'Origin':'http://mirdb.org',
        'Referer':'http://mirdb.org/index.html',
        'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/90.0.4430.85 Safari/537.36'
        }
        payload = {'species':'Human','geneChoice':'symbol','searchBox':f'{gene_symbol}', 'submitButton':'Go', 'searchType':'gene'}

        session = requests.Session()
        res = session.post('http://mirdb.org/cgi-bin/search.cgi',headers=headers,data=payload)
        print(f"STATUS CODE:{res.status_code}")
        #print(f"CONTENT: {res.text}")
        df_list = pd.read_html(res.text) # this parses all the tables in webpages to a list
        df = df_list[0]
        new_header = df.iloc[0]
        df = df[1:]
        df.columns = new_header
        for index, row in df.iterrows():
            if gene_symbol not in dict_gene:
                dict_gene[gene_symbol] = {}
            mirna = row['miRNA Name']   
            dict_gene[gene_symbol][mirna] = {
                'score': row['Target Score'],
                'miRNA': mirna
            }

    def get_mirnas_for_all_genes(self, gene_list):
        gene_dict = {}
        for gene_symbol in gene_list:
            self.get_mirnas_for_gene(gene_symbol, gene_dict)
        return gene_dict

