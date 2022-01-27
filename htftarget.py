import csv
import os

class hTFTargetDb:
    def __init__(self, folder_name):
        self.folder_name = folder_name
        self.tf_to_genes = {}
        self.gene_to_tfs = {}
        file_name = os.path.join(folder_name, "hTFTarget.txt")
        
        self.load(file_name)
    
    def load(self, file_name):
        return_dict = {}
        with open(file_name, encoding="utf8") as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter="\t")
            for idx, row in enumerate(tsv_reader):
                if idx == 0:
                    continue
                
                tf = row[0]
                gene = row[1]
                
                if tf not in self.tf_to_genes:
                    self.tf_to_genes[tf] = []
                
                if gene not in self.gene_to_tfs:
                    self.gene_to_tfs[gene] = []
                self.tf_to_genes[tf].append(gene)
                self.gene_to_tfs[gene].append(tf)
        return return_dict
    def get_gene_tfs(self, gene):
        if gene not in self.gene_to_tfs:
            return [] 
        return self.gene_to_tfs[gene]
   
def test():
    db = hTFTargetDb("data")
    print(db.get_gene_tfs("AEBP2"))
   