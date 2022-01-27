import csv
import os

class LincsRegulateDb:
    def __init__(self, folder_name):
        self.folder_name = folder_name
        self.db = {}
        down_file = os.path.join(folder_name, "gene_set_library_dn_crisp.gmt")
        up_file = os.path.join(folder_name, "gene_set_library_up_crisp.gmt")
        self.db["down"] = self.load_gmt(down_file)
        self.db["up"] = self.load_gmt(up_file)
        self.db["reverse"] = self.load_reverse_db()
    
    def _add_perturb(self, db, type):
        for drug, genes in self.db[type].items():
            for gene in genes:
                if gene not in db:
                    db[gene] = {}
                else:
                    if type not in db[gene]:
                        db[gene][type] = [drug]
                    else:
                        db[gene][type].append(drug)

    def load_reverse_db(self,):
        db = {}
        self._add_perturb(db, "down")
        self._add_perturb(db, "up")
        return db

    def load_gmt(self, file_name):
        return_dict = {}
        self.drugs = set([])
        with open(file_name, encoding="utf8") as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter="\t")
            for idx, row in enumerate(tsv_reader):
                drug = row[0]
                # row[1] is the description, we don't care about it
                genes = row[2:]
                return_dict[drug] = genes
                self.drugs.add(drug) 
        return return_dict
    
    def get_genes_drug(self, drug):
        gene_dict = {}
        for gn in self.db["up"][drug]:
            gene_dict[gn] = 1
        for gn in self.db["down"][drug]:
            gene_dict[gn] = -1
        return gene_dict

    def get_drugs_with_genes(self, genes):
        ret = []
        drgs = list(self.drugs)
        for drg in drgs:
            drg_gns = set()
            if drg in self.db["up"]:
                lst = self.db["up"][drg]
                drg_gns = drg_gns.union(set(lst))
            if drg in self.db["down"]:
                lst = self.db["down"][drg]
                drg_gns = drg_gns.union(set(lst))
            gset = set(genes)
            inter = gset.intersection(drg_gns)
            if len(inter) != 0:
                ret.append(drg)
        return ret

    def get_perturb_drugs(self, gene_name, perturbation_type):
        if perturbation_type != 'down' and perturbation_type != "up":
            return []
        search_db = self.db["reverse"]
        
        if gene_name not in search_db:
            return []
        gene_set = search_db[gene_name]
        
        if perturbation_type not in gene_set:
            return []
        
        return gene_set[perturbation_type]
        
def test():
    db = LincsRegulateDb("data/LINCS_L1000")
    gene_name_up = 'COL1A2'
    gene_name_down = 'GATM'
    up = db.get_perturb_drugs(gene_name_up, "up")
    down = db.get_perturb_drugs(gene_name_down, "down")
    print(f"drugs that perturbe the gene {gene_name_down} down {down}")
    print("--------------------")
    print(f"drugs that perturbe the gene {gene_name_up} up {up}")
