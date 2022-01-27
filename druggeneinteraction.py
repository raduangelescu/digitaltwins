import csv

class DrugGeneInteraction:
    def __init__(self, filename):
        self.filename = filename
        db_drug_genes = {}
        db_gene_drugs = {}
        DRUG_NAME_IDX = 7
        GENE_NAME_IDX = 0
        DRUG_CLAIM_NAME_IDX = 6
        with open(self.filename, encoding="utf8") as tsvfile:
            tsvreader = csv.reader(tsvfile, delimiter="\t")
            for idx, line in enumerate(tsvreader):
                if idx == 0:
                    continue
                drug = line[DRUG_NAME_IDX]
                gene = line[GENE_NAME_IDX]
                if drug.isspace() or drug == '':
                    drug = line[DRUG_CLAIM_NAME_IDX]
                if drug not in db_drug_genes:
                    db_drug_genes[drug] = []
                if gene not in db_gene_drugs:
                    db_gene_drugs[gene] = []
                db_gene_drugs[gene].append(drug)
                db_drug_genes[drug].append(gene)
        self.db_gene_drugs = db_gene_drugs
        self.db_drug_genes = db_drug_genes

    def get_all_drugs(self):
        return [*self.db_drug_genes]
        
    def get_drugs_for_gene(self, gene_name):
        if gene_name in self.db_gene_drugs:
            return self.db_gene_drugs[gene_name]
        return []

    def get_genes_for_drug(self, drug_name):
        if drug_name in self.db_drug_genes:
            return self.db_drug_genes[drug_name]
        return []