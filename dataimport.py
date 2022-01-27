import numpy as np
import pandas as pd
import numpy as np
import os
import pickle

class DataImport:
    def __init__(self):
        self.folder = 'data'
    
    def _get_data_from_XLSX(self, filename):
        read_file = os.path.join(self.folder, filename)
        df = pd.read_excel(read_file)
        gene_names = df['Gene'].to_list()
        # 0 6 12 24 48 h
        DBTRG = self._get_from_prefix(df, 'D',1, 5) #sensitive cells
        LN18 = self._get_from_prefix(df, 'L',1, 5) #resistant cells
        U87MG = self._get_from_prefix(df, 'U',1, 5) #unknown
        return gene_names, DBTRG, LN18, U87MG

    