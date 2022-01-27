# Using Digital twins to fight drug resistance support material
These are all the scripts used in the paper "Using Digital twins to fight drug resistance"
To reproduce our experiments you need install Python 3.9.7.

## How to use
- Open up a terminal in your prefered folder and run:

    `git clone git@github.com:raduangelescu/digitaltwins.git .`
- Install all python requirements by running (you may use a virtual environment): 
    
    `pip install -r requirements.txt`
- Run the startup script that downloads and unpacks experiment data:

    `python startup.py`

- To run all experiments you just need to:

    `python digitaltwins.py`

The output will reside in the .log file and in the results folder. Note that the 2-drug combination searcher is disabled by default because it is computationally costly, to change that you just need to add the True parameter at the diff_model function call:
    
        self.diff_model( sens_data, res_from_sens, sens_genes, diff_params, **True**) 

