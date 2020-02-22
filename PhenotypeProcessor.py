
from concurrent.futures import ProcessPoolExecutor
import os
from SodalisEvolver import save_phenotypes
import numpy as np
import collections
import csv
import multiprocessing



if __name__ == "__main__":
    results_directories =['iRHx_Blood']
    models = ['iRHx_blood.json']
    # test_results_directories = ['blood','iRH826_glucose']
    # test_models = ['iRH826_blood_v3.json','iRH826_v3.json']

    with ProcessPoolExecutor() as executor:
        for path, run in zip(results_directories, executor.map(save_phenotypes, results_directories,models)):
            pass
