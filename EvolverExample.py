
import numpy as np
import time as time
from concurrent.futures import ProcessPoolExecutor
from SodalisEvolver import CobraFBAEvolver
import os
from itertools import repeat
import sys

'''
Example of running the CobraFBAEvolver on mulitple cores
'''

'''
The function to be parallelised 

Arguments:
results_directory--filepath where the results will be saved
model-- filepath to cobrapy json FBA model
num_gens -- number of generations the algorithm is to be run for
knockout-- Reaction name to be knocked out, default no reactions to be knocked

'''
def work(results_directory,model, num_gens, knockout):

    gen_record = range(0, num_gens, 50)

    # print('Beginning FBA Evolution for ' + str(num_runs[-1]) + ' runs of ' + str(gen_record[-1]) + ' generations')
    print(results_directory)
    num_run = int(results_directory.split('/')[-1])
    print(model)
    pop = None
    for gens in gen_record:
        evolver = CobraFBAEvolver(model, results_directory,knockout=knockout)
        pop, log = evolver.run_nsga2evo(num_run, gens, population=pop,ngen=51)


if __name__ == "__main__":  # Allows for the safe importing of the main module
    # print("There are %d CPUs on this machine" % multiprocessing.cpu_count())

    results_directory  = sys.argv[1]

    model = sys.argv[2]
    model_path = os.path.join('../models',model)

    #Name of reaction knockout,
    if len(sys.argv) > 3:
        knockout = sys.argv[3]
    else:
        knockout= False
    times = list()

    #Number of cores to be used
    cores = 10
    cores_list = np.arange(1,cores,1)

    #Full results path
    results_path = os.path.join('../results',results_directory)

    results_paths = list()

    #Construct results directory paths for passing to each worker process
    for core in cores_list:
        try:
            os.makedirs(os.path.join(results_path,str(core)))
        except FileExistsError:
            print('File Exists')
        results_paths.append(os.path.join(results_path,str(core)))


    #number_processes = cores
    para_start_time = time.time()
    with ProcessPoolExecutor() as executor:
        #For each results directory create a worker process
        for path, run in zip(results_paths, executor.map(work, results_paths, repeat(model_path),repeat(knockout))):
            pass

    finish_time = time.time()
    execution_time = finish_time - para_start_time
    print('Execution Time: ' + str(execution_time))

