
from concurrent.futures import ProcessPoolExecutor
import os
import numpy as np
import cobra
import random
import os
import multiprocessing
from itertools import repeat
from deap import base
from deap import creator
from deap import tools
import collections
import csv
import multiprocessing

"""
Return list of named tuple Phenotype collections
growth-- biomass output of model
used_reactions -- list of reaction names of the active non essential reactions
removed_reactions -- list of reaction names of the inactive non essential reactions
essential_reactions -- list of reaction names deemed essential (biomass < 0.001 when knocked out)

Arguments:
filename -- full file path to results .csv file
evolver -- FBAEvolver object
"""
def get_phenotypes(filename,evolver):
    with open(filename, 'r') as csv_file:

        phenotypes = list()
        reader = csv.reader(csv_file, delimiter=',')
        for row in reader:
            removed_reactions = list()
            used_reactions = list()
            essential_reactions = list()
            int_list = [int(num) for num in row]
            queue = multiprocessing.Queue()
            growth, reactions = evolver.fitness_function(int_list,queue)
            for reaction, ind in zip(evolver.fba.non_essential_reactions, int_list):
                if ind:
                    used_reactions.append(reaction.id)
                else:
                    removed_reactions.append(reaction.id)
            for reaction in evolver.fba.essential_reactions:
                essential_reactions.append(reaction.id)

            Phenotype = collections.namedtuple('Phenotype', ['growth','used_reactions', 'removed_reactions',
                                                             'essential_reactions'])
            phenotypes.append(Phenotype(growth,used_reactions,removed_reactions,essential_reactions))
    return phenotypes

"""
Get the list of generations that have been saved for each repeat of the algorithm

Arugments:
results_directory -- filepath to the FBAEvolver results
num_gens: maxmium number of generations the model has been evolved for
"""
def collect_generation_file(results_directory,num_gens):

    generations = range(0, num_gens, 50)
    generations_list = list()

    for generation in generations:
        file_list = list()
        for directory in os.listdir(results_directory):
            try:
                print(directory)

                print(os.path.join(results_directory, directory))
                for file in os.listdir(os.path.join(results_directory, directory)):
                    print(file)
                    filename = file.split('.')[0]
                    gen_number = filename.split('_')[-1]
                    print(gen_number)
                    if generation == int(gen_number):
                        file_list.append(os.path.join(results_directory, directory, file))

            except NotADirectoryError:
                pass
        generations_list.append(file_list)
        print(generations_list)

"""
Convert the raw binary data output from the FBAEvolver to numpyz files 
"""

def save_phenotypes(results,model):
    results_directory = os.path.join('../results',results)
    model_path = os.path.join('../models', model)
    evolver = CobraFBAEvolver(model_path, '/results')

    generations_list, generations = collect_generation_file(results_directory)

    for files, generation in zip(generations_list, generations):
        file_phenotypes = list()
        growth = list()
        used_reactions = list()
        used_reaction_names = list()
        for file in files:
            file_phenotypes.append(get_phenotypes(file, evolver))
        phenotypes = [item for sublist in file_phenotypes for item in sublist]
        for phenotype in phenotypes:
            growth.append(phenotype.growth)
            used_reactions.append(len(phenotype.used_reactions))
            used_reaction_names.append(phenotype.used_reactions)
        np.savez(os.path.join(results_directory, str(generation)), growth=growth, used_reactions=used_reactions,used_reaction_names=used_reaction_names)

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]
'''
CobraFBA class contains the cobrafba model and evaluates the lists of essential and non essential reactions
'''
class CobraFBA():
    def __init__(self, input, knockout=False):

        # import iaf1260
        self.model = cobra.io.load_json_model(input)
        print(self.model.objective)
        self.knockout = knockout
        self.model.objective = "Biomass"
        if knockout:
            # print('Knocked Out: ' + str(knockout))
            self.model.reactions.get_by_id(knockout).knock_out()
        self.non_essential_reactions, self.essential_reactions, self.exchanges, self.exchange_bounds, self.modifiable_exchanges = self.evaluate_essential_reactions(
            self.model)
        self.reaction_names, self.initial_reaction_bounds = self.find_reaction_intial_bounds(
            self.non_essential_reactions)
        self.essential_reaction_names, essential_reaction_bounds = self.find_reaction_intial_bounds(
            self.essential_reactions)

    def find_reaction_intial_bounds(self, reactions):
        reaction_names = list()
        exchange_names = list()
        exchange_bounds = list()
        reaction_bounds = list()
        for reaction in reactions:
            if 'EX' in reaction.id[:2]:
                exchange_names.append(reaction.id)
                exchange_bounds.append(reaction.upper_bound)

            else:
                reaction_names.append(reaction.id)
                reaction_bounds.append([reaction.lower_bound, reaction.upper_bound])

        return reaction_names, dict(zip(reaction_names, reaction_bounds))

    def set_reaction_bounds(self, individual):

        for reaction_active, reaction in zip(individual, self.reaction_names):
            if reaction_active:
                bounds = self.initial_reaction_bounds[reaction]
                self.model.reactions.get_by_id(reaction).lower_bound = bounds[0]
                self.model.reactions.get_by_id(reaction).upper_bound = bounds[1]
            else:
                self.model.reactions.get_by_id(reaction).lower_bound = 0
                self.model.reactions.get_by_id(reaction).upper_bound = 0

    def set_modifiable_exchanges(self, individual):
        for exchange_active, exchange in zip(individual, self.modifiable_exchanges):
            if exchange_active:
                self.model.reactions.get_by_id(exchange.id).lower_bound = -1000
            else:
                self.model.reactions.get_by_id(exchange.id).lower_bound = 0

    def run_fba(self, queue=False):
        growth = self.model.slim_optimize()
        if self.knockout:
            self.model.reactions.get_by_id(self.knockout).knock_out()
        if queue:
            queue.put(growth)
        return growth

    def evaluate_essential_reactions(self, model):
        reactions_list = list()
        reaction_bounds = list()
        exchange_list = list()
        modifiable_exchange_list = list()
        exchange_bounds = list()

        for reaction in model.reactions:
            #Add all exchange reactions to exchange list
            if 'EX' in reaction.id[:2]:
                exchange_list.append(reaction)
                exchange_bounds.append(reaction.lower_bound)

                if reaction.lower_bound != 0:
                    pass
                else:
                    modifiable_exchange_list.append(reaction)
            else:
                if np.abs(reaction.lower_bound) + reaction.upper_bound != 0:
                    reactions_list.append(reaction)
                    reaction_bounds.append([reaction.lower_bound, reaction.upper_bound])


        non_essential_reactions = list()
        essential_reactions = list()

        #For each reaction test if growth is greater than 0.001, if so reaction is non essential
        for reaction, reaction_bound in zip(reactions_list, reaction_bounds):

            model.reactions.get_by_id(reaction.id).lower_bound = 0
            model.reactions.get_by_id(reaction.id).upper_bound = 0
            knockout_growth = model.slim_optimize()

            if knockout_growth > 0.001:
                non_essential_reactions.append(reaction)
            else:
                essential_reactions.append(reaction)

            model.reactions.get_by_id(reaction.id).lower_bound = reaction_bound[0]
            model.reactions.get_by_id(reaction.id).upper_bound = reaction_bound[1]

        print('Number of essential reactions: ' + str(len(essential_reactions)))
        print('Number of non essential reactions: ' + str(len(non_essential_reactions)))
        print('Number of unset reactions: ' + str(len(modifiable_exchange_list)))

        return non_essential_reactions, essential_reactions, exchange_list, dict(
            zip(exchange_list, exchange_bounds)), modifiable_exchange_list

    def run_full_fba(self):
        if self.knockout:
            self.model.reactions.get_by_id(self.knockout).knock_out()
        solution = self.model.optimize()
        return solution

'''
Base class for FBAEvolver classes delcaring the evolutionary algorithm and some basic figntess functions
'''
class CobraFBABase():
    def __init__(self):
        self.fba = None

    def create_evo(self, num_reactions):
        creator.create("FitnessMin", base.Fitness, weights=(1, -1))
        creator.create("Individual", list, fitness=creator.FitnessMin)
        self.toolbox = base.Toolbox()
        # Attribute generator
        self.toolbox.register("attr_bool", random.randint, 1, 1)
        # Structure initializers
        self.toolbox.register("individual", tools.initRepeat, creator.Individual,
                              self.toolbox.attr_bool, num_reactions)
        # pool = multiprocessing.Pool()
        # self.toolbox.register("map", pool.starmap)

        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)
        self.toolbox.register("individual_restart", self.initInd, creator.Individual)
        self.toolbox.register("population_restart", self.initPop, list, self.toolbox.individual_restart,
                              self.population_file)

        self.toolbox.register("mate", tools.cxTwoPoint)
        self.toolbox.register("mutate", tools.mutFlipBit, indpb=0.0005)
        self.toolbox.register("select", tools.selNSGA2)
        self.toolbox.register("HallOfFame", tools.HallOfFame)

    def fitness_function_full(self, individual, fba):
        fba.set_reaction_bounds(individual)
        solution = fba.run_full_fba()
        return solution

    def fitness_function_no_queue(self, individual):
        self.fba.set_reaction_bounds(individual)
        growth = self.fba.run_fba()

        if growth > 0.001:
            reactions = sum(individual)
        else:
            reactions = len(self.fba.non_essential_reactions)
        return growth, reactions

    '''
    Main nsga2 fitness function which maxmimises the biomass output while minimising the number of active non essential
    reactions. The use of multiprocessing.Process is a crude fix for a bug which caused the algorithm to hang after 
    >500 generations. Specifically the algorithm hangs during solving an evolved FBA model at a random time number of generations
    The model in question can be successfully solved if saved and run again separately so there is nothing intrinsically
    wrong with the model.
     
    A solution to this is to check if the FBA model takes more than 3 seconds to converge on a solution it kills the process, and sets
    the biomass output to 0 and maximises the number of reactions to try to remove this solution from the population.
    As far as I'm aware from multiple logs however a process has never been killed by this code so it makes me think that 
    there is some kind of threading conflict that is generated by cobrapy. As I'm unable to consistently reproduce this bug
    I did not raise this with cobrapy and instead continued with this fix. If anyone has any thoughts on why this might be 
    happening I would very much appreciate you getting in touch!
    
    Arguments:
    inidividual -- member of the closne population to be evaluated
    queue -- python Queue object to get the output of the python process used for the aforementioned bug
    
    
    '''
    def fitness_function(self, individual, queue):

        self.fba.set_reaction_bounds(individual)
        p = multiprocessing.Process(target=self.fba.run_fba, args=(queue,))
        p.start()
        p.join(3)

        if p.is_alive():
            # print('This is taking too long...')
            p.terminate()
            p.join()
            growth = 0
        else:
            growth = queue.get()

        if growth > 0.001:
            reactions = sum(individual)
        else:
            reactions = len(self.fba.non_essential_reactions)
        return growth, reactions

    '''
    Fitness function for generating minimal media
    '''
    def exchange_fitness_function(self, individual):

        with open('current_ind.csv', 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(individual)

        self.fba.set_modifiable_exchanges(individual)

        growth = self.fba.run_fba()
        if growth > 0.001:
            reactions = sum(individual)
        else:
            reactions = len(self.fba.modifiable_exchanges)
        return growth, reactions

    def initInd(self, icls, content):
        return icls(content)

    def initPop(self, pcls, ind_init, filename):
        with open(filename, 'r') as csv_file:
            initial_list = list()
            reader = csv.reader(csv_file, delimiter=',')
            for ind in reader:
                if ind:
                    initial_list.append([int(i) for i in ind])
            return pcls(ind_init(ind) for ind in initial_list)



'''
An evolver class that calculates the fitness of a model run in two conditions
The fitness function is g1*t + g2*(1-t) where g1 is the growth in media 1, 
g2 the growth in media 2 and t is the number of proportion of generations in media 1
'''
class JamieFBAEvolver(CobraFBABase):
    def __init__(self, model_1, model_2, results_folder, knockout=False):
        super().__init__()
        self.population_file = 'false'
        self.results_folder = results_folder
        self.fba_1 = CobraFBA(model_1, knockout=knockout)
        self.fba_2 = CobraFBA(model_2, knockout=knockout)

        # Get all reactions and reaction bounds to avoid individuals being different lengths
        # Different media has different essential reactions
        self.fba_1.reaction_names, self.fba_1.initial_reaction_bounds = self.fba_1.find_reaction_intial_bounds(
            self.fba_1.non_essential_reactions + self.fba_1.essential_reactions)
        self.fba_2.reaction_names, self.fba_2.initial_reaction_bounds = self.fba_2.find_reaction_intial_bounds(
            self.fba_2.non_essential_reactions + self.fba_2.essential_reactions)
        self.create_evo(len(self.fba_1.reaction_names))
        self.toolbox.register("evaluate", self.jamie_fitness_function)

    def run_cycling_nsga2evo(self, num_run, gen_record, blood, famine, seed=None, population=None,
                             population_file=None, ngen=51):
        random.seed(seed)

        NGEN = ngen
        MU = 100
        self.population_file = population_file

        stats = tools.Statistics(lambda ind: ind.fitness.values)
        # stats.register("avg", numpy.mean, axis=0)
        # stats.register("std", numpy.std, axis=0)
        stats.register("min", np.min, axis=0)
        stats.register("max", np.max, axis=0)

        logbook = tools.Logbook()
        logbook.header = "gen", "evals", "std", "min", "avg", "max"

        if population == None:
            try:
                pop = self.toolbox.population_restart()
            except FileNotFoundError:
                pop = self.toolbox.population(n=MU)
        else:
            pop = population

        # Evaluate the individuals with an invalid fitness
        queue = multiprocessing.Queue()
        invalid_ind = [ind for ind in pop if not ind.fitness.valid]
        fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind, repeat(blood), repeat(famine), repeat(queue),
                                     repeat(queue))
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # This is just to assign the crowding distance to the individuals
        # no actual selection is done
        pop = self.toolbox.select(pop, len(pop))

        record = stats.compile(pop)
        logbook.record(gen=0, evals=len(invalid_ind), **record)
        # print(logbook.stream)

        # Begin the generational process
        for condition in range(1, ngen):
            # Vary the population

            # print('Cloning Population')
            offspring = tools.selTournamentDCD(pop, len(pop))
            offspring = [self.toolbox.clone(ind) for ind in offspring]
            # print('Cloning Completed')

            # print('Mutating Offspring')
            for ind1 in offspring:
                self.toolbox.mutate(ind1)
                del ind1.fitness.values
            # print('Mutation Completed')

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            # print('Evaluating Fitness')

            fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind, repeat(blood), repeat(famine),
                                         repeat(queue), repeat(queue))
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit
            # print("Fitness Evaluation Completed")

            # Select the next generation population
            # print('Selecting Population')
            pop = self.toolbox.select(pop + offspring, MU)
            # print('Selection Completed')
            record = stats.compile(pop)
            logbook.record(evals=len(invalid_ind), **record)

        # print('Evolution Completed')
        filename = str(num_run) + '_' + str(gen_record) + '.csv'
        filepath = os.path.join(self.results_folder, filename)
        with open(filepath, 'w') as csvfile:
            # print('Writing to ' + filename)
            writer = csv.writer(csvfile, delimiter=',')
            for p in pop:
                writer.writerow(p)
        # print("Final population hypervolume is %f" % hypervolume(pop, [11.0, 11.0]))

        return pop, logbook

'''
Evolver class that cycles between different media conditions for a given ratio of generations
'''
class CyclingFBAEvolver(CobraFBABase):
    def __init__(self, model_1, model_2, results_folder, knockout=False):
        super().__init__()
        self.population_file = 'false'
        self.results_folder = results_folder
        self.fba_1 = CobraFBA(model_1, knockout=knockout)
        self.fba_2 = CobraFBA(model_2, knockout=knockout)

        # Get all reactions and reaction bounds to avoid individuals being different lengths
        # Different media has different essential reactions
        self.fba_1.reaction_names, self.fba_1.initial_reaction_bounds = self.fba_1.find_reaction_intial_bounds(
            self.fba_1.non_essential_reactions + self.fba_1.essential_reactions)
        self.fba_2.reaction_names, self.fba_2.initial_reaction_bounds = self.fba_2.find_reaction_intial_bounds(
            self.fba_2.non_essential_reactions + self.fba_2.essential_reactions)
        self.create_evo(len(self.fba_1.reaction_names))
        self.toolbox.register("evaluate", self.cycling_fitness_function)

    def run_cycling_nsga2evo(self, num_run, gen_record, media_conditions, seed=None, population=None,
                             population_file=None, ngen=51):
        random.seed(seed)
        print(media_conditions)

        NGEN = ngen
        MU = 100
        self.population_file = population_file

        stats = tools.Statistics(lambda ind: ind.fitness.values)
        # stats.register("avg", numpy.mean, axis=0)
        # stats.register("std", numpy.std, axis=0)
        stats.register("min", np.min, axis=0)
        stats.register("max", np.max, axis=0)

        logbook = tools.Logbook()
        logbook.header = "gen", "evals", "std", "min", "avg", "max"

        if population == None:
            try:
                pop = self.toolbox.population_restart()
            except FileNotFoundError:
                pop = self.toolbox.population(n=MU)
        else:
            pop = population

        # Evaluate the individuals with an invalid fitness
        queue = multiprocessing.Queue()
        invalid_ind = [ind for ind in pop if not ind.fitness.valid]
        fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind, repeat(True), repeat(queue))
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # This is just to assign the crowding distance to the individuals
        # no actual selection is done
        pop = self.toolbox.select(pop, len(pop))

        record = stats.compile(pop)
        logbook.record(gen=0, evals=len(invalid_ind), **record)
        # print(logbook.stream)

        # Begin the generational process
        for condition in media_conditions:
            # Vary the population

            # print('Cloning Population')
            offspring = tools.selTournamentDCD(pop, len(pop))
            offspring = [self.toolbox.clone(ind) for ind in offspring]
            # print('Cloning Completed')

            # print('Mutating Offspring')
            for ind1 in offspring:
                self.toolbox.mutate(ind1)
                del ind1.fitness.values
            # print('Mutation Completed')

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            # print('Evaluating Fitness')
            if condition:
                fba = self.fba_1
            else:
                fba = self.fba_2

            fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind, fba, repeat(queue))
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit
            # print("Fitness Evaluation Completed")

            # Select the next generation population
            # print('Selecting Population')
            pop = self.toolbox.select(pop + offspring, MU)
            # print('Selection Completed')
            record = stats.compile(pop)
            logbook.record(evals=len(invalid_ind), **record)
            # print(logbook.stream)

        # print('Evolution Completed')
        filename = str(num_run) + '_' + str(gen_record) + '.csv'
        filepath = os.path.join(self.results_folder, filename)
        with open(filepath, 'w') as csvfile:
            # print('Writing to ' + filename)
            writer = csv.writer(csvfile, delimiter=',')
            for p in pop:
                writer.writerow(p)
        # print("Final population hypervolume is %f" % hypervolume(pop, [11.0, 11.0]))

        return pop, logbook

    def cycling_fitness_function(self, individual, fba, queue):

        fba.set_reaction_bounds(individual)
        p = multiprocessing.Process(target=fba.run_fba, args=(queue,))
        p.start()
        p.join(3)

        if p.is_alive():
            # print('This is taking too long...')
            p.terminate()
            p.join()
            growth = 0
        else:
            growth = queue.get()

        if growth > 0.1:
            reactions = sum(individual)
        else:
            reactions = len(self.fba_1.non_essential_reactions) + len(self.fba_1.essential_reactions)
        return growth, reactions


def set_media_conditions(high_media=0, low_media=0, n_gens=2, record_num=2):
    high = np.ones(high_media, dtype=bool)
    low = np.zeros(low_media, dtype=bool)

    # records = int(n_gens / record_num)
    total = high_media + low_media
    repeats = int(n_gens / total)
    ratio_array = np.concatenate([high, low])
    record = int(n_gens / record_num)

    x = np.tile(ratio_array, int(record_num / total))
    return np.tile(x, (int(n_gens / record_num), 1))


'''
Main FBA Evolver that sets up the cobrapy model
The algorithm is setarted using the run_nsga2evo method

Arguments:
input -- filepath to the json fba model to be evovled
results -- filepath for the results directory
knockout -- reaction name to be knocked out, default is false

'''
class CobraFBAEvolver(CobraFBABase):
    def __init__(self, input, results_folder, knockout=False):
        super().__init__()
        self.population_file = 'false'
        self.results_folder = results_folder
        self.fba = CobraFBA(input, knockout=knockout)
        self.create_evo(len(self.fba.non_essential_reactions))
        self.toolbox.register("evaluate", self.fitness_function)

    '''
    The main nsga2 EA method that is run to evolve the model
    
    Arguments:
    num_run -- Integer identify of the replicate number
    gen_record -- list of integers to save the population
    ngen -- number of generations to run the algorithm for, default 51
    seed -- seed number for initialising pseudo random number generation, to be set for testing
    population: passing a previously saved population of individual for continuing evolution
    
    
    '''
    def run_nsga2evo(self, num_run, gen_record, ngen=51,seed=None, population=None, ):
        random.seed(seed)

        NGEN = ngen
        MU = 100


        stats = tools.Statistics(lambda ind: ind.fitness.values)
        # stats.register("avg", numpy.mean, axis=0)
        # stats.register("std", numpy.std, axis=0)
        stats.register("min", np.min, axis=0)
        stats.register("max", np.max, axis=0)

        logbook = tools.Logbook()
        logbook.header = "gen", "evals", "std", "min", "avg", "max"

        if population == None:
            try:
                pop = self.toolbox.population_restart()
            except FileNotFoundError:
                pop = self.toolbox.population(n=MU)
        else:
            pop = population

        # Evaluate the individuals with an invalid fitness
        queue = multiprocessing.Queue()
        invalid_ind = [ind for ind in pop if not ind.fitness.valid]
        fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind, repeat(queue))
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # This is just to assign the crowding distance to the individuals
        # no actual selection is done
        pop = self.toolbox.select(pop, len(pop))

        record = stats.compile(pop)
        logbook.record(gen=0, evals=len(invalid_ind), **record)
        # print(logbook.stream)

        # Begin the generational process
        for gen in range(1, NGEN):
            # Vary the population

            # print('Cloning Population')
            offspring = tools.selTournamentDCD(pop, len(pop))
            offspring = [self.toolbox.clone(ind) for ind in offspring]
            # print('Cloning Completed')

            # print('Mutating Offspring')
            for ind1 in offspring:
                self.toolbox.mutate(ind1)
                del ind1.fitness.values
            # print('Mutation Completed')

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            # print('Evaluating Fitness')
            fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind, repeat(queue))
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit
            # print("Fitness Evaluation Completed")

            # Select the next generation population
            # print('Selecting Population')
            pop = self.toolbox.select(pop + offspring, MU)
            # print('Selection Completed')
            record = stats.compile(pop)
            logbook.record(gen=gen, evals=len(invalid_ind), **record)
            # print(logbook.stream)

        # print('Evolution Completed')
        filename = str(num_run) + '_' + str(gen_record) + '.csv'
        filepath = os.path.join(self.results_folder, filename)
        with open(filepath, 'w') as csvfile:
            # print('Writing to ' + filename)
            writer = csv.writer(csvfile, delimiter=',')
            for p in pop:
                writer.writerow(p)
        # print("Final population hypervolume is %f" % hypervolume(pop, [11.0, 11.0]))

        return pop, logbook