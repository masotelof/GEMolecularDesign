import numpy as np
import pandas as pd
from IA import Item, binary_to_decimal
from itertools import count
from grammatical import Grammar
from chemistry.Components import Components
from os import mkdir, path, chmod, stat, remove, listdir, rename
from sys import exit
import subprocess as sp
import shutil
import tarfile
import glob
import stat as stt
from time import sleep
from datetime import datetime
import json

parameters = {
    "MaxTime": None,
    "Babel": None,
    "Optimizer": None,
    "Gaussian": None,
    "Autogrid": None,
    "Autodock": None,
    "QID": None,
    "User": None,
    "Processors": None,
    "QueueLevel": None,
    "Memory": None,
    "Options": None,
    "fileName": None,
    "Python2": None,
    "Python3": None,
    "grammarFile": None,
    "ScriptPath": None,
    "MGLToosPath": None,
    "Receptor": None,
    "ScriptLigand": None,
    "ScriptGPF": None,
    "ScriptDPF": None,
    "ScriptSummarize": None,
    "ResultsFile": None,
    "ResultDir": None,
    "Incomplete": None,
    "PopulationSize": None,
    "Iterations": None,
    "ItemSize": None,
    "AttributeMin": None
}

# load the basic configuration
with open('covid.config') as fp:
    parameters = json.load(fp)

if sum([1 if parameters[key]==None else 0 for key in parameters]) != 0:
    print("Hacen falta algunos parámetros en el archivo de configuración")
    exit(1)
    
if not path.isfile(parameters["Receptor"]):
    print(f'The receptor file {parameters["Receptor"]} doesn\'t exist')
    exit(1)

if path.isdir(parameters["Incomplete"]):
    print(f'Cleaning {parameters["Incomplete"]} directory')
    shutil.rmtree(parameters["Incomplete"])
mkdir(parameters["Incomplete"])

if path.isdir(parameters["ResultDir"]):
    print(f'Cleaning {parameters["ResultDir"]} directory')
    shutil.rmtree(parameters["ResultDir"])
mkdir(parameters["ResultDir"])

if path.isfile(parameters["ResultsFile"]):
    print(f'Backup the {parameters["ResultsFile"]} file')
    rename(parameters["ResultsFile"], f'{parameters["ResultsFile"]}.{datetime.today().strftime("%Y%m%d")}')

def conv_hex(value, maxitems):
    value_hex = hex(value)
    size = len(value_hex)
    return value_hex if size == maxitems else value_hex.replace("0x", f'0x{"0"*(maxitems-size)}')

def seconds(hour):
    try:
        values = hour.split(":")
        return (int(values[0]) * 60 + int(values[1])) * 60
    except:
        return 0

def execute(commandLine):
    try:
        p = sp.Popen(commandLine, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        p.wait()
    except:
        raise Exception(f'Can\'t execute the command {commandLine}')

def compress(id):
    try:
        with tarfile.open(f'{id}.tar.gz', 'w:gz') as fp:
            fp.add(f'{id}.log')
            for f in listdir(id):
                if path.isfile(path.join(id, f)):
                    fp.add(path.join(id, f))

    except Exception as err:
        print(err)

def delete_tmp(items):
    try:
        # copiar los temporales del mejor
        best = min(items)
        if not path.isfile(path.join(parameters["ResultDir"], f'{best.id}.tar.gz')):
            compress(best.id)
            shutil.move(f'{best.id}.tar.gz', path.join(parameters["ResultDir"], f'{best.id}.tar.gz'))

        for item in items:
            if path.isdir(item.id):
                print(f'Removing directory {item.id}')
                shutil.rmtree(item.id)
                
            for file in glob.glob(f'*{item.id}*'):
                if path.isfile(file):
                    print(f'Removing file {file}')
                    remove(file)
    except Exception as err:
        print(err)

def obtain_results(item):
    try:
        resume = pd.read_csv(f'{path.join(item.id, item.id)}.csv', names=["#dlgfn","#incluster","#LE","#rmsd","#ats","#tors","#h_ats","#lig_eff"],skiprows=1)
        mejorconformacion=resume.iloc[0]
        return mejorconformacion[parameters["AttributeMin"]]
    except:
        return float("inf")

def generate_jobfile(item):
    mkdir(item.id)

    shutil.copyfile(parameters["Receptor"], path.join(item.id, parameters["Receptor"]))
    execute(f'{parameters["Babel"]} -:"{item.phenotype}" -o mol2 -O {path.join(item.id, item.id)}.mol2 -h --gen3d')

    with open(f'{item.id}.job', "w") as fp:
        fp.write(f'#PBS -N {parameters["QID"]}{item.id}\n')
        fp.write(f'#PBS -j oe\n')
        fp.write(f'#PBS -q {parameters["QueueLevel"]}\n')
        fp.write(f'#PBS -S /bin/bash\n')
        fp.write(f'#PBS -l nodes=1:ppn={parameters["Processors"]}\n')
        fp.write(f'#\n')
        fp.write(f'cd $PBS_O_WORKDIR\n')
        fp.write(f'export LD_LIBRARY_PATH=/opt/glibc-2.14/lib:$LD_LIBRARY_PATH\n')
        fp.write(f'{parameters["Python3"]} execute_operations.py {item.id} > {item.id}.log\n')

    st = stat(f'{item.id}.job')
    chmod(f'{item.id}.job', st.st_mode | stt.S_IEXEC)
    execute(f'qsub {item.id}.job')
    print(f'executing {item.id}.job')

def evaluate_individuals(items):
    [generate_jobfile(item) for item in items]

    maxtime = seconds(parameters["MaxTime"])

    while True:
        sleep(20)
        p = sp.Popen(f'qstat -u {parameters["User"]} | grep {parameters["QID"]}', shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        p.wait()
        lines = p.stdout.readlines()
        total = len(lines)
        print(f'Executing {total} indivudal')
        if total==0:
            break

        # qdel $(qstat -u msotelo | grep msotelo | awk '{print $1}')
        # cat 0x002.log | grep -c Error
        for line in lines:
            values = [val for val in line.decode('utf8').split(" ") if val is not ""]
            itemtime = seconds(values[-1])
            if itemtime > maxtime:
                print(f'Killing {values[0]}')
                execute(f'qdel {values[0]}')
                line_id=values[3].replace(parameters["QID"], "")
                compress(line_id)
                shutil.move(f'{line_id}.tar.gz', path.join(parameters["Incomplete"], f'{line_id}.tar.gz'))

    for item in items:
        item.fitness = obtain_results(item)

if __name__ == "__main__":
    # parameters
    population_size = parameters["PopulationSize"]
    iterations = parameters["Iterations"]
    item_size = parameters["ItemSize"]

    max_items_size = len(hex((population_size * iterations)+1))
    
    grammar = Grammar(file=parameters["grammarFile"], start_count=4, replace_chart="&")
    component = Components(grammar)

    pop = []
    cont = 1
    for _ in range(population_size):
        item = Item(id=conv_hex(cont, max_items_size))
        item.initialize(item_size)
        item.evaluate(component)
        pop.append(item)
        cont += 1
    evaluate_individuals(pop)
    print(min(pop))
    delete_tmp(pop)

    with open(parameters["ResultsFile"], "w") as fp:
        fp.write(str(min(pop)))
        fp.write("\n")


    for i in range(iterations):
        sons = []
        # selection by binary tournament
        pos = np.arange(population_size)
        parents_pos1 = np.random.choice(pos, size=(population_size, 2))
        parents_pos2 = np.random.choice(pos, size=(population_size, 2))
        parents_pos1 = np.array([ind[0] if pop[ind[0]]<pop[ind[1]] else ind[1] for ind in parents_pos1])
        parents_pos2 = np.array([ind[0] if pop[ind[0]]<pop[ind[1]] else ind[1] for ind in parents_pos2])
        parents_pos = np.concatenate((parents_pos1, parents_pos2)).reshape((population_size, 2))

        # crossover by single point
        crossover_pos = np.random.randint(item_size-1, size = population_size)
        for i, ind, pos in zip(count(), parents_pos, crossover_pos):
            ind1, ind2 = ind
            sons.append((
                Item(values=np.concatenate((pop[ind1].values[:pos], pop[ind2].values[pos:])), id=conv_hex(cont, max_items_size)), 
                Item(values=np.concatenate((pop[ind2].values[:pos], pop[ind1].values[pos:])), id=conv_hex(cont+1, max_items_size))
            ))
            cont += 2

        mutation_pos = np.random.randint(item_size, size=(len(sons),2))
        
        sons1 = []
        sons2 = []
        for son, pos in zip(sons, mutation_pos):
            son1, son2 = son
            pos1, pos2 = pos
            son1.values[pos1] = 1 if son1.values[pos1] == 0 else 0
            son2.values[pos2] = 1 if son1.values[pos2] == 0 else 0
            son1.evaluate(component)
            son2.evaluate(component)
            sons1.append(son1)
            sons2.append(son2)
        
        evaluate_individuals(sons1)
        evaluate_individuals(sons2)
        new_pop = []
        for son1, son2 in zip(sons1, sons2):
            new_pop.append(son1 if son1<son2 else son2)

        with open(parameters["ResultsFile"], "a+") as fp:
            fp.write(str(min(new_pop)))
            fp.write("\n")

        pop = [min(pop)] + sorted(new_pop)[:population_size-1]
        print(min(pop))
        delete_tmp(pop)
        delete_tmp(sons1)
        delete_tmp(sons2)

        with open(parameters["ResultsFile"], "a+") as fp:
            fp.write(str(min(pop)))
            fp.write("\n")
