import numpy as np
import pandas as pd
import subprocess as sp
from time import time_ns
#import os
import stat as stt
from grammatical import Grammar
from os import path, stat, chmod, mkdir, chdir, getcwd
import shutil

parameters = {  "Babel" : "obabel",
                "Optimizer" : "obminimize",  
                "Gaussian" : "g09", 
                "Autogrid": "autogrid4",
                "Autodock": "/home/autodoc/adsuite/src100/autodock/self/autodock4",
                "QID" : "MSOTPR", 
                "Processors": "16",
                "QueueLevel": "large",
                "Memory": "64",
                "Options": "PM6 opt",
                "fileName" :"",
                "Python2" : "python2",
                "grammarFile" : "grammar_siesta.bnf",
                "ScriptPath" :"/home/mgltools/mglt_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24",
                "MGLToosPath" : "/home/mgltools/mglt_1.5.6/bin",
                "Receptor" : "6y2fsl.pdbqt",
                "ScriptLigand": "prepare_ligand4.py",
                "ScriptGPF": "prepare_gpf4.py",
                "ScriptDPF": "prepare_dpf4.py",
                "ScriptSummarize": "summarize_results4.py"}


def execute(commandLine):
    print(commandLine)
    try:
        p = sp.Popen(commandLine, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        p.wait()
    except:
        raise Exception(f'Can\'t execute the command {commandLine}')

def convert_Smiles_Gaussian(repr, fileName="output"):
    try:
        execute(f'{parameters["Babel"]} -:"{repr}" -o mol2 -O {fileName}.mol2 -h --gen3d')
        execute(f'{parameters["Optimizer"]} -ff UFF {fileName}.mol2 > {fileName}.tmp')
        execute(f'{parameters["Babel"]} -ipdb {fileName}.tmp -o com -O {fileName}.com')
    except Exception as err:
        raise Exception(f'Error converting from SMILE to Gaussian {err}')

def generate_GaussianFile(fileName):
    with open(f'{fileName}.inp', "w") as gfp:
        gfp.write(f'%nprocshared={parameters["Processors"]}\n')
        gfp.write(f'%mem={parameters["Memory"]}GB\n')
        gfp.write(f'%chk={fileName}.chk\n')
        gfp.write(f'#{parameters["Options"]}\n')
        gfp.write('\nCOVID19_GE: Generated Input\n\n')
        with open(f'{fileName}.com', "r") as cfg:
            gfp.write("".join(str(line) for line in cfg.readlines()[4:]))

def execute_Gaussian(fileName):
    try:
        execute(f'{parameters["Gaussian"]} < {fileName}.inp  >& {fileName}.log')
    except Exception as err:
        raise Exception(f'Error executing the Autogrind {err}')

def generate_ClusterJob(fileName):
    with open(f'{fileName}.job', "w") as fp:
        fp.write(f'#PBS -N {parameters["QID"]}\n')
        fp.write(f'#PBS -j oe\n')
        fp.write(f'#PBS -q {parameters["QueueLevel"]}\n')
        fp.write(f'#PBS -S /bin/bash\n')
        fp.write(f'#PBS -l nodes=1:ppn={parameters["Processors"]}\n')
        fp.write(f'#\n')
        fp.write(f'cd $PBS_O_WORKDIR\n')
        fp.write(f'time g09 < {fileName}.inp  >& {fileName}.log\n')

    st = stat(f'{fileName}.job')
    chmod(f'{fileName}.job', st.st_mode | stt.S_IEXEC)

def execute_ClusterJob(fileName):
    try:
        execute(f'qsub {fileName}.job')
    except Exception as err:
        raise Exception(f'Error executing the Cluster Job {err}')

def generate_Ligand(fileName):
    try:
        execute(f'formchk {fileName}.chk {fileName}.fchk')
        execute(f'{parameters["Babel"]} -ifchk {fileName}.fchk -o pdb -O {fileName}.pdb')
        execute(f'sh {path.join(parameters["MGLToosPath"], "pythonsh")} {path.join(parameters["ScriptPath"], parameters["ScriptLigand"])} -B amide -l {fileName}.pdb -o {fileName}.pdbqt')
    except Exception as err:
        raise Exception(f'Error generating the Ligand {err}')

def execute_Autogrid(fileName):
    try:
        execute(f'sh {path.join(parameters["MGLToosPath"], "pythonsh")} {path.join(parameters["ScriptPath"], parameters["ScriptGPF"])} -l{fileName}.pdbqt -r {parameters["Receptor"]} -o {fileName}.gpf -p npts="90,126,100" -p gridcenter="11.687,-1.471,20.252" -p spacing=0.141666666667')
        execute(f'{parameters["Autogrid"]} -p {fileName}.gpf -l {fileName}.glg')
    except Exception as err:
        raise Exception(f'Error executing the Autogrind {err}')

def execute_Autodock(fileName):
    try:
        execute(f'sh {path.join(parameters["MGLToosPath"], "pythonsh")} {path.join(parameters["ScriptPath"], parameters["ScriptDPF"])} -l{fileName}.pdbqt -r {parameters["Receptor"]} -o {fileName}.dpf')
        execute(f'{parameters["Autodock"]} -p {fileName}.dpf -l {fileName}.dlg')
    except Exception as err:
        raise Exception(f'Error executing the Autodock {err}')

def generate_Results(fileName):
    try:
        execute(f'sh {path.join(parameters["MGLToosPath"], "pythonsh")} {path.join(parameters["ScriptPath"], parameters["ScriptSummarize"])} -d./ -o {fileName}.cvs')
        resume = pd.read_csv(f'{fileName}.cvs', names=["#dlgfn","#incluster","#LE","#rmsd","#ats","#tors","#h_ats","#lig_eff"],skiprows=1)
        mejorconformacion=resume.iloc[0]
        eficiencia=mejorconformacion[7]
        return eficiencia
    except:
        return float("inf")

if __name__ == "__main__":
    old_path = getcwd()
    try:
        codons = np.random.randint(0, 10, 100)
        grammar = Grammar(file=parameters["grammarFile"], start_count=4, replace_chart="&")

        repr = grammar.mapping(codons).replace(" ","") 

        print(repr)

        parameters["fileName"] = str(time_ns())
        mkdir(f'{parameters["fileName"]}')
        shutil.copyfile(parameters["Receptor"], path.join(parameters["fileName"], parameters["Receptor"]))
        chdir(parameters["fileName"])

        convert_Smiles_Gaussian(repr, parameters["fileName"])
        generate_GaussianFile(parameters["fileName"])
        execute_Gaussian(parameters["fileName"])
        generate_ClusterJob(parameters["fileName"])
        #execute_ClusterJob(parameters["fileName"])
        while not path.exists(f'{parameters["fileName"]}.log'): pass

        generate_Ligand(parameters["fileName"])
        execute_Autogrid(parameters["fileName"])
        execute_Autodock(parameters["fileName"])
        generate_Results(parameters["fileName"])
    except Exception as err:
        print(err)
    finally:
        chdir(old_path)
