# Grammatical Evolution Molecular Design (GEMolecularDesign)
The **Grammatical Evolution Molecular Design** (GEMolecularDesign) is an evolutionary framework to evolves a molecule used as inhibitor of the SARS-CoV2 main protease.

GEMolecularDesign uses a Grammatical Evolution to proposed a molecule candidate based on Grammar that includes chemical fragments in the SMILES Language. The proposed molecule is optimized by Gaussian 09 using a semi-empirical optimization. Is used the Autodock to calculate the energy between the proposed molecule optimized and the SARS-CoV2 main protease.


## Usage

To run **GEMolecularDesign** only is necessary to indicate the config file.

```
python3 covid-py <config file>
```

## Config file

It is a json file with the following parameters:
- MaxTime: max time used to generate a proposal molecule.
- Babel: babel program path.
- Optimizer: optimizer program path.
- Gaussian: gaussian program path.
- Autogrid: autogrid program path.
- Autodock: autodock program path.
- QID: queue id.
- User: queue user .
- Processors: processors used by each proposal molecule.
- QueueLevel: queue level.
- Memory: memory used by each proposal molecule.
- Options: optimization level used by Gaussian.
- Python2: python2 program path.
- Python3: python3 program path.
- grammarFile: grammar file with chemical fragments.
- ScriptPath: autodock scripts path.
- MGLToosPath: MGLTools path.
- Receptor: reseptor file.
- ScriptLigand: script to prepare the ligand. 
- ScriptGPF: script to prepare gpf.
- ScriptDPF: script to prepare dpf.
- ScriptSummarize: script used to summarize the results.
- ResultsFile: result file.
- ResultDir: results dir path.
- Incomplete: incomplete dir path.
- PopulationSize: population size used in Grammatical Evolution.
- Iterations: iterations used in Grammatical Evolution.
- ItemSize: item size used in Grammatical Evolution.
- AttributeMin: attributed to be minimize from the summarize result file.

## License
This software was developed for the contribution https://doi.org/10.1039/D1CP04159B. The data, usage,
and files needed for its execution are under the Open Source
license.
