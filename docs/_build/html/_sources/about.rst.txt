What is the pyGOMoDo?
#######################

The **pyGOMoDo** library is divided into two main modules: the modeling and docking modules. 
The modeling module follows the same framework pipeline as the original web-service: 

1. Protein sequence alignment and hidden markov model profile generation with hh-suite(`Steinegger et al. 2019 <https://doi.org/10.1186/s12859-019-3019-7>`_).
2. Homology modeling with the Modeller software(`Webb and Sali 2016 <https://doi.org/10.1002/cpbi.3>`_).
3. Model quality assessment with QMEANBrane(`Studer, Biasini, and Schwede 2014 <https://doi.org/10.1093/bioinformatics/btu457>`_) through the Swiss-Model web-service API(`Waterhouse et al. 2018 <https://doi.org/10.1093/nar/gky427>`_).

The HMM profiles are generated against the most updated version of the UniRef30(`Mirdita et al. 2017 <https://doi.org/10.1093/nar/gkw1081>`_) 
database, and the templates are obtained by aligning the generated HMM profile against a database of HMM profiles 
for the GPCRs whose structures are known. This latter database was built by us from scratch by gathering all the human GPCRs 
structures deposited on the Protein Data Bank (`rcsb.org <https://www.rcsb.org/>`_). All the structures were parsed, cleaned and 
uniformed, and their HMM profiles were built also against the UniRef30 database(`Mirdita et al. 2017 <https://doi.org/10.1093/nar/gkw1081>`_). For each pdb structure 
we have calculated its coverage - percentage of resolved sequence amino acids - and retrieved, when possible, the conformational 
state of the receptor. Thus, the choice for the best templates for modeling can be done according to the “resolved coverage” of 
the structure and its conformational state (when possible). 

Differently from the original GOMoDo web-service, the docking module of pyGOMoDo works independently from the modeling module. 
The user can perform molecular docking over a structure modeled with the modeling module, another software or web-server, or even 
over an experimental structure. 
The docking module implements the Autodock Vina(`Trott and Olson 2010 <https://doi.org/10.1002/jcc.21334>`_) and rDOCK(`Ruiz-Carmona et al. 2014 <https://doi.org/10.1371/journal.pcbi.1003571>`_) docking programs in 3 
different docking protocols: 

1. A typical Autodock Vina protocol but with automatic selection of the orthosteric binding site (all the autodock Vina setting values can be toggled by the user), 
2. Docking of molecules on the binding site of crystallographic ligand with rDock
3. Tethered docking of molecules on a crystallographic ligand also with rDOCK. The docking module also implements the ProLIF (Protein-Ligand Interaction Fingerprints)(`Bouysset and Fiorucci 2021 <https://doi.org/10.1186/s13321-021-00548-6>`_) tool that allows the analysis of all the ligand-protein interactions.

The pyGOMoDO library was specifically designed having in mind its usage through jupyter notebooks providing the sweet-spot between
performance, flexibility and visualization. Tables visualization are rendered by the `pandas <https://pandas.pydata.org/>`_ library and protein 
visualization by the `py3Dmol <https://pypi.org/project/py3Dmol/>`_ (`Rego and Koes 2015 <https://doi.org/10.1093/bioinformatics/btu829>`_) library.