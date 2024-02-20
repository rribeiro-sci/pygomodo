#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2021 Rui Ribeiro
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
This python wrap-up of the GOMoDO webserver
"""
__author__ = "Rui Ribeiro"
__email__ = "rui.ribeiro@univr.it"

import os, subprocess, datetime, sys, warnings, shutil

warnings.filterwarnings('ignore')

homedirectory=os.path.dirname(__file__)
sys.path.append(homedirectory)


if 'HHSUITEPATH'  in os.environ.keys(): pass
else: os.environ['HHSUITEPATH']= os.path.join(os.path.dirname(os.path.dirname(shutil.which('hhblits'))),'bin')
if 'HHSUITESCRIPTS' in os.environ.keys(): pass
else: os.environ['HHSUITESCRIPTS']=os.path.join(homedirectory,'opt/scripts')

from Utils import utils, get

class Modeller:
    """Hommology Modeling with Modeller (https://salilab.org/modeller/)."""

    def __init__(self, **kwargs):
        """
        Create a Modeller instance:

        :parameter jobname:     optional (str): create a directory 
        :parameter ncpus:       optional (int): number of cpus (default: 1)

        """  
        #PREPARE ENVIRONMENT
        try:
            self._hhsuitePath=os.environ['HHSUITEPATH']
            self._hhsuitescripts=os.environ['HHSUITESCRIPTS']
        except:
            print('Warning: HH-suite unknown')
        
        self._HMMdatabase=os.path.join(homedirectory,'databases/GPCR_profiles_db/GPCR_profiles')
        self._TemplatesDatabase = os.path.join(homedirectory,'databases/GPCR_templates/hmm_profiles/GPCR_templates')
        self._processedPDB=os.path.join(homedirectory,'databases/GPCR_templates/pdb_clean_structures')
        #self._TemplatesDatabase = os.path.join(homedirectory,'databases/costum_db/wo_isoforms/gpcr_db_wo_isoforms')
        #self._processedPDB=os.path.join(homedirectory,'databases/ProcessedPdbs_wo_isoforms')
        self._GPCR_templates = os.path.join(homedirectory,'databases/GPCR_templates_dataframe.csv')
        
        ##DEFINING JOBNAME
        if 'jobname' in kwargs:
            self._jobname = kwargs.pop('jobname')
        else: 
            #£import datetime
            self._jobname = 'pyGOMODO_'+datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
        self._jobdir = os.path.join(os.getcwd(), self._jobname)
        
        #create the job folder  if not exists
        if not os.path.isdir(self._jobdir):
            subprocess.call('mkdir ' + self._jobdir, shell=True)
        
        #create log file is not exists
        self._filename='sequence.seq'
        self._name = self._filename[:self._filename.find('.')]
        
        self._logfilePath = os.path.join(self._jobdir, 'output.log')
        if not os.path.isfile(self._logfilePath):
            subprocess.call('touch '+ self._logfilePath, shell=True)

        #DEFINING SHARE VARIABLES
        if 'ncpus' in kwargs:
            self._ncpus = kwargs.pop('ncpus')
        else: self._ncpus=1    
        self._uniprotID=None
        self._fasta=None
        
        self._cwd = os.getcwd()+"/"
        self._humanDB = True
        self._rounds=1
        
    def CreateHMMProfile(self, **kwargs):
        """
        Creates HMM profiles with hhblits (https://github.com/soedinglab/hh-suite).

        :parameter uniprot:     optional (str): Uniprot ID (required if sequence not given)
        :parameter filename:    optional (str): filename or path+filename (default: "sequence.seq")
        :parameter sequence:    optional (str): sequence of in FASTA format (required if UniprotID not given)
        :parameter ncpus:       optional (int): number of cpus (default: 1)
        :parameter rounds:      optional (int): number of HHBLITS runs (default: 1)
        :parameter databse:     optional (str): dir path of the uniref30 database
        """       
        
        from IPython.display import clear_output, display

        #DIFINING SHARED VARIABLES
        if 'uniprotID' in kwargs:
            self._uniprotID = kwargs.pop('uniprotID')
        if 'sequence' in kwargs:
            self._fasta = kwargs.pop('sequence')
        if self._uniprotID.strip():
            try:
                self._fasta =get.seq_from_uniprot(uniprotID=self._uniprotID)
            except: pass
            with open(os.path.join(self._jobdir,'sequence.seq'), 'w') as f:
                f.write(self._fasta)
        elif not self._uniprotID.strip() and self._fasta:
            with open(os.path.join(self._jobdir,'sequence.seq'), 'w') as f:
                f.write(self._fasta)
        else: raise ValueError('UniprotID or Sequence unknown.')

        if 'ncpus' in kwargs: self._ncpus=kwargs.pop('ncpus')
        if 'rounds' in kwargs:
            self._rounds=kwargs.pop('rounds')
        if 'database' in kwargs: 
            self._HMMdatabase=kwargs.pop('database')
        


        #is sequence in fasta format? Save in .fas file
        output_fas = os.path.join(self._jobdir, 'sequence.fas')
        self._fasta80 = utils.searchFasta(self._fasta, output_fas)

        #Run HHBLITS
        print('Running...\n')
        utils.hhblits(self._hhsuitePath, output_fas, os.path.join(self._jobdir, "query.a3m"), self._logfilePath, self._HMMdatabase, self._ncpus, self._rounds)
    
        def parsehhr():
            import pandas as pd
            

            filename = os.path.join(self._jobdir,'sequence.hhr')
            col_names=utils.hhr_name_parser(filename)
            df = pd.DataFrame(columns=col_names)
            df = utils.hhr_summaryhit_parser(filename, df, col_names)
            pd.set_option("display.max_rows", None, "display.max_columns", None)
            #print(df.to_csv(r'results.txt', header=col_names, index=None, sep='	', mode='w'), '\nDone! See results.txt')
            
            return df
        df = parsehhr()    
        self._HMMprofiles=df
        clear_output(wait=True)

        header = {'selector': 'th:not(.index_name)', 'props': [
                ('background-color', 'white'), ('font-size', '13px'),
                ('color', 'black'), ('border', '2px solid white')]}

        poses = {'selector': 'th.col_heading.level0', 'props': [
                ('font-size', '13px'),('color', 'white'), 
                ('background-color', '#5D884E'),
                ("border", "2px solid white")]}

        row = {'selector': '.l0', 'props': 'color:blue;'}
        
        df_display = df.head(20).style.set_table_styles([header, poses, row]).hide_index()
        df_return = df.style.set_table_styles([header, poses, row]).hide_index()
        display(df_display)
        return df_return

    def SearchTemplates(self, **kwargs):
        """
        Searches for homologous templates with hhblits (https://github.com/soedinglab/hh-suite).
        :parameter hmmprofile: optional (str): user profile path+filename (.a3m file)
        :parameter ncpus:       optional (int): number of cpus (default: 1)
        :parameter rounds:      optional (int): number of HHBLITS runs (default: 1)
        :parameter database     optional (str): dir path of the HHblists database
        """   
        from IPython.display import clear_output, display
    
        if 'hmmprofile' in kwargs: 
            self._hmmprofile = kwargs.pop('hmmprofile')
        else: self._hmmprofile=None    
        
        if 'ncpus' in kwargs: self._ncpus=kwargs.pop('ncpus')
        if 'rounds' in kwargs:
            self._rounds=kwargs.pop('rounds')
    
        if 'database' in kwargs: 
            self._TemplatesDatabase=kwargs.pop('database')
    
        #Run HHBLITS
        if self._hmmprofile:
            print('Running...\n')
            utils.hhblits(self._hhsuitePath, self._hmmprofile, os.path.join(self._jobdir + "templates.hhr"), self._logfilePath, self._TemplatesDatabase, self._ncpus, self._rounds)
        else:
            print('Running...\n')
            utils.hhblits(self._hhsuitePath,os.path.join(self._jobdir, "query.a3m"), os.path.join(self._jobdir, "templates.hhr"), self._logfilePath, self._TemplatesDatabase, self._ncpus, self._rounds)

        df =  utils.getConformation(self._jobdir, self._GPCR_templates)
        self._templates=df.sort_values(by=['Sequence identity %', 'E-value', 'Coverage %'], ascending=[False, False, False]).reset_index().drop(columns=['index'])
        self._templates['Template']=self._templates.index
        self._templates=self._templates[['Template','No','Hit','Prob','E-value','P-value','Cols','Query HMM','Template HMM','Sequence identity %','Coverage %','State']]
        clear_output(wait=True)

        header = {'selector': 'th:not(.index_name)', 'props': [
                ('background-color', 'white'), ('font-size', '13px'),
                ('color', 'black'), ('border', '2px solid white')]}

        poses = {'selector': 'th.col_heading.level0', 'props': [
                ('font-size', '13px'),('color', 'white'), 
                ('background-color', '#5D884E'),
                ("border", "2px solid white")]}

        row = {'selector': '.l0', 'props': 'color:green;'}
        
        df_display = self._templates.head(20).style.set_table_styles([header, poses, row]).hide_index()
        display(df_display)

        return self._templates

    def MakeModels(self, **kwargs):
        """
        Runs Modeller
        
        :parameter templates: list of templates
        :parameter nmodels: (int) number of models
        :parameter loop: (boolean) loop refinement (default False)
        :parameter nloop: (int)  number of loop models computed (Default 2)
        :parameter ncpus: (int) number of cpus (default 1)
        
        """
        from IPython.display import display
 
        #DEFINING VARIABLES
        self._loop = False      
        if 'loop' in kwargs:
            loop = kwargs.pop('loop')
            if loop == True:
                self._loop=True
            else:pass
        if 'trim' in kwargs:
            trim=kwargs.pop('trim')
            if trim == True:
                self._trim=True
        else: self._trim = False
        if 'nloops' in kwargs:
            self._nloops=kwargs.pop('nloops')
        else:
            self._nloops=2
        if 'nmodels' in kwargs:
            self._nmodels = kwargs.pop('nmodels')
        else: self._nmodels = 100
        if 'ncpus' in kwargs:
            self._ncpus = kwargs.pop('ncpus')
    
        if 'humanDB' in kwargs:
            self._humanDB = kwargs.pop('humanDB')
        if 'hhsuitePath' in kwargs:
            self._hhsuitePath = kwargs.pop('hhsuitePath')
        if 'processedPDB' in kwargs:
            self._processedPDB = kwargs.pop('processedPDB')
        
        #templates choosen by the user
        if 'templates' in kwargs:
            self._selected_templates = kwargs.pop('templates')
            ntemplates=list(self._templates[self._templates.Hit == x].No.values[0] for x in self._selected_templates)
            selected_templates = list(x.replace(":", "_") for x in self._selected_templates)
        else: raise ValueError('Please select the templates')

        #create the folder "templates" if not exists
        templates_dir = os.path.join(self._jobdir, 'templates')
        if not os.path.isdir(templates_dir):
            subprocess.call('mkdir ' + templates_dir, shell=True)
        
        #copy the pdb (if does not exists yet) in the "templates" folder
        for i in range(0,len(self._selected_templates)):
            #pdb_hello, chain = utils.divide_code(self._selected_templates[i])
            
            pdb = selected_templates[i].split('_')[1]
            
            template = selected_templates[i]
            if not os.path.isfile(os.path.join(templates_dir, pdb)):
                if self._humanDB:
                    subprocess.call("cp "+os.path.join(self._processedPDB, template+'.pdb.gz') +' '+ os.path.join(templates_dir, template+'.pdb.gz'), shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                    #print("cp "+os.path.join(self._processedPDB, template+'.pdb') +' '+ os.path.join(templates_dir, template+'.pdb'))
                else:
                    utils.download_pdb(pdb, self._cwd + "templates")
        
        pir_name = os.path.join(templates_dir,  "mytemplate.pir")
    
        #myseq = sharedFunctions.extract_seqname(os.path.join(self._jobdir , "templates.hhr"))
        # write target sequence in pir format        
        utils.hhmakemodel(self._hhsuitescripts,self._logfilePath, os.path.join(self._jobdir, "templates.hhr"), self._jobdir, 'sequence', templates_dir, ntemplates)


        #modeling (MODELLER)
        if self._loop == True:
            #print('*'*100 + '\n loop')
            models = utils.modellerLoop('sequence', self._selected_templates, pir_name, self._cwd, templates_dir, self._ncpus, self._nmodels, self._nloops)
        else:
            #print('*'*100 + '\n automodel')
            models = utils.modeller('sequence', self._selected_templates, pir_name, self._jobdir, templates_dir, self._ncpus, self._nmodels)
        self._models = models


        header = {'selector': 'th:not(.index_name)', 'props': [
                ('background-color', 'white'), ('font-size', '13px'),
                ('color', 'black'), ('border', '2px solid white')]}

        poses = {'selector': 'th.col_heading.level0', 'props': [
                ('font-size', '13px'),('color', 'white'), 
                ('background-color', '#5D884E'),
                ("border", "2px solid white")]}

        row = {'selector': '.l0', 'props': 'color:blue;'}
        
        df_display = models.head(20).style.set_table_styles([header, poses, row]).hide_index()
        display(df_display)

        return models
    
    def ViewModels(self, **kwargs):
        """3D visualization of the models with py3Dmol."""
        import py3Dmol, ipywidgets
        
        
        #difining visualization
        def vismol(**kwargs):   
            model = kwargs.pop('model')          
            mol_view = py3Dmol.view()
            #complxvis(mol_view,'Model.'+model+'.pdb')
            mol_view.addModels(open('Model.'+model+'.pdb', 'r').read(),'pdb')
            mol_view.setStyle({'cartoon':{'arrows':True, 'tubes':False, 'style':'oval', 'color':'white'}})
            mol_view.setBackgroundColor('0xeeeeee')
            mol_view.zoomTo({'model':0})  
            mol_view.rotate(90, {'x':1,'y':0,'z':1},viewer=(0))
            mol_view.show()


        options=[x.split('.')[1] for x in list(self._models.name)]
        
        ipywidgets.interact(vismol, model=ipywidgets.Dropdown(options=options,value=options[0],description='Model:',disabled=False))  

        return
    
    def Qmeanbrane(self, **kwargs):
        """Submits a job to the Qmean web-service (https://swissmodel.expasy.org/qmean)
        
        :parameter models: required (list) list of models
        :parameter email: required (str) An email adress is required by the swissmodel webserver.
        
        .. note:: The e-mail provided is not stored, used or shared with third-parties excepting with swissmodel.org. For details regarding how Sisswmodel uses personal data please visit the Swissmodel website: https://swissmodel.expasy.org
        """
        
        import json
        import requests
        import time
        qmean_url = "https://swissmodel.expasy.org/qmean/submit/"

        if 'email' in kwargs:
            email = kwargs.pop('email')
        else: raise ValueError('Email unknown. An email adress is required by the swissmodel webserver.')

        if 'models' in kwargs:
            models = kwargs.pop('models')
            if isinstance(models, list):
                if len(models) > 5: raise ValueError('Too many models. (Max. 5 models)')
                else: 
                    import tarfile
                    tar = tarfile.open("models.tar.gz", "w:gz")
                    for name in models:
                        tar.add(name)
                    tar.close()
            else: raise ValueError('Please introduce a list of models.')
            
        else: raise ValueError('Please introduce a list of models.')
    
        def getstatus(response):
            current_status = requests.get(response.json()["results_json"])
            status = current_status.json()['status']
            return status
        
        def getresultspage(response):
            current_status = requests.get(response.json()["results_json"])
            page = current_status.json()['results_page']
            return page
        
        response = requests.post(url=qmean_url, data={"email":email ,"method":"qmeanbrane"},files={"structure": open('models.tar.gz', 'rb')})
        
        
        status = getstatus(response)
        if status == 'QUEUEING':
            print(status+"\n")
        while status == 'QUEUEING':
            #print(status+"\n")
            time.sleep(5)
            status=getstatus(response)
            
        if status == 'RUNNING':
            print(status+"\n")    
        while status == 'RUNNING':
            #print(status+"\n")
            time.sleep(5)
            status=getstatus(response)
        
        if status == 'COMPLETED':
            print(status+"\n")
            page = getresultspage(response)
            print(page)
        return
        
