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

from turtle import heading
import warnings

warnings.filterwarnings('ignore')

from asyncio.constants import SENDFILE_FALLBACK_READBUFFER_SIZE
import os, re, subprocess
from sre_constants import SRE_FLAG_IGNORECASE
from pydoc import doc
from tkinter import E
from xml.dom.xmlbuilder import DocumentLS
import pandas as pd
from IPython.display import clear_output
import ipywidgets
homedirectory=os.path.dirname(__file__)

from src.utils import utils

import platform
if platform.system() == 'Linux':
    vina_path = os.path.join(homedirectory,'opt/vina1.1.2/bin/vina')
elif platform.system() == 'Darwin':
    vina_path = os.path.join(homedirectory,'opt/vina1.1.2_mac/bin/vina')
else:
    raise TypeError('Platform unknown! The pygomodo was tested in Linux and Darwin (Mac) platforms.')

mypython='python'


class get:
    def seq_from_uniprot(uniprotID, **kwargs):
        """
        Get Fasta sequence from Uniprot.

        :parameter uniprot:     required (str): Uniprot ID

        Example: seq_from_uniprot(uniprotid='P2314')
        """     

        from bioservices import UniProt
        u = UniProt(verbose=False)
    
        return u.get_fasta(uniprotID)
    
class gomodo:

    def __init__(self, **kwargs):
        """
        Create a GOMODO instance:

        :parameter jobname:     optional (str): create a directory 
        :parameter ncpus:       optional (int): number of cpus (default: 1)

        """  
        #PREPARE ENVIRONMENT
        try:
            self._hhsuitePath=os.environ['HHSUITEPATH']
            self._hhsuitescripts=os.environ['HHSUITESCRIPTS']
        except:
            print('Warning: HH-suite unknown')
        
        self._HHMdatabase=os.path.join(homedirectory,'databases/pyGOMODO_db/pyGOMODO')
        self._TemplatesDatabase = os.path.join(homedirectory,'databases/costum_db/wo_isoforms/gpcr_db_wo_isoforms')
        self._processedPDB=os.path.join(homedirectory,'databases/ProcessedPdbs_wo_isoforms')
        self._GPCR_Ref = os.path.join(homedirectory,'databases/GPCR_Ref.sqlite3')

        ##DEFINING JOBNAME
        if 'jobname' in kwargs:
            self._jobname = kwargs.pop('jobname')
        else: 
            import datetime
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
        
    def createHHMprofile(self, **kwargs):
        """
        Run hhblits

        :parameter uniprot:     optional (str): Uniprot ID (required if sequence not given)
        :parameter filename:    optional (str): filename or path+filename (default: "sequence.seq")
        :parameter sequence:    optional (str): sequence of in FASTA format (required if UniprotID not given)
        :parameter ncpus:       optional (int): number of cpus (default: 1)
        :parameter rounds:      optional (int): number of HHBLITS runs (default: 1)
        :parameter databse:     optional (str): dir path of the uniref30 database

        Example:
        """       
        #DIFINING SHARED VARIABLES
        if 'uniprotID' in kwargs:
            self._uniprotID = kwargs.pop('uniprotID')
        if 'sequence' in kwargs:
            self._fasta = kwargs.pop('sequence')
        if self._uniprotID.strip():
            self._fasta =get.seq_from_uniprot(uniprotID=self._uniprotID)
            with open(os.path.join(self._jobdir,'sequence.seq'), 'w') as f:
                f.write(self._fasta)
        elif not self._uniprotID.strip() and self._fasta:
            with open(os.path.join(self._jobdir,'sequence.seq'), 'w') as f:
                f.write(self._fasta)
        else: raise TypeError('UniprotID or Sequence unknown.')

        if 'ncpus' in kwargs: self._ncpus=kwargs.pop('ncpus')
        if 'rounds' in kwargs:
            self._rounds=kwargs.pop('rounds')
        if 'database' in kwargs: 
            self._HHMdatabase=kwargs.pop('database')
        


        #is sequence in fasta format? Save in .fas file
        output_fas = os.path.join(self._jobdir, 'sequence.fas')
        self._fasta80 = utils.searchFasta(self._fasta, output_fas)

        #Run HHBLITS
        print('Running...\n')
        utils.hhblits(self._hhsuitePath, output_fas, os.path.join(self._jobdir, "query.a3m"), self._logfilePath, self._HHMdatabase, self._ncpus, self._rounds)
       

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
        self._HHMprofiles=df
        print('Done!\n')
        return df

    def searchTemplates(self, **kwargs):
        """
        Run hhblits
        :parameter hhmprofile: optional (str): user profile path+filename (.a3m file)
        :parameter ncpus:       optional (int): number of cpus (default: 1)
        :parameter rounds:      optional (int): number of HHBLITS runs (default: 1)
        :parameter database     optional (str): dir path of the HHblists database

        Example:
        """       
        if 'hhmprofile' in kwargs: 
            self._hhmprofile = kwargs.pop('hhmprofile')
        else: self._hhmprofile=None    
        
        if 'ncpus' in kwargs: self._ncpus=kwargs.pop('ncpus')
        if 'rounds' in kwargs:
            self._rounds=kwargs.pop('rounds')
       
        if 'database' in kwargs: 
            self._TemplatesDatabase=kwargs.pop('database')
        

        #Run HHBLITS
        if self._hhmprofile:
            print('Running...\n')
            utils.hhblits(self._hhsuitePath, self._hhmprofile, os.path.join(self._jobdir + "templates.hhr"), self._logfilePath, self._TemplatesDatabase, self._ncpus, self._rounds)
        else:
            print('Running...\n')
            utils.hhblits(self._hhsuitePath,os.path.join(self._jobdir, "query.a3m"), os.path.join(self._jobdir, "templates.hhr"), self._logfilePath, self._TemplatesDatabase, self._ncpus, self._rounds)

        df =  utils.getConformation(self._jobdir, self._GPCR_Ref)
        self._templates=df
        print('Done!\n')
        return df

    def makeModels(self, **kwargs):
        #DEFINING VARIABLES
        self._loop = False      
        if 'loop' in kwargs:
            loop = kwargs.pop('loop')
            if loop == 'Yes':
                self._loop=True
            else: pass
        else: self._loop=False
        if 'trim' in kwargs:
            self._trim=kwargs.pop('trim')
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

        #create the folder "templates" if not exists
        templates_dir = os.path.join(self._jobdir, 'templates')
        if not os.path.isdir(templates_dir):
            subprocess.call('mkdir ' + templates_dir, shell=True)

        pir_name = os.path.join(templates_dir,  "mytemplate.pir")
     
        #myseq = sharedFunctions.extract_seqname(os.path.join(self._jobdir , "templates.hhr"))
        # write target sequence in pir format
        utils.hhmakemodel(self._hhsuitescripts,self._logfilePath, os.path.join(self._jobdir, "templates.hhr"), self._jobdir, 'sequence', templates_dir)


        #templates choosen by the user
        if 'templates' in kwargs:
            self._templates = kwargs.pop('templates')
        else: raise TypeError('Please select the templates')
        
    
        #copy the pdb (if does not exists yet) in the "templates" folder
        for i in range(0,len(self._templates)):
            pdb, chain = utils.divide_code(self._templates[i])
            if not os.path.isfile(os.path.join(templates_dir, pdb + ".pdb")):
                if self._humanDB:
                    subprocess.call("cp "+os.path.join(self._processedPDB, pdb + "_proc.pdb ") + os.path.join(templates_dir, pdb + ".pdb"), shell=True)
                    #print("cp "+os.path.join(self._processedPDB, pdb + "_proc.pdb ") + os.path.join(templates_dir, pdb + ".pdb"))
                else:
                    utils.download_pdb(pdb, self._cwd + "templates")

        #modeling (MODELLER)
        if self._loop:
            #print('*'*100 + '\n loop')
            models = utils.modellerLoop('sequence', self._templates, pir_name, self._cwd, templates_dir, self._ncpus, self._nmodels, self._nloops)
        else:
            #print('*'*100 + '\n automodel')
            models = utils.modeller('sequence', self._templates, pir_name, self._jobdir, templates_dir, self._ncpus, self._nmodels)
        self._models = models
        return models
    
    def vizModels(self, **kwargs):
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
    
    def qmeanbrane(self, **kwargs):
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
        #############################################################
        # To upload a local file found at /path/to/my_structure.pdb
        # ('rb' is recommended to allow zip file upload)
        response = requests.post(url=qmean_url, data={"email":email ,"method":"qmeanbrane"},files={"structure": open('models.tar.gz', 'rb')})
        ##############################################################
        
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
        
    def gdrive(self, **kwargs):
        import subprocess
        if 'path' in kwargs:
            path=kwargs.pop('path')
            command='cp -r '+os.path.join(self._cwd,self._jobname)+' '+path
            subprocess.call(command, shell=True)
        else: raise ValueError('Destination path unknown.')
        return

class docking:
    class vina:
        def __init__(self, **kwargs):  
            self._receptor_family ='GPCR' ## Review
            self._receptor_file=None
            self._receptor_name = None
            self._boxcenter_x=None
            self._boxcenter_y=None
            self._boxcenter_z=None
            self._boxsize_x=None
            self._boxsize_y=None
            self._boxsize_z=None
            self._center_box_vector=None
            self._exhaustiveness=8
            self._energy_range=4
            self._num_modes=20
            self._ncpu=1
            self._outpath=os.path.split(homedirectory)[0]
            self._results=None
            self._reference_files_dir = os.path.join(homedirectory,'autogrids/reference_files/') ###REVIEW##
            if 'receptor' in kwargs:
                self._receptor_file = kwargs.pop('receptor')
                self._receptor_name = os.path.splitext(self._receptor_file)[0]
            if 'ligands' in kwargs:
                self._ligands = kwargs.pop('ligands')
            if 'cpus' in kwargs:
                self._ncpu = kwargs.pop('cpus')
            if 'outpath' in kwargs:
                self._outpath=kwargs.pop('outpath')
            return
     
        def box(self, **kwargs):
            """
            perform an automatic detection of the binding pocket center coordinates, 
               for specific proteins families (also works for heteromultimer)
            """  

            from Bio.PDB.vectors import Vector

            if 'receptor' in kwargs:
                self._receptor_file = kwargs.pop('receptor')
                self._receptor_name = os.path.splitext(self._receptor_file)[0]
            
            if 'outpath' in kwargs:
                self._outpath=kwargs.pop('outpath')
            if self._receptor_file: pass
            else: raise ValueError('Unknown receptor. Please select a receptor file.')
  
            print('Pre-calculating grid...\n')
            refe = os.path.join(self._reference_files_dir, 'GPCr_Ref.pdb') ###REVIEW###
            r_chain='A' ###REVIEW###
            if 'chain' in kwargs:
                m_chain = kwargs.pop('chain')
            else: 
                m_chain = 'A'
            aligned_receptor, rot_tmp, tran_tmp, rms_tmp, score_tmp = utils.align3D(refe,r_chain, self._receptor_file, m_chain, output_directory=self._outpath)
            aligned_receptor_tmp, rot, tran, rms, score = utils.align3D(refe,r_chain, os.path.basename(aligned_receptor), m_chain)

            #create a dictionary containing the coordinates of the reference grids (as vectors) which can be
            #modified during the execution in order to not to modify the original reference grids file,
            #and a dictionary that will contain the transformed coordinates of the grids.

            ref_grids_dic = {}
            with open(os.path.join(self._reference_files_dir, 'GPCr_Ref_grids.txt'), "r") as f:
                for line in f:
                    toks = line.split()
                    chain_id = toks[0]
                    coor = Vector(toks[2],toks[3],toks[4])
                    ref_grids_dic.update({chain_id : coor})

            targ_grids_dic = {}

            #set a threshold for the alignment score to judge the goodness of the calculated alignment
            if score < 80 and rms > 7:

                print('Error: no good structural alignment found ')
                print(score,rms)
            else:
                #read the reference coordinates file and transform the relative coordinates currently stored in the dictionary,
                #then write them to the output target's grids coordinates file.
                #NB: this is important , because all the transformations performed of the reference chains change the
                #position of the chains themselves, but not the relative grids coordinates!
                ref_grid_coor = ref_grids_dic[r_chain]
                targ_grid_coor = Vector.right_multiply(ref_grid_coor, rot) + tran
                
                #store the transformed coordinates in the target grids dictionary
                targ_grids_dic.update({m_chain : targ_grid_coor})

            self._center_box_vector=targ_grid_coor
            
            ### PART II
            import py3Dmol
            
            #defining box function
            def visbox(objeto, bxi, byi, bzi, bxf, byf, bzf):
                opacity=1.0
                color='red'
                linewidth=0.3

                tx1 = bxi-(bxf/2)
                tx2 = bxi+(bxf/2)
                ty1 = byi-(byf/2)
                ty2 = byi+(byf/2)
                tz1 = bzi-(bzf/2)
                tz2 = bzi+(bzf/2)

                r1 = (bxi, ty1, tz1)
                r2 = (bxi, ty1, tz2)
                r3 = (bxi, ty2, tz1)
                r4 = (bxi, ty2, tz2)
                r5 = (tx1, byi, tz1)
                r6 = (tx1, byi, tz2)
                r7 = (tx2, byi, tz1)
                r8 = (tx2, byi, tz2)
                r9 = (tx1, ty1, bzi)
                r10 = (tx1, ty2, bzi)
                r11 = (tx2, ty1, bzi)
                r12 = (tx2, ty2, bzi)

                objeto.addBox({'center':{'x':r1[0],'y':r1[1],'z':r1[2]},'dimensions': {'w':bxf,'h':linewidth,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r2[0],'y':r2[1],'z':r2[2]},'dimensions': {'w':bxf,'h':linewidth,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r3[0],'y':r3[1],'z':r3[2]},'dimensions': {'w':bxf,'h':linewidth,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r4[0],'y':r4[1],'z':r4[2]},'dimensions': {'w':bxf,'h':linewidth,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r5[0],'y':r5[1],'z':r5[2]},'dimensions': {'w':linewidth,'h':byf,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r6[0],'y':r6[1],'z':r6[2]},'dimensions': {'w':linewidth,'h':byf,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r7[0],'y':r7[1],'z':r7[2]},'dimensions': {'w':linewidth,'h':byf,'d':linewidth},'opacity':opacity , 'color':color})
                objeto.addBox({'center':{'x':r8[0],'y':r8[1],'z':r8[2]},'dimensions': {'w':linewidth,'h':byf,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r9[0],'y':r9[1],'z':r9[2]},'dimensions': {'w':linewidth,'h':linewidth,'d':bzf},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r10[0],'y':r10[1],'z':r10[2]},'dimensions': {'w':linewidth,'h':linewidth,'d':bzf},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r11[0],'y':r11[1],'z':r11[2]},'dimensions': {'w':linewidth,'h':linewidth,'d':bzf},'opacity':opacity , 'color':color})
                objeto.addBox({'center':{'x':r12[0],'y':r12[1],'z':r12[2]},'dimensions': {'w':linewidth,'h':linewidth,'d':bzf},'opacity': opacity, 'color':color})

            #defining protein import and style
            def protvis(objeto, protein_name):
                mol1 = open(protein_name,'r').read()
                objeto.addModel(mol1, 'pdb')
                objeto.setStyle({'cartoon': {'color':'white'}})

            #difining visualization
            def vismol(**kwargs): 
                
                self._boxcenter_x=kwargs.pop('bxi')
                self._boxcenter_y=kwargs.pop('byi')
                self._boxcenter_z=kwargs.pop('bzi')
                self._boxsize_x=kwargs.pop('bxf')
                self._boxsize_y=kwargs.pop('byf')
                self._boxsize_z=kwargs.pop('bzf')

                mol_view = py3Dmol.view(width=1080, height=500,viewergrid=(1,2))
                visbox(mol_view,self._boxcenter_x,self._boxcenter_y,self._boxcenter_z,self._boxsize_x,self._boxsize_y, self._boxsize_z)
               
                protvis(mol_view,self._receptor_name+'_grid.pdb')
                mol_view.setBackgroundColor('0xeeeeee')
                mol_view.rotate(180, {'x':0,'y':1,'z':0},viewer=(0,1))
                mol_view.rotate(90, {'x':1,'y':0,'z':0},viewer=(0,0))
                mol_view.zoomTo()  
                mol_view.show()

            vector=self._center_box_vector
            clear_output(wait=True)
            ipywidgets.interact(vismol, 
                    bxi=ipywidgets.FloatSlider(description='Center X',value=vector[0],min=vector[0]-5,max=vector[0]+5, step=0.2) ,
                    byi=ipywidgets.FloatSlider(description='Center Y',value=vector[1],min=vector[1]-5,max=vector[1]+5, step=0.2),
                    bzi=ipywidgets.FloatSlider(description='Center Z',value=vector[2],min=vector[2]-5,max=vector[2]+5, step=0.2),
                    bxf=ipywidgets.FloatSlider(description='Size X',value=22.5, min=22.5,max=35, step=0.1),
                    byf=ipywidgets.FloatSlider(description='Size Y',value=22.5, min=22.5,max=35, step=0.1),
                    bzf=ipywidgets.FloatSlider(description='Size Z',value=22.5, min=22.5,max=35, step=0.1))
            
            return 
     
        def run(self, **kwargs):
            import subprocess  
      
            if self._center_box_vector:
                self._boxcenter_x = self._center_box_vector[0]
                self._boxcenter_y = self._center_box_vector[1]
                self._boxcenter_z = self._center_box_vector[2]
            if 'size_x' in kwargs:
                self._boxsize_x = kwargs.pop('size_x')
            if 'size_y' in kwargs:
                self._boxsize_y = kwargs.pop('size_y')
            if 'size_z' in kwargs:
                self._boxsize_z = kwargs.pop('size_z')
            if 'center_x' in kwargs:
                self._boxcenter_x = kwargs.pop('center_x')
            if 'center_y' in kwargs:
                self._boxcenter_y = kwargs.pop('center_y')
            if 'center_z' in kwargs:
                self._boxcenter_z = kwargs.pop('center_z')
            if 'exhaustiveness' in kwargs:
                self._exhaustiveness = kwargs.pop('exhaustiveness')
            if 'energy_range' in kwargs:
                self._energy_range = kwargs.pop('energy_range')
            if 'num_modes' in kwargs:
                self._num_modes = kwargs.pop('num_modes')

            #colocar sizes e center no self._box
            if not self._boxcenter_x or not self._boxcenter_z or not self._boxcenter_z:
                raise ValueError('Vina box undefined. Please introduce box coordinates or run grid() for automatic detection.')            

            print('Docking...\n')
            docking_files=[]
            #for receptor in self._receptors:
                
                #receptor_name, receptor_ext = os.path.splitext(receptor)
            if self._center_box_vector:
                receptor_aligned = self._receptor_name+'_grid.pdb'
                                    
            receptor4dock = utils.pdb2pdbqt(receptor_aligned, self._outpath, hyd=True, receptor=True)
        
            for ligand in self._ligands:
                ligand_name, ligand_ext = os.path.splitext(ligand)

                ligand_file = utils.pdb2pdbqt(ligand,self._outpath, hyd=True)

                vina_command = vina_path+' --receptor '+receptor4dock+' --ligand '+ligand_file+' --center_x '+str(self._boxcenter_x)+' --center_y '+str(self._boxcenter_y)+' --center_z '+str(self._boxcenter_z)+ ' --size_x '+str(self._boxsize_x)+ ' --size_y '+str(self._boxsize_y)+ ' --size_z '+str(self._boxsize_z)+ ' --exhaustiveness '+str(self._exhaustiveness)+ ' --energy_range '+str(self._energy_range)+ ' --num_modes '+str(self._num_modes)+ ' --cpu '+str(self._ncpu)+ ' --out '+os.path.join(self._outpath, self._receptor_name+'_'+ligand_name+'_vina_out.pdbqt')+' --log '+os.path.join(self._outpath,self._receptor_name+'_'+ligand_name+'_vina_out.log')
                
                print('\tDocking', ligand_name,'\n')
                subprocess.call(vina_command, shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                docking_files.append(os.path.join(self._outpath, self._receptor_name+'_'+ligand_name+'_vina_out.pdbqt'))
            
            print('\nAnalysing...')
            vina_dict={}
            self._modefiles={}
            for dockfile in docking_files:
                root, ext = os.path.splitext(dockfile)

                head, tail = os.path.split(root)

                with open(root+'.log', 'r') as vina_log:
                    lines = vina_log.read()
                    lines_splitted = lines.split('\n')
                vina_log_data = []
                for line in lines_splitted[24:-2]:
                    mode = line.strip().split()
                    vina_log_data.append(mode)
                vina_log_df = pd.DataFrame(vina_log_data, columns=['mode', 'affinity (kcal/mol)', 'rmsd l.b.','rmsd u.b.'])
                vina_log_df['mode'] = pd.to_numeric(vina_log_df['mode'])

                vina_log = vina_log_df[vina_log_df['mode'] <= 10]
                self._rui = vina_log
                #import ast
                #vina_dict[tail] = ast.literal_eval(self._resultsframe.to_json(orient='records'))
                vina_dict[tail] = vina_log.to_dict(orient='records')
                
                #Prepare vina files
                split_command = vina_path+'_split --input '+dockfile
                subprocess.call(split_command, shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                
                mode_files=[]
                for i in range(len(vina_log)):
                    if i <9:
                        mode_files.append(tail+'_ligand_0'+str(i+1)+'.pdbqt')
                        utils.pdbqt2pdb(tail+'_ligand_0'+str(i+1)+'.pdbqt', self._outpath)
                    else: 
                        mode_files.append(tail+'_ligand_'+str(i+1)+'.pdbqt')
                        utils.pdbqt2pdb(tail+'_ligand_'+str(i+1)+'.pdbqt', self._outpath)
                self._modefiles[tail.split('_')[1]]=mode_files

                vinabestresult = utils.vina_extract_best_result(dockfile, root+'_best.pdbqt')
                vinabestresultpdb = utils.pdbqt2pdb(vinabestresult, self._outpath)
                receptorpdb = utils.pdbqt2pdb(receptor4dock, self._outpath)
                #utils.mergePDBs(os.path.basename(receptorpdb), os.path.basename(vinabestresultpdb), os.path.splitext(os.path.basename(vinabestresultpdb))[0]+'_complex.pdb')

            print('\nDone!')
            self._results=vina_dict
            return 

        def results(self):
            import ipywidgets
            
            def inception(ligand):
                ligand=ligand
                def scores(mode):
                    def df(ligand):
                        from IPython.core.display import display
                        import pandas as pd
                        data = self._results[self._receptor_name+'_'+ligand+'_vina_out']
                        return display(pd.DataFrame.from_dict(data))
                    df(ligand)
                    
                    utils.viz(height=300, 
                        rmol=self._receptor_name+'_grid.pdb',
                        lmol=self._receptor_name+'_'+ligand+'_vina_out_ligand_'+mode+'.pdb',
                        bxi=self._boxcenter_x, 
                        byi=self._boxcenter_y, 
                        bzi=self._boxcenter_z, 
                        bxf=self._boxsize_x, 
                        byf=self._boxsize_y, 
                        bzf=self._boxsize_z)
                    return
                
                options2=[x.split('_')[5].split('.')[0] for x in self._modefiles[ligand]]                    
                ipywidgets.interact(scores, mode=ipywidgets.Dropdown(options=options2,value=options2[0], description='Mode:', disabled=False))
                return 
                
            options=[x.split('_')[1] for x in self._results]
            
            ipywidgets.interact(inception, ligand=ipywidgets.Dropdown(options=options,value=options[0], description='Ligand:', disabled=False))
            return 
        
        def scores(self):
            import pandas as pd
            import ipywidgets
            
            def df(ligand):
                data = self._results[self._receptor_name+'_'+ligand+'_vina_out']
                return pd.DataFrame.from_dict(data)
            options=[x.split('_')[1] for x in self._results]
            ipywidgets.interact(df, 
                    ligand=ipywidgets.Dropdown(
                                        options=options,
                                        value=options[0],
                                        description='Ligand:',
                                        disabled=False))
            return

        def viz(self):
            import py3Dmol, ipywidgets

            #defining box function
            def visbox(objeto, bxi, byi, bzi, bxf, byf, bzf):
                opacity=1.0
                color='red'
                linewidth=0.3

                tx1 = bxi-(bxf/2)
                tx2 = bxi+(bxf/2)
                ty1 = byi-(byf/2)
                ty2 = byi+(byf/2)
                tz1 = bzi-(bzf/2)
                tz2 = bzi+(bzf/2)

                r1 = (bxi, ty1, tz1)
                r2 = (bxi, ty1, tz2)
                r3 = (bxi, ty2, tz1)
                r4 = (bxi, ty2, tz2)
                r5 = (tx1, byi, tz1)
                r6 = (tx1, byi, tz2)
                r7 = (tx2, byi, tz1)
                r8 = (tx2, byi, tz2)
                r9 = (tx1, ty1, bzi)
                r10 = (tx1, ty2, bzi)
                r11 = (tx2, ty1, bzi)
                r12 = (tx2, ty2, bzi)

                objeto.addBox({'center':{'x':r1[0],'y':r1[1],'z':r1[2]},'dimensions': {'w':bxf,'h':linewidth,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r2[0],'y':r2[1],'z':r2[2]},'dimensions': {'w':bxf,'h':linewidth,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r3[0],'y':r3[1],'z':r3[2]},'dimensions': {'w':bxf,'h':linewidth,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r4[0],'y':r4[1],'z':r4[2]},'dimensions': {'w':bxf,'h':linewidth,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r5[0],'y':r5[1],'z':r5[2]},'dimensions': {'w':linewidth,'h':byf,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r6[0],'y':r6[1],'z':r6[2]},'dimensions': {'w':linewidth,'h':byf,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r7[0],'y':r7[1],'z':r7[2]},'dimensions': {'w':linewidth,'h':byf,'d':linewidth},'opacity':opacity , 'color':color})
                objeto.addBox({'center':{'x':r8[0],'y':r8[1],'z':r8[2]},'dimensions': {'w':linewidth,'h':byf,'d':linewidth},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r9[0],'y':r9[1],'z':r9[2]},'dimensions': {'w':linewidth,'h':linewidth,'d':bzf},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r10[0],'y':r10[1],'z':r10[2]},'dimensions': {'w':linewidth,'h':linewidth,'d':bzf},'opacity': opacity, 'color':color})
                objeto.addBox({'center':{'x':r11[0],'y':r11[1],'z':r11[2]},'dimensions': {'w':linewidth,'h':linewidth,'d':bzf},'opacity':opacity , 'color':color})
                objeto.addBox({'center':{'x':r12[0],'y':r12[1],'z':r12[2]},'dimensions': {'w':linewidth,'h':linewidth,'d':bzf},'opacity': opacity, 'color':color})

            #defining protein import and style
            def complxvis(objeto, protein_name, ligand_file):
                mol1 = open(protein_name,'r').read()
                mol2 = open(ligand_file,'r').read()
                objeto.addModel(mol1, 'pdb')
                objeto.setStyle({'cartoon': {'color':'white'}})
                objeto.addModel(mol2,'pdb')
                objeto.setStyle({'model':1},{'stick':{'colorscheme':'cyanCarbon','radius':0.2}})


            #difining visualization
            def vismol(**kwargs): 
                ligand=kwargs.pop('ligand')
                
                mol_view = py3Dmol.view(width=1080, height=500)
                visbox(mol_view,self._boxcenter_x,self._boxcenter_y,self._boxcenter_z,self._boxsize_x,self._boxsize_y, self._boxsize_z)
                ligand_file= self._receptor_name+'_'+ligand+'_vina_out_best.pdb'
                complxvis(mol_view,self._receptor_name+'_grid.pdb', ligand_file)
                mol_view.setBackgroundColor('0xeeeeee')
                mol_view.zoomTo({'model':1})  
                mol_view.show()

            vector=self._center_box_vector

            options=[x.split('_')[1] for x in self._results]
            ipywidgets.interact(vismol, ligand=ipywidgets.Dropdown(options=options,value=options[0],description='Ligand:',disabled=False))

        def analyse(self, all=True):
            self._interactions_table={}
            receptor_file = self._receptor_name+'_grid.pdb'
            
            if all:
                for ligand in self._modefiles:
                    for mode in self._modefiles[ligand]:
                        ligand_file = mode.split('.')[0]+'.pdb'
                        self._interactions_table[ligand+'_'+mode.split('_')[-1].split('.')[0]] = docking.interactions.prolif(receptor_file, ligand_file)
                self.__best=True
            else:
                for ligand in self._modefiles:
                    ligand_file = self._receptor_name+'_'+ligand+'_vina_out_ligand_01.pdb'
                    self._interactions_table[ligand] = docking.interactions.prolif(receptor_file, ligand_file)
                self.__best=False
            return

        def interactionsMap(self, **kwargs):
            import ipywidgets
            
            if 'opacity' in kwargs: 
                opacity=kwargs.pop('opacity')
            else: opacity=0.65

            if self.__best:
                #print('Rui Ribeiro')

                def inception(ligand):
                    ligand=ligand
                    def network_interations(mode):
                        mode=mode
                        mol3D = docking.interactions.vismol(self._receptor_name+'_grid.pdb',self._receptor_name+'_'+ligand+'_vina_out_ligand_'+mode+'.pdb',self._interactions_table[ligand+'_'+mode].copy(), opacity=opacity)
                        mol2D = docking.interactions.lignet(self._interactions_table[ligand+'_'+mode],self._receptor_name+'_'+ligand+'_vina_out_ligand_'+mode+'.pdb' )
                        return mol2D.display(width=600) #changeme
                    
                    options2=[x.split('_')[5].split('.')[0] for x in self._modefiles[ligand]] 
                    ipywidgets.interact(network_interations, mode=ipywidgets.Dropdown(options=options2,value=options2[0], description='Mode:', disabled=False))
                    return
                
                options=[x.split('_')[1] for x in self._results]
                ipywidgets.interact(inception, ligand=ipywidgets.Dropdown(options=options,value=options[0],description='Ligand:',disabled=False))

            
            
            else:
                def network_interations(**kwargs):
                    ligand=kwargs.pop('ligand')
                    mol3D = docking.interactions.vismol(self._receptor_name+'_grid.pdb',self._receptor_name+'_'+ligand+'_vina_out_ligand_01.pdb',self._interactions_table[ligand].copy(), opacity=opacity)
                    mol2D = docking.interactions.lignet(self._interactions_table[ligand],self._receptor_name+'_'+ligand+'_vina_out_ligand_01.pdb' )
                    return mol2D.display(width=600) #changeme
                
                options=[x.split('_')[1] for x in self._results]
                ipywidgets.interact(network_interations, ligand=ipywidgets.Dropdown(options=options,value=options[0],description='Ligand:',disabled=False))

    class interactions:

        def __init__(self, **kwargs):
         
            return

        def prolif(fi1, fi2):
            
            import prolif as plf
            from rdkit import Chem

            mol = Chem.MolFromPDBFile(fi1, removeHs=False)
            prot = plf.Molecule(mol)
            mol = Chem.MolFromPDBFile(fi2, removeHs=False)
            lig = plf.Molecule(mol)
            fp = plf.Fingerprint()
            fp.run_from_iterable([lig], prot, progress=False)
            df = fp.to_dataframe(return_atoms=True)
            return df

        def lignet(df, ligand):
            from prolif.plotting.network import LigNetwork
            from rdkit import Chem
            import prolif as plf
            mol = Chem.MolFromPDBFile(ligand, removeHs=False)
            lmol = plf.Molecule(mol)
            net = LigNetwork.from_ifp(df, lmol,
                          # replace with `kind="frame", frame=0` for the other depiction
                          kind="aggregate", threshold=.3,
                          rotation=270)
            return net

        def vismol(fi1, fi2, df, **kwargs):
            from rdkit import Chem, Geometry
            import py3Dmol, ipywidgets
            import prolif as plf
            def get_ring_centroid(mol, index):
                # find ring using the atom index
                Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
                ri = mol.GetRingInfo()
                for r in ri.AtomRings():
                    if index in r:
                        break
                else:
                    raise ValueError("No ring containing this atom index was found in the given molecule")
                # get centroid
                coords = mol.xyz[list(r)]
                ctd = plf.utils.get_centroid(coords)
                return Geometry.Point3D(*ctd)

            #ligand=kwargs.pop('ligand')
            #fi1 = self._receptor_name+'_grid.pdb'
            #fi2 = self._receptor_name+'_'+ligand+'_vina_out_best.pdb'
            mol = Chem.MolFromPDBFile(fi1, removeHs=False)
            pmol = plf.Molecule(mol)
            mol = Chem.MolFromPDBFile(fi2, removeHs=False)
            lmol = plf.Molecule(mol)
            #df = self._interactions_table[self._receptor_name+'_'+ligand+'_vina_out'].copy()
            colors = {
                    "HBAcceptor": "cyan",
                    "HBDonor": "cyan",
                    "Cationic": "green",
                    "PiStacking": "purple",
                    "Hydrophobic":"lime",
                    "Anionic":"gray",
                    "CationPi":"gray",
                    "EdgeToFace":"gray",
                    "FaceToFace":"gray",
                    "Interaction":"gray",
                    "MetalAcceptor":"gray",
                    "MetalDonor":"gray",
                    "PiCation":"gray",
                    "VdWContact":"gray",
                    "XBAcceptor":"gray",
                    "XBDonor":"gray",
                    "_BaseCationPi":"gray",
                    "_BaseHBond":"gray",
                    "_BaseIonic":"gray",
                    "_BaseMetallic":"gray",
                    "_BaseXBond":"gray",
                    "_Distance":"gray"
                }

            # JavaScript functions
            resid_hover = """function(atom,viewer) {{
                if(!atom.label) {{
                    atom.label = viewer.addLabel('{0}:'+atom.atom+atom.serial,
                        {{position: atom, backgroundColor: 'mintcream', fontColor:'black'}});
                }}
            }}"""
            hover_func = """
            function(atom,viewer) {
                if(!atom.label) {
                    atom.label = viewer.addLabel(atom.interaction,
                        {position: atom, backgroundColor: 'black', fontColor:'white'});
                }
            }"""
            unhover_func = """
            function(atom,viewer) {
                if(atom.label) {
                    viewer.removeLabel(atom.label);
                    delete atom.label;
                }
            }"""

            #v = py3Dmol.view(650, 600)
            v = py3Dmol.view(width=500, height=600) #changeme
            v.removeAllModels()

            models = {}
            mid = -1
            for i, row in df.T.iterrows():
                lresid, presid, interaction = i
                lindex, pindex = row[0]
                lres = lmol[0]
                pres = pmol[presid]
                # set model ids for reusing later
                for resid, res, style in [(lresid, lres, {"colorscheme": "cyanCarbon"}),
                                        (presid, pres, {})]:
                    if resid not in models.keys():
                        mid += 1
                        v.addModel(Chem.MolToMolBlock(res), "sdf")
                        model = v.getModel()
                        model.setStyle({}, {"stick": style})
                        # add residue label
                        model.setHoverable({}, True, resid_hover.format(resid), unhover_func)
                        models[resid] = mid
                # get coordinates for both points of the interaction
                if interaction in ["PiStacking", "EdgeToFace", "FaceToFace", "PiCation"]:
                    p1 = get_ring_centroid(lres, lindex)
                else:
                    p1 = lres.GetConformer().GetAtomPosition(lindex)
                if interaction in ["PiStacking", "EdgeToFace", "FaceToFace", "CationPi"]:
                    p2 = get_ring_centroid(pres, pindex)
                else:
                    p2 = pres.GetConformer().GetAtomPosition(pindex)
                # add interaction line
                v.addCylinder({"start": dict(x=p1.x, y=p1.y, z=p1.z),
                            "end":   dict(x=p2.x, y=p2.y, z=p2.z),
                            "color": colors[interaction],
                            "radius": .1,
                            "dashed": True,
                            "fromCap": 1,
                            "toCap": 1,
                            })
                # add label when hovering the middle of the dashed line by adding a dummy atom
                c = Geometry.Point3D(*plf.utils.get_centroid([p1, p2]))
                modelID = models[lresid]
                model = v.getModel(modelID)
                model.addAtoms([{"elem": 'Z',
                                "x": c.x, "y": c.y, "z": c.z,
                                "interaction": interaction}])
                model.setStyle({"interaction": interaction}, {"clicksphere": {"radius": .5}})
                model.setHoverable(
                    {"interaction": interaction}, True,
                    hover_func, unhover_func)

            # show protein
            if 'opacity' in kwargs: 
                opacity=kwargs.pop('opacity')
            else: opacity=0.65
            mol = Chem.RemoveAllHs(pmol)
            pdb = Chem.MolToPDBBlock(mol, flavor=0x20 | 0x10)
            v.addModel(pdb, "pdb")
            model = v.getModel()
            model.setStyle({}, {"cartoon": {"style":"edged", "opacity":opacity}})
            v.zoomTo({"model": list(models.values())})
            
            return v.show()

        