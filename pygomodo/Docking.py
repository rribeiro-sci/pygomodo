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

import warnings
warnings.filterwarnings('ignore')

##### Libraries
import os, re, subprocess
import pandas as pd
from IPython.display import clear_output, display
import ipywidgets, py3Dmol 

homedirectory=os.path.dirname(__file__)

if 'LD_LIBRARY_PATH'  in os.environ.keys(): pass
else: os.environ['LD_LIBRARY_PATH']='/opt/rDock/lib'
if 'RBT_ROOT'  in os.environ.keys(): pass
else: os.environ['RBT_ROOT']='/opt/rDock/'
if 'PATH'  in os.environ.keys(): pass
else: os.environ['PATH']='/opt/venv/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/rDock/bin:/opt/hh-suite/bin:/opt/hh-suite/scripts'

import sys
sys.path.append(homedirectory)

from Utils import utils, Interactions

import platform
if platform.system() == 'Linux':
    vina_path = os.path.join(homedirectory,'opt/vina1.1.2/bin/vina')
elif platform.system() == 'Darwin':
    vina_path = os.path.join(homedirectory,'opt/vina1.1.2_mac/bin/vina')
else:
    raise ValueError('Platform unknown! The pygomodo was tested in Linux and Darwin (Mac) platforms.')


smina_path = os.path.join(homedirectory,'opt/smina/smina.static')

class VINA:
    """Molecular docking with AutoDock Vina (https://vina.scripps.edu/)."""
    
    def __init__(self, **kwargs):
        """
        Create a VINA instance
        """  
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
        self._outpath=os.getcwd()
        self._results=None
        self._reference_files_dir = os.path.join(homedirectory,'autogrids/reference_files/') ###REVIEW##
        if 'receptor_file' in kwargs:
            self._receptor_file = kwargs.pop('receptor_file')
            self._receptor_name = os.path.splitext(self._receptor_file)[0]
        if 'ligand_files' in kwargs:
            self._ligands = kwargs.pop('ligand_files')
        if 'ncpus' in kwargs:
            self._ncpu = kwargs.pop('ncpus')
        
        # ##DEFINING JOBNAME
        # if 'jobname' in kwargs:
        #     self._jobname = kwargs.pop('jobname')
        # else: 
        #     import datetime
        #     self._jobname = 'pyGOMODO_VINA_'+datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
        # self._jobdir = os.path.join(os.getcwd(), self._jobname)
        
        # #create the job folder  if not exists
        # if not os.path.isdir(self._jobdir):
        #     subprocess.call('mkdir ' + self._jobdir, shell=True)
        
        self._interactions_table=None
        return
    
    def Grid(self, **kwargs):
        """
        Performs an automatic detection of the binding pocket center coordinates, 
        for specific proteins families (also works for heteromultimer).

        .. warning: This funtion works exclusively with GPCRs.
        """  

        from Bio.PDB.vectors import Vector

        if 'receptor' in kwargs:
            self._receptor_file = kwargs.pop('receptor')
            self._receptor_name = os.path.splitext(self._receptor_file)[0]
        
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
    
    def Run(self, **kwargs):
        """vina docking engine itself.
        
        :parameter num_modes: exhaustiveness of the global search (roughly proportional to time) (default 20)
        :parameter size_x: size in the X dimension (Angstroms)
        :parameter size_y: size in the Y dimension (Angstroms)
        :parameter size_z: size in the Z dimension (Angstroms)
        :parameter center_x: X coordinate of the center
        :parameter center_Y: Y coordinate of the center
        :parameter center_Z: Z coordinate of the center
        :parameter exhaustiveness: exhaustiveness of the global search (roughly proportional to time) (default 8)
        :parameter energy_range: exhaustiveness of the global search (roughly proportional to time) (default 4)
        :parameter cpus: he number of CPUs to use (default 1)
        
        .. note:: Size and Center parameters can be automatica aquired by running Vina.Grid() funtion. 
                    Those parameters will be passed automatically to the function.
        """
        import openbabel
        openbabel.obErrorLog.StopLogging()
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
            raise ValueError('Vina box undefined. Please introduce box coordinates or run vina.Grid() function for automatic grid detection.')            

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
            self._modefiles[tail.split('_')[-3]]=mode_files

            #vinabestresult = utils.vina_extract_best_result(dockfile, root+'_best.pdbqt')
            #vinabestresultpdb = utils.pdbqt2pdb(vinabestresult, self._outpath)
            receptorpdb = utils.pdbqt2pdb(receptor4dock, self._outpath)

        print('\nDone!')
        self._results=vina_dict
        return 

    def ViewPoses(self,surface=False, fancy=False):
        """3D visualization of the docking poses with py3Dmol."""
        #£import ipywidgets
        
        def inception(ligand):
            ligand=ligand
            def scores(mode):
                #def df(ligand):
                #£from IPython.core.display import display
                #£import pandas as pd
                    #data = self._results[self._receptor_name+'_'+ligand+'_vina_out']
                    #return display(pd.DataFrame.from_dict(data))
                
                
                df_scores = pd.DataFrame.from_dict(self._results[self._receptor_name+'_'+ligand+'_vina_out'])
                header = {'selector': 'th:not(.index_name)', 'props': [
                        ('background-color', 'white'), ('font-size', '13px'),
                        ('color', 'black'), ('border', '2px solid white')]}

                poses = {'selector': 'th.col_heading.level0', 'props': [
                        ('font-size', '13px'),('color', 'white'), 
                        ('background-color', '#5D884E'),
                        ("border", "2px solid white")]}

                row = {'selector': '.l0', 'props': 'color:blue;'}
                slice_ = df_scores.index[int(mode)-1]
                df_display = df_scores.style.set_table_styles([header, poses, row]).set_properties(**{'background-color': 'lightgreen', 'font-weight':'bold'}, subset=slice_).hide_index()
                display(df_display)
                
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
            
            options2=[x.split('_')[-1].split('.')[0] for x in self._modefiles[ligand]]                    
            ipywidgets.interact(scores, mode=ipywidgets.Dropdown(options=options2,value=options2[0], description='Mode:', disabled=False))
            return 
            
        options=[x for x in self._modefiles.keys()]

        ipywidgets.interact(inception, ligand=ipywidgets.Dropdown(options=options,value=options[0], description='Ligand:', disabled=False))
        return 

    def AnalyseInteractions(self):
        
        """Protein-Ligand interaction analysis with ProLIF (https://prolif.readthedocs.io/en/latest/index.html)."""

        from openbabel import pybel  

        self._interactions_table={}
        receptor_file = self._receptor_name+'_grid.pdb'
        receptorMol = next(pybel.readfile('pdb', receptor_file)).write('pdb') 
        
        for ligand in self._modefiles:
            for mode in self._modefiles[ligand]:
                ligand_file = mode.split('.')[0]+'.pdb'
                poseMol = next(pybel.readfile('pdb', ligand_file)).write('pdb')
                self._interactions_table[ligand+'_'+mode.split('_')[-1].split('.')[0]] = Interactions.Prolif(receptorMol, poseMol)
        return

    def ViewInteractionsMap(self, map3D=True, map2D=True,surface=False, fancy=False, **kwargs):
        """Protein-Ligand interaction visualization.
        
        :parameter map3D: (True) If true 3D iteraction map is displayed
        :parameter map2D: (True) If true 2D interaction map id displayed
        :parameter opacity: (default 0.65) opacity of protein cartoon (map3D must be True)
        :parameter surface: (default False) add protein surface (map3D must be True)
        :parameter surface_opacity: (default 0.50) opacity of protein surface (map3D and surface must be True)
        :parameter fancy: (default Flase) outline color black (map3D must be True)
        """
        #£import ipywidgets
        from openbabel import pybel
        
        if 'opacity' in kwargs: 
            opacity=kwargs.pop('opacity')
        else: opacity=0.65
        if 'surface_opacity' in kwargs:
            surface_opacity=kwargs.pop('surface_opacity')
        else: surface_opacity=0.50

        if not self._interactions_table:
            VINA.AnalyseInteractions(self)
        
        def inception(ligand):
            ligand=ligand
            def network_interations(mode):
                lmol = next(pybel.readfile('pdb',self._receptor_name+'_'+ligand+'_vina_out_ligand_'+mode+'.pdb')).write('pdb')
                rmol = next(pybel.readfile('pdb', self._receptor_name+'_grid.pdb')).write('pdb')
                
                mode=mode
                if map3D is True:
                    mol3D = Interactions.Vismol(rmol, lmol, self._interactions_table[ligand+'_'+mode].copy(), 
                                                opacity=opacity,fancy=fancy, surface=surface, surface_opacity=surface_opacity)
                if map2D is True:
                    mol2D = Interactions.Lignet(self._interactions_table[ligand+'_'+mode],lmol)
                    return mol2D.display(width=600)
                                
            options2=[x.split('_')[-1].split('.')[0] for x in self._modefiles[ligand]] 
            ipywidgets.interact(network_interations, mode=ipywidgets.Dropdown(options=options2,value=options2[0], description='Mode:', disabled=False))
            return
        
        #options=[x.split('_')[1] for x in self._results]
        options=[x for x in self._modefiles.keys()]

        ipywidgets.interact(inception, ligand=ipywidgets.Dropdown(options=options,value=options[0],description='Ligand:',disabled=False))
        
        return

class RDOCK:
    """Molecular docking with rDock (http://rdock.sourceforge.net/)."""
    
    def __init__(self, **kwargs):
        """
        Create a rDock instance
        """ 
        self._outpath=os.getcwd()
        if 'title' in kwargs: self._TITLE = kwargs.pop('title')
        else: self._TITLE = 'TITLE' #resolve

        if 'receptor_file' in kwargs: 
            self._RECEPTOR_FILE = utils.MolConvert(kwargs.pop('receptor_file'), self._outpath, 'mol2', hyd=True, receptor=True)
        else: raise ValueError('Receptor file unknown.')

        if 'receptor_flex' in kwargs: self._RECEPTOR_FLEX = kwargs.pop('receptor_flex')
        else: self._RECEPTOR_FLEX = '3.0'
        self._interactions_table = None
        return

    def Cavity(self, **kwargs):
        """
        rbcavity – Cavity mapping and preparation of docking site 

        :parameter recerence mol: Required (kwargs str): reference molecule path
        :parameter REF_MOL:       Optional
        :parameter SITE_MAPPER:   Optional (str): default 'RbtLigandSiteMapper'
        :parameter RADIUS:        Optional (int): default '6.0' 
        :parameter SMALL_SPHERE:  Optional (int): default '1.0'
        :parameter MIN_VOLUME:    Optional (int): default '100'
        :parameter MAX_CAVITIES:  Optional (int): default '1'
        :parameter VOL_INCR:      Optional (int): default '0.0'
        :parameter GRIDSTEP:      Optional (int): default '0.5'
        :parameter SCORING_FUNCTION: Optional (str): default 'RbtCavityGridSF'
        :parameter WEIGHT:        Optional (int): default '1.0'
        
        .. note:: For a detailed description of the parameters please see rDock documentation http://rdock.sourceforge.net/.
        
        """
        if 'reference_mol' in kwargs: self._REF_MOL = utils.MolConvert(kwargs.pop('reference_mol'), self._outpath, 'sdf', hyd=True, receptor=True)

        if 'SITE_MAPPER' in kwargs: self._SITE_MAPPER = kwargs.pop('SITE_MAPPER')
        else: self._SITE_MAPPER = 'RbtLigandSiteMapper'
        if 'REF_MOL' in kwargs: self._REF_MOL = utils.MolConvert(kwargs.pop('REF_MOL'), self._outpath, 'sdf', hyd=True, receptor=True)

        if self._REF_MOL:pass
        else: raise ValueError('reference_mol unknown')

        if 'RADIUS' in kwargs: self._RADIUS = kwargs.pop('RADIUS')
        else: self._RADIUS = '6.0'
        if 'SMALL_SPHERE' in kwargs: self._SMALL_SPHERE = kwargs.pop('SMALL_SPHERE')
        else: self._SMALL_SPHERE = '1.0'
        if 'MIN_VOLUME' in kwargs: self._MIN_VOLUME = kwargs.pop('MIN_VOLUME')
        else: self._MIN_VOLUME = '100'
        if 'MAX_CAVITIES' in kwargs: self._MAX_CAVITIES = kwargs.pop('MAX_CAVITIES')
        else: self._MAX_CAVITIES = '1'
        if 'VOL_INCR' in kwargs: self._VOL_INCR = kwargs.pop('VOL_INCR')
        else: self._VOL_INCR = '0.0'
        if 'GRIDSTEP' in kwargs: self._GRIDSTEP = kwargs.pop('GRIDSTEP')
        else: self._GRIDSTEP = '0.5'
        if 'SCORING_FUNCTION' in kwargs: self._SCORING_FUNCTION = kwargs.pop('SCORING_FUNCTION')
        else: self._SCORING_FUNCTION = 'RbtCavityGridSF'
        if 'WEIGHT' in kwargs: self._WEIGHT = kwargs.pop('WEIGHT')
        else: self._WEIGHT = '1.0'
        

        self.__prm = """RBT_PARAMETER_FILE_V1.00
TITLE {}

RECEPTOR_FILE {}
RECEPTOR_FLEX {}

##################################################################
### CAVITY DEFINITION: REFERENCE LIGAND METHOD
##################################################################
SECTION MAPPER
SITE_MAPPER {}
REF_MOL {}
RADIUS {}
SMALL_SPHERE {}
MIN_VOLUME {}
MAX_CAVITIES {}
VOL_INCR {}
GRIDSTEP {}
END_SECTION

#################################
#CAVITY RESTRAINT PENALTY
#################################
SECTION CAVITY
SCORING_FUNCTION {}
WEIGHT {}
END_SECTION
            
        """.format(self._TITLE, self._RECEPTOR_FILE, self._RECEPTOR_FLEX, self._SITE_MAPPER,
        self._REF_MOL, self._RADIUS, self._SMALL_SPHERE, self._MIN_VOLUME, self._MAX_CAVITIES, 
        self._VOL_INCR, self._GRIDSTEP, self._SCORING_FUNCTION, self._WEIGHT)
        
        #print(prm)
        with open('docking.prm', 'w') as f:
            f.write(self.__prm)


        rbcavity = '/opt/rDock/bin/rbcavity -was -d -r docking.prm > rbcavity.log'
        subprocess.call(rbcavity, shell=True)
        
        #show output
        with open('rbcavity.log', 'r') as flog:
            for line in flog.readlines()[-5:]:
                print(line)
        return
    
    def Run(self, tethered=False,  smina_minimize=False, smina_score_only=False, **kwargs):
        """rbdock – the rDock docking engine itself.
        
        :parameter target_mol: Required (list or str): path of target molecule
        :parameter nruns: Optional (int): number of docking poses (default 10)
        :parameter output_name: Optional (str): name of the docking files (default 'docking_poses')
        :parameter tethered:    (False) If true it will run the tethered docking
        :parameter smina_minimize: (False) If true the docking solutions will be minimized with smina (https://sourceforge.net/projects/smina/)
        :parameter smina_score_only: (False) If true the docking solutions will be rescored with smina (https://sourceforge.net/projects/smina/)
        :parameters tethered_parameters: {'TRANS_MODE':'TETHERED',
                                            'ROT_MODE:'TETHERED',
                                            'DIHEDRAL_MODE='FREE ',
                                            'MAX_TRANS':1.0
                                            'MAX_ROT':30.0}
        
        .. note:: For a detailed description of the tethered parameters please see rDock documentation http://rdock.sourceforge.net/.
        """
        from openbabel import pybel
        from rdkit import Chem
        from rdkit.Chem import rdMolAlign

        if 'output_name' in kwargs: 
            output_name = kwargs.pop('output_name')
        else: output_name = 'docking_poses'
        
        if 'target_mol' in kwargs: 
            target_mol = kwargs.pop('target_mol')
            if isinstance(target_mol, list): pass
            else: target_mol=[target_mol]
        else: raise ValueError('molecules to dock unknown.') 

        if 'nruns' in kwargs: nruns = kwargs.pop('nruns')
        else: nruns = 10

        #TETHERED
        TRANS_MODE='TETHERED'
        ROT_MODE='TETHERED'
        DIHEDRAL_MODE='FREE'
        MAX_TRANS=1.0
        MAX_ROT=30.0 
        if 'tethered_parameters' in kwargs:
            tethered_parameters = kwargs.pop('tethered_parameters')
            if isinstance(tethered_parameters, dict):
                if 'TRANS_MODE' in tethered_parameters.keys(): TRANS_MODE=tethered_parameters.pop('TRANS_MODE')
                if 'ROT_MODE' in tethered_parameters.keys(): ROT_MODE=tethered_parameters.pop('ROT_MODE')
                if 'DIHEDRAL_MODE' in tethered_parameters.keys(): DIHEDRAL_MODE=tethered_parameters.pop('DIHEDRAL_MODE')
                if 'MAX_TRANS' in tethered_parameters.keys(): MAX_TRANS=tethered_parameters.pop('MAX_TRANS')
                if 'MAX_ROT' in tethered_parameters.keys(): MAX_ROT=tethered_parameters.pop('MAX_ROT')
            else: raise AttributeError('"tethered_parameters is not a dict"')
    
        if tethered is True:
            params='''TETHERED PARAMETERS:

TRANS_MODE {}
ROT_MODE {}
DIHEDRAL_MODE {}
MAX_TRANS {}
MAX_ROT {}'''.format(TRANS_MODE,ROT_MODE, DIHEDRAL_MODE, MAX_TRANS,MAX_ROT)
            print(params)
            tethered_prm_text=params="""
################################
# TETHERED SCAFFOLF
################################
SECTION LIGAND
TRANS_MODE {}
ROT_MODE {}
DIHEDRAL_MODE {}
MAX_TRANS {}
MAX_ROT {}
END_SECTION
""".format(TRANS_MODE,ROT_MODE, DIHEDRAL_MODE, MAX_TRANS,MAX_ROT)
            
            with open('docking.prm', 'w') as prm_tethered_file:
                prm_tethered_file.write(self.__prm)
                prm_tethered_file.write(tethered_prm_text)


        self._dockposes = {}
        self._scores = {}
        print('Running...\n')
        for target in target_mol:

            #Convert ligand
            target=os.path.basename(utils.MolConvert(target, self._outpath, 'sdf', hyd=True, receptor=True))
            target_name, target_ext = os.path.splitext(target)
            self._output_name = target_name+'_'+output_name
            self._output_name_sorted = '{}_sorted.sd'.format(self._output_name)
          
            #RUN rbDock
            print('\tdocking', os.path.basename(target_name),'\n')

            rbdock_cmd = '/opt/rDock/bin/rbdock -i {} -o {} -r docking.prm -p dock.prm -n {}'.format(target, self._output_name, nruns)
            subprocess.call(rbdock_cmd, shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            

            ### ANALYSIS
            if smina_minimize is False:
                #sort values
                if 'by' in kwargs: by=kwargs.pop('by')
                else: by='SCORE'
                sdsort_cmd = '/opt/rDock/bin/sdsort -n -f{} {}.sd > {}'.format(by,self._output_name, self._output_name_sorted)
                subprocess.call(sdsort_cmd, shell=True)

                #Split de .sd file
                #create a dict with pose id and sd string
                mols = list(pybel.readfile(target_ext.strip('.'), self._output_name_sorted))
                self._dockposes[target_name] = dict(zip(range(1,len(mols)+1),list(mols[idx].write('sd') for idx in range(len(mols)))))
            
                #create a dict with scores
                self._scores[target_name] = dict(zip(range(1,len(mols)+1), [dict() for i in range(len(mols))]))
            
                
                for i in range(len(mols)):
                    self._scores[target_name][i+1]['RDOCK Score']=mols[i].data['SCORE']

                #RMSD conformations
                with open(self._output_name_sorted, 'r') as f:
                    contents = f.read()
                contents = contents.replace(re.findall('\slibRbt.so/\S*\s\S*', contents)[0],'')
                conformations = contents.split("$$$$")[:-1]
                conf0 = Chem.MolFromMolBlock(conformations[0].replace('\n  rDOCK(R)','\n\n  rDOCK(R)'))

                rmsd=[]
                rmsd.append(0.00)
                for conf in range(1,len(conformations)):
                    rmsd.append(rdMolAlign.CalcRMS(Chem.MolFromMolBlock(conformations[conf]), conf0))
                for i in range(len(mols)):
                    self._scores[target_name][i+1]['RMSD']=rmsd[i]

            if smina_minimize is True:
                smina_command = '{} --receptor {} --ligand {}.sd --out {}_minimized.sd --minimize'.format(smina_path, self._RECEPTOR_FILE,self._output_name, self._output_name )
                smina_log = subprocess.check_output(smina_command, shell=True)
                smina_log = smina_log.decode('utf-8').splitlines()
                #Sort
                sdsort_cmd = '/opt/rDock/bin/sdsort -n -f{} {}_minimized.sd > {}_minimized_sorted.sd'.format('minimizedAffinity',self._output_name, self._output_name)
                subprocess.call(sdsort_cmd, shell=True)

                #Split de .sd file
                #create a dict with pose id and sd string
                mols = list(pybel.readfile(target_ext.strip('.'), '{}_minimized_sorted.sd'.format(self._output_name)))
                self._dockposes[target_name] = dict(zip(range(1,len(mols)+1),list(mols[idx].write('sd') for idx in range(len(mols)))))
            
                #create a dict with scores
                self._scores[target_name] = dict(zip(range(1,len(mols)+1), [dict() for i in range(len(mols))]))
            
                smina_score_command = '{} --receptor {} --ligand {}_minimized_sorted.sd --score_only'.format(smina_path, self._RECEPTOR_FILE,self._output_name)
                smina_score_log = subprocess.check_output(smina_score_command, shell=True)
                smina_score_log = smina_score_log.decode('utf-8').splitlines()
                affinity_values=[]
                intramolecular_values = []
                for line in smina_score_log:
                    try:
                        affinity = re.findall('^Affinity:\s*(\S*)\s*', line)[0]
                        affinity_values.append(affinity)
                    except: pass
                    try: 
                        intramolecular = re.findall('^Intramolecular energy:\s*(\S*)', line)[0]
                        intramolecular_values.append(intramolecular)
                    except:pass
                for i in range(len(mols)):
                    self._scores[target_name][i+1]['Affinity (Kcal/mol)']=affinity_values[i]
                    self._scores[target_name][i+1]['Intramolecular Energy (Kcal/mol)']=intramolecular_values[i]
                    self._scores[target_name][i+1]['RMSD minimized']=mols[i].data['minimizedRMSD']

                #RMSD conformations
                with open('{}_minimized_sorted.sd'.format(self._output_name), 'r') as f:
                    contents = f.read()
                conformations = contents.split("$$$$")[:-1]
                conf0 = Chem.MolFromMolBlock(conformations[0])

                rmsd=[]
                rmsd.append(0.00)
                for conf in range(1,len(conformations)):
                    rmsd.append(Chem.rdMolAlign.CalcRMS(Chem.MolFromMolBlock(conformations[conf].replace('\n\n\n\n','\n\n\n')), conf0))
                for i in range(len(mols)):
                    self._scores[target_name][i+1]['RMSD']=rmsd[i]
                
            
            if smina_minimize is False and smina_score_only is True:
                smina_command = '{} --receptor {} --ligand {}_sorted.sd --score_only'.format(smina_path, self._RECEPTOR_FILE,self._output_name)
                smina_score_log = subprocess.check_output(smina_command, shell=True)
                smina_score_log = smina_score_log.decode('utf-8').splitlines()
                affinity_values=[]
                intramolecular_values = []
                for line in smina_score_log:
                    try:
                        affinity = re.findall('^Affinity:\s*(\S*)\s*', line)[0]
                        affinity_values.append(affinity)
                    except: pass
                    try: 
                        intramolecular = re.findall('^Intramolecular energy:\s*(\S*)', line)[0]
                        intramolecular_values.append(intramolecular)
                    except:pass
                
                for i in range(len(mols)):
                    self._scores[target_name][i+1]['Affinity (Kcal/mol)']=affinity_values[i]
                    self._scores[target_name][i+1]['Intramolecular Energy (Kcal/mol)']=intramolecular_values[i]
        
        print('\nDONE!')
        return

    def ViewPoses(self, ref_mol=False,surface=False, fancy=False, surface_opacity=0.65, **kwargs):
        """3D visualization of the docking poses with py3Dmol.
        
        :parameter ref_mol: Optional (str) Reference ligand
        """
        #£import py3Dmol
        from openbabel import pybel
        #£import ipywidgets
        if 'opacity' in kwargs:
            opacity=kwargs.pop('opacity')
        else: opacity=0.65
        def inception(molname):
            
            #£from IPython.core.display import display
            #£import pandas as pd
                    
        
            def vismol(pose):
                df_scores = pd.DataFrame.from_dict(self._scores[molname], orient='index')

                header = {'selector': 'th:not(.index_name)', 'props': [
                        ('background-color', 'white'), ('font-size', '13px'),
                        ('color', 'black'), ('border', '2px solid white')]}

                poses = {'selector': 'th.col_heading.level0', 'props': [
                        ('font-size', '13px'),('color', 'white'), 
                        ('background-color', '#5D884E'),
                        ("border", "2px solid white")]}

                row = {'selector': '.l0', 'props': 'color:blue;'}
                slice_ = df_scores.index[pose-1]
                df_display = df_scores.style.set_table_styles([header, poses, row]).set_properties(**{'background-color': 'lightgreen', 'font-weight':'bold'}, subset=slice_).hide_index()

                display(df_display)

                pose=pose
                molview=py3Dmol.view(height=500)
                if fancy is True: molview.setViewStyle({'style':'outline','color':'black','width':0.1})
                mol1=next(pybel.readfile(os.path.splitext(self._RECEPTOR_FILE)[1].strip('.'), self._RECEPTOR_FILE)).write('pdb')
                mol2=self._dockposes[molname][pose]
                molview.addModel(mol1, 'pdb')
                molview.setStyle({'cartoon': {'color':'white'}})
                
                if surface is True: molview.addSurface(py3Dmol.VDW,{'opacity':surface_opacity,'color':'white'})
                molview.addModel(mol2, 'sd')
                molview.setStyle({'model':1},{'stick':{'colorscheme':'cyanCarbon','radius':0.2}})

                if ref_mol is True:
                    molview.addModel(next(pybel.readfile(os.path.splitext(self._REF_MOL)[1].strip('.'),self._REF_MOL)).write('pdb'), 'pdb')
                    molview.setStyle({'model':2},{'stick':{'colorscheme':'whiteCarbon','radius':0.2, 'opacity':opacity}})
                molview.setBackgroundColor('0xeeeeee')
                molview.zoomTo({'model':1})
                molview.show()
                return

            molname=molname
            options2=[x for x in self._dockposes[molname].keys()]
            ipywidgets.interact(vismol, pose=ipywidgets.Dropdown(options=options2, value=options2[0], description='Pose:', disable=False))

            return
        

        options1=[x for x in self._scores.keys()]
        ipywidgets.interact(inception, molname=ipywidgets.Dropdown(options=options1, value=options1[0], description='Ligand:', disable=False))

        return

    def AnalyseInteractions(self):
        """Protein-Ligand interaction analysis with ProLIF (https://prolif.readthedocs.io/en/latest/index.html)."""
        from openbabel import pybel

        self._interactions_table = {}
        receptorMol = next(pybel.readfile(os.path.splitext(self._RECEPTOR_FILE)[1].strip('.'), self._RECEPTOR_FILE)).write('pdb')
        
        for lig in self._dockposes.keys():
            self._interactions_table[lig]={}
            for pose in self._dockposes[lig]:
                poseMol = pybel.readstring('sd', self._dockposes[lig][pose]).write('pdb')
                self._interactions_table[lig][pose] = Interactions.Prolif(receptorMol, poseMol)            
        return

    def ViewInteractionsMap(self,map3D=True, map2D=True,surface=False, fancy=False, **kwargs):
        """Protein-Ligand interaction visualization.
        
        :parameter map3D: (True) If true 3D iteraction map is displayed
        :parameter map2D: (True) If true 2D interaction map id displayed
        :parameter opacity: (default 0.65) opacity of protein cartoon (map3D must be True)
        :parameter surface: (default False) add protein surface (map3D must be True)
        :parameter surface_opacity: (default 0.50) opacity of protein surface (map3D and surface must be True)
        :parameter fancy: (default Flase) outline color black (map3D must be True)
        """
        #£import ipywidgets
        from openbabel import pybel

        if 'opacity' in kwargs:
            opacity=kwargs.pop('opacity')
        else: opacity=0.65
        if 'surface_opacity' in kwargs:
            surface_opacity=kwargs.pop('surface_opacity')
        else: surface_opacity=0.50
        
        if not self._interactions_table:
            RDOCK.AnalyseInteractions(self)

        if self._interactions_table:
            def inception(molname):
                def network_interactions(pose):
                    lmol = pybel.readstring('sd',self._dockposes[molname][pose]).write('pdb')
                    rmol = next(pybel.readfile(os.path.splitext(self._RECEPTOR_FILE)[1].strip('.'), self._RECEPTOR_FILE)).write('pdb')
                    if map3D is True:
                        mol3D = Interactions.Vismol(rmol, lmol, self._interactions_table[molname][pose].copy(), 
                                                    opacity=opacity, fancy=fancy, surface=surface, surface_opacity=surface_opacity)
                    if map2D is True:
                        mol2D = Interactions.Lignet(self._interactions_table[molname][pose], lmol)
                        return mol2D.display(width=600)
                options2=[x for x in self._dockposes[molname].keys()]
                ipywidgets.interact(network_interactions, pose=ipywidgets.Dropdown(options=options2, value=options2[0], description='Pose: ', disabled=False))
                return
            options=[x for x in self._scores.keys()]
            ipywidgets.interact(inception, molname=ipywidgets.Dropdown(options=options, value=options[0], description='Ligand: ', disabled=False))

        return
    
