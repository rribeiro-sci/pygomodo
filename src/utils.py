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

import os, codecs, re, subprocess
import pandas as pd
from IPython.display import display
from ipywidgets import FileUpload
homedirectory=os.path.dirname(__file__)
mypython='python'

class utils:
    def pdbqt2pdb(fi, outpath):
        import openbabel
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats('pdbqt', 'pdb')
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol,fi)
        obConversion.WriteFile(mol, os.path.join(outpath, os.path.splitext(fi)[0]+'.pdb'))
        return os.path.join(outpath, os.path.splitext(fi)[0]+'.pdb')

    def pdb2pdbqt(fi, outpath, hyd=False, receptor=False):
        import openbabel
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats('pdb', 'pdbqt')
        if receptor==True: 
            obConversion.SetOptions('x', openbabel.OBConversion.OUTOPTIONS)
            obConversion.SetOptions('r', openbabel.OBConversion.OUTOPTIONS)
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol,fi)
        if hyd==True:
            mol.AddPolarHydrogens()
        else:
            mol.DeleteHydrogens()
        obConversion.WriteFile(mol, os.path.join(outpath, os.path.splitext(fi)[0]+'.pdbqt'))
        return os.path.join(outpath, os.path.splitext(fi)[0]+'.pdbqt')
    
    def vina_extract_best_result(pdbqt_result_filename, pdbqt_best_result_filename):

        #extract the first results from vina pdbqt result file
        with open(pdbqt_result_filename, "r") as file:
            lines = []
            for line in file:
                lines.append(line)
                columns=line.split()
                if columns[0] == 'ENDMDL': break

        with open(pdbqt_best_result_filename, "w") as out:
            out.writelines(lines)

        return pdbqt_best_result_filename

    def mergePDBs(fi1, fi2,fout):
        from Bio import PDB
        
        parser = PDB.PDBParser()
        sfi1 = parser.get_structure('receptor', fi1)
        sfi2 = parser.get_structure('ligand', fi2)
        
        sfi1_chains = list(sfi1.get_chains())
        sfi2_chains = list(sfi2.get_chains())
        
        sfi2_chains[0].id='B'
        sfi2_chains[0].detach_parent()
        
        sfi1[0].add(sfi2_chains[0])
        
        pdb_io = PDB.PDBIO()
        pdb_io.set_structure(sfi1)
        pdb_io.save(fout)

    def showMolecularStructures(mols):
        import rdkit
        from rdkit import Chem
        
        from rdkit.Chem import Draw
        from rdkit.Chem.Draw import IPythonConsole
        molecules=[]
        for m in mols:
            mol=Chem.MolFromPDBFile(m, removeHs=True)
            mol.RemoveAllConformers()
            molecules.append(mol)
               
        return Draw.MolsToGridImage(molecules, molsPerRow=4, subImgSize=(300, 300),legends=[x.split('.')[0] for x in mols])

    def viz(**kwargs):
        import py3Dmol

        if 'width' in kwargs: width = kwargs.pop('width')
        if 'height' in kwargs: height = kwargs.pop('height')
        if 'bxi' in kwargs: bxi = kwargs.pop('bxi')
        if 'byi' in kwargs: byi = kwargs.pop('byi')
        if 'bzi' in kwargs: bzi = kwargs.pop('bzi')
        if 'bxf' in kwargs: bxf = kwargs.pop('bxf')
        if 'byf' in kwargs: byf = kwargs.pop('byf')
        if 'bzf' in kwargs: bzf = kwargs.pop('bzf')
        if 'lmol' in kwargs: lmol = kwargs.pop('lmol')
        if 'rmol' in kwargs: rmol = kwargs.pop('rmol')


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
        #ligand=kwargs.pop('ligand')
        mol_view = py3Dmol.view(height=height)
        visbox(mol_view,bxi, byi, bzi, bxf, byf, bzf)
        #ligand_file= self._receptor_name+'_'+ligand+'_vina_out_best.pdb'
        complxvis(mol_view,rmol, lmol)
        mol_view.setBackgroundColor('0xeeeeee')
        mol_view.zoomTo({'model':1})  
        #mol_view.zoomTo()
        mol_view.show()


    def charPerLine(t,n):
        """adjust number of char per line
           this modifies a string t in a string seq with n char per line"""

        seq = ""
        for i in range(1,int((len(t)/n))+1):
            seq = seq + t[len(seq)-i+1:i*n] + "\n"
        if (len(t)%n!=0):
            seq = seq + t[len(seq)-int(len(t)/n)+1:len(t)] + "\n"
        return seq
    
    def searchFasta(sequence, output):
        """search in .seq file if there is a FASTA sequence
           check also for number of chars in each line (max 80)
           save it in a .fas file"""
        
        try:
            if re.match(">[^\n]*\n([GALMFWKQESPVICYHRNDT]{1,80}\n)+$", sequence) != None:
                with open(output, 'w') as f: 
                    f.write(sequence)
                return sequence
            else:
                if (sequence.find(">") == -1):
                    raise TypeError("""input value is not in FASTA""")
                else:
                    for idx in re.finditer(">", sequence):
                        pos = idx.span()[0]
                        seq = sequence[pos:]
                        #separe defline from sequence
                        pos = seq.find("\n")
                        defline = seq[:pos+1]
                        newseq = seq[pos+1:].upper()
                        newseq = newseq.replace("\t", "").replace(" ", "")
                        #adjust char per line
                        newseq = newseq.replace("\n", "")
                        newseq = utils.charPerLine(newseq,80)
                        final = defline + newseq
                        if re.match(">[^\n]*\n([GALMFWKQESPVICYHRNDT]{1,80}\n)+$", final) != None:
                            with open(output, 'w') as f: f.write(final)
                            return final
        except: raise TypeError("""input value is not in FASTA""")

    def hhblits(hhsuitePath,seqfile, output, logfilePath, hhblitsDbsPath, ncpus, rounds):
        pos = output.find(".")
        if output[pos:]==".a3m":
            ext = " -oa3m "
        else:
            ext = " -o "
        f = open(output, "w+")
        command = hhsuitePath+"hhblits -cpu " + str(ncpus) + " -i " + seqfile + " -d "+ hhblitsDbsPath  + ext + output  + " -n " + str(rounds)
        subprocess.call(command, shell=True)
        f.close()
        utils.update_logfile(logfilePath, command, "hhblits")
        return

    def update_logfile(logname, command, step):
        f = open(logname, "a+")
        f.write("\n" + step.upper() + ": " + command + "\n")
        f.close()
        return

    def getConformation(dir, database):
        filename = os.path.join(dir, 'templates.hhr')
        GPCR_Ref = database

        col_names=utils.hhr_name_parser(filename)
        df = pd.DataFrame(columns=col_names)
        df = utils.hhr_summaryhit_parser(filename, df, col_names)
        # print(df)
        df['Sequence identity'] = utils.extract_seq_identity(filename)
        df['State'] = utils.retrieve_state(GPCR_Ref, df)

        pd.set_option("display.max_rows", None, "display.max_columns", None)

        col_names.append('Sequence identity')
        col_names.append('State')
        #print(df.to_csv(r'results.txt', header=col_names, index=None, sep='	', mode='w'), '\nDone! See results.txt')
        
        return df

    def hhr_name_parser(filename):
        """extract the column names from the .hhr file 
           (summary of hits, thus without sequence identity and active/inactive state)"""
        f = open(filename)
        for i, line in enumerate(f):
            if i == 8:
                a = ' '.join(line.split())
                a = a.split(' ')
                a[8 : 10] = [' '.join(a[8 : 10])]
                a[9 : 11] = [' '.join(a[9 : 11])]
        f.close()
        return a

    def hhr_summaryhit_parser(filename, df, col_names):
        """extract for .hhr file the values of the summary hit list 
           store them in a dataframe (previously created) with column names in col_names"""
        f = open(filename)
        for i, line in enumerate(f):
            if i > 8:
                if line in ['\n','\n\r']:
                    break
                else:
                    a = ' '.join(line.split())
                    row = a.split(' ')
                    values_to_add = {}
                    if len(row[1])!=13:
                        #13 is the fixed len of the hit name from custom GPCR db
                        #if the results come from UniRef30/pdb70 
                        #the hit name is composed by more than one string after the split
                        #this consider only the first part (pdb code and chain)
                        row = row[0:2] + row [len(row)-9:len(row)]
                    for j in range(len(col_names)):
                        if len(row) == len(col_names):
                            values_to_add[col_names[j]] = row[j]
                        else: 
                            #cases in which in the .hhr file the last column (Template HMM)
                            #is separed by a blankspace. Thus in the list 'row' we will have 1 extra element.
                            #Join the 2 last elements of the list
                            row[9 : 11] = [' '.join(row[9 : 11])]
                            values_to_add[col_names[j]] = row[j]
                    row_to_add = pd.Series(values_to_add, name = row[0])
                    df = df.append(row_to_add)
        f.close()
        return df

    def extract_seq_identity(filename):
        """extract sequence identity between target and hit sequence from the .hhr file"""
        f = open(filename)
        line = f.readline()
        si = []
        while line:
            if line.startswith('Probab='):
                line = ' '.join(line.split())
                line = line.split(' ')
                value = line[4]
                pos = value.find('=')
                value = value[pos+1:]
                si.append(value)
            line = f.readline()
        return si
    
    def retrieve_state(db_name, df):
        """for each hit in df retrives the conformation state of the protein from 'GPCR_Ref.sqlite3' db 
           returns a list"""
        import sqlite3
        hit_names=df['Hit'].tolist()
        state_list = []

        conn = sqlite3.connect(db_name)
        cur = conn.cursor()
        for j in range(len(hit_names)):
            PDB_ID = hit_names[j][:4]
            cur.execute("select * from gpcr_pdb where pdb_id is '" + PDB_ID + "'")
            state = cur.fetchall()
            if len(state)==0:
                state_list.append('NaN')
            else:
                for i in state:
                    state = i[2]
                state_list.append(state)
        cur.close()

        return state_list

    def extract_seqname(filename):
        """extract the name of the target sequence from a .hhr file (result of hhblits over pdb70)"""
        f = open(filename, "r")
        line = f.readline().strip()
        pos = line.find(" ")
        line = line[pos+1:].strip()
        if (line.find(" ") != -1):
            pos_end = line.find(" ")
            line = line[:pos_end]
        seqname = line.replace("|", "_")
        return seqname.strip()

    def hhmakemodel(hhsuitescriptspath, logfilePath, input_file, mydir, seqname, templates_dir, nmodel = 1):
        """create a .ali file (pir format - necessary as Modeller input)"""
        command = mypython+" "+os.path.join(hhsuitescriptspath, "hhmakemodel.py ")+input_file+" "+templates_dir+" "+os.path.join(mydir,seqname + ".ali")+" ./ -m "+str(nmodel)
        
        subprocess.call(command, shell=True)
        utils.update_logfile(logfilePath, command, "hhmakemodel")

        return

    def divide_code(code):
        """obtain pdb and chain from the code extracted with extract_template_codes()"""
        if len(code)>4:
            pdb = code[:-2]
            chain = code[len(code)-1:].upper()
        else:
            pdb = code
            chain = ""
        return pdb, chain

    def extract_uknp_name(filename):
        """extract the name of the unknown protein from the .pir file"""
        f = open(filename, "r")
        line = f.readline()
        while line:
            if line.find(">") != -1:
                pos = line.find(";")
                name = line[pos+1:].strip()
                break
            else:
                line = f.readline()
        return name
    
    def download_pdb(code, templates_dir, ext = ".pdb"):
        """download pdb file (or cif - default ".pdb")"""
        command = "wget -P " + templates_dir + " http://files.rcsb.org/download/" + code + ext
        #print(command)
        subprocess.call(command, shell=True)
        return

    def modellerLoop(myseq, pdb_codes, pir_name, mydir, templates_dir, cpu, nmodels = 5, nloops = 2):
        """call Modeller to:
           1. align the target to the template sequence (using pir file format)
           2. create the models - 5 by default - with loop refinement (results in .pdb files)"""

        import modeller
        import modeller.parallel
        import modeller.automodel 
        os.chdir(mydir)

        j = modeller.parallel.job()
        for ncpu in range(cpu):
            j.append(modeller.parallel.local_slave())

        env = modeller.environ()
        env.io.atom_files_directory = [".", templates_dir]
        env.libs.topology.read(file = "$(LIB)/top_heav.lib")
        modeller.log.none()
        aln = modeller.alignment(env)
        for i in range(0,len(pdb_codes)):
            code, chain = utils.divide_code(pdb_codes[i])
            mdl = modeller.model(env, file = code, model_segment = ("FIRST:" + chain, "LAST:" + chain))
            aln.append_model(mdl, atom_files = code, align_codes = code + chain)

        for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                            ((1., 0.5, 1., 1., 1., 0.), False, True),
                                            ((1., 1., 1., 1., 1., 0.), True, False)):
            aln.salign(rms_cutoff = 3.5, normalize_pp_scores = False,
                    rr_file = "$(LIB)/as1.sim.mat", overhang = 30,
                    gap_penalties_1d = (-450, -50),
                    gap_penalties_3d = (0, 3), gap_gap_score = 0, gap_residue_score = 0,
                    dendrogram_file = myseq + "_templates.tree",
                    alignment_type = "tree",   # If "progresive", the tree is not
                                               # computed and all structues will be
                                               # aligned sequentially to the first
                    feature_weights = weights, # For a multiple sequence alignment only
                                               # the first feature needs to be non-zero
                    improve_alignment = True, fit = True, write_fit = write_fit,
                    write_whole_pdb = whole, output = "ALIGNMENT QUALITY")

        aln.write(file = myseq + "_temp.pir", alignment_format = "PIR")

        aln.salign(rms_cutoff = 1.0, normalize_pp_scores = False,
                rr_file = "$(LIB)/as1.sim.mat", overhang = 30,
                gap_penalties_1d = (-450, -50), gap_penalties_3d = (0, 3),
                gap_gap_score = 0, gap_residue_score = 0, dendrogram_file = "1is3A.tree",
                alignment_type = "progressive", feature_weights = [0]*6,
                improve_alignment = False, fit = False, write_fit = True,
                write_whole_pdb = False, output = "QUALITY")

        # # Read aligned structure(s):
        # aln = alignment(env)
        # aln.append(file="TvLDH_temp.pir", align_codes="all")
        aln_block = len(aln)

        # Read aligned sequence(s):
        # .ali file comes from hhmakemodel.py 
        uknp_name = utils.extract_uknp_name(myseq + ".ali")
        aln.append(file = myseq + ".ali", align_codes = uknp_name)

        # Structure sensitive variable gap penalty sequence-sequence alignment:
        aln.salign(output = "", max_gap_length = 20,
                gap_function = True,   # to use structure-dependent gap penalty
                alignment_type = "PAIRWISE", align_block = aln_block,
                feature_weights = (1., 0., 0., 0., 0., 0.), overhang = 0,
                gap_penalties_1d = (-450, 0),
                gap_penalties_2d = (0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
                similarity_flag = True)

        aln.write(file = pir_name, alignment_format = "PIR")

        codes = []

        for i in range(0,len(pdb_codes)):
            if len(pdb_codes[i])>6:
                codes.append(pdb_codes[i].replace(":", ""))
            else:
                codes.append(pdb_codes[i].replace("_", ""))


        m = modeller.automodel.loopmodel(env, alnfile= pir_name,
                    knowns=codes,
                    sequence = uknp_name,
                    assess_methods = (modeller.automodel.assess.DOPE,
                                      modeller.automodel.assess.normalized_dope,
                                      modeller.automodel.assess.GA341), 
                    root_name='Model')
        m.starting_model= 1                 # index of the first model
        m.ending_model  = nmodels           # index of the last model
                                            # (determines how many models to calculate)
        m.loop.starting_model= 1            # index of the first loop model 
        m.loop.ending_model  = nloops       # index of the last loop model 
                                            #(namely, if m.loop.starting_model = 1, it is the number of loop models computed)
        m.loop.md_level = modeller.automodel.refine.slow       # loop refinement method; this yields
                                            # models quickly but of low quality;
                                            # use refine.slow for better models
        m.final_malign3d = True
        allow_alternates=True
        m.use_parallel_job(j)
        m.make()
        #modeller.info.time_mark()

        results = []
        for x in m.outputs:
            results.append((x['name'], x['molpdf'], x['DOPE score'], x['GA341 score'][1], x['Normalized DOPE score']))
        
        df = pd.DataFrame(results, columns=['name', 'molpdf', 'DOPE score', 'GA341 score', 'Normalized DOPE score'])
        return df

    def modeller(myseq, pdb_codes, pir_name, mydir, templates_dir, cpu, nmodels = 5):
        """call Modeller to:
           1. align the target to the template sequence (using pir file format)
           2. create the models - 5 by default (results in .pdb files)"""

        #log.verbose()
        import modeller
        import modeller.parallel
        import modeller.automodel
        
        os.chdir(mydir)

        j = modeller.parallel.Job()
        for ncpu in range(cpu):
            j.append(modeller.parallel.LocalWorker())
        
        
        env = modeller.environ()
        env.io.atom_files_directory = [".", templates_dir]
        env.libs.topology.read(file = "$(LIB)/top_heav.lib")
        modeller.log.none()
        aln = modeller.alignment(env)
        for i in range(0,len(pdb_codes)):
            code, chain = utils.divide_code(pdb_codes[i])
            mdl = modeller.model(env, file = code, model_segment = ("FIRST:" + chain, "LAST:" + chain))
            aln.append_model(mdl, atom_files = code, align_codes = code + chain)

        for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                            ((1., 0.5, 1., 1., 1., 0.), False, True),
                                            ((1., 1., 1., 1., 1., 0.), True, False)):
            aln.salign(rms_cutoff = 3.5, normalize_pp_scores = False,
                    rr_file = "$(LIB)/as1.sim.mat", overhang = 30,
                    gap_penalties_1d = (-450, -50),
                    gap_penalties_3d = (0, 3), gap_gap_score = 0, gap_residue_score = 0,
                    dendrogram_file = myseq + "_templates.tree",
                    alignment_type = "tree", # If "progresive", the tree is not
                                            # computed and all structues will be
                                            # aligned sequentially to the first
                    feature_weights = weights, # For a multiple sequence alignment only
                                                # the first feature needs to be non-zero
                    improve_alignment = True, fit = True, write_fit = write_fit,
                    write_whole_pdb = whole, output = "ALIGNMENT QUALITY")

        aln.write(file = myseq + "_temp.pir", alignment_format = "PIR")

        aln.salign(rms_cutoff = 1.0, normalize_pp_scores = False,
                rr_file = "$(LIB)/as1.sim.mat", overhang = 30,
                gap_penalties_1d = (-450, -50), gap_penalties_3d = (0, 3),
                gap_gap_score = 0, gap_residue_score = 0, dendrogram_file = "1is3A.tree",
                alignment_type = "progressive", feature_weights = [0]*6,
                improve_alignment = False, fit = False, write_fit = True,
                write_whole_pdb = False, output = "QUALITY")

        # # Read aligned structure(s):
        # aln = alignment(env)
        # aln.append(file="TvLDH_temp.pir", align_codes="all")
        aln_block = len(aln)

        # Read aligned sequence(s):
        # .ali file comes from hhmakemodel.py 
        uknp_name = utils.extract_uknp_name(myseq + ".ali")
        aln.append(file = myseq + ".ali", align_codes = uknp_name)

        # Structure sensitive variable gap penalty sequence-sequence alignment:
        aln.salign(output = "", max_gap_length = 20,
                gap_function = True,   # to use structure-dependent gap penalty
                alignment_type = "PAIRWISE", align_block = aln_block,
                feature_weights = (1., 0., 0., 0., 0., 0.), overhang = 0,
                gap_penalties_1d = (-450, 0),
                gap_penalties_2d = (0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
                similarity_flag = True)

        aln.write(file = pir_name, alignment_format = "PIR")

        codes = []

        for i in range(0,len(pdb_codes)):
            if len(pdb_codes[i])>6:
                codes.append(pdb_codes[i].replace(":", ""))
            else:
                codes.append(pdb_codes[i].replace("_", ""))

        a = modeller.automodel.automodel(env, alnfile = pir_name,
                    knowns = codes, 
                    sequence = uknp_name,
                    assess_methods = (modeller.automodel.assess.DOPE,
                                    modeller.automodel.assess.normalized_dope,
                                    modeller.automodel.assess.GA341),
                    root_name='Model')
        a.md_level = modeller.automodel.refine.very_slow
        a.starting_model = 1
        a.ending_model = nmodels
        a.final_malign3d = True
        allow_alternates=True
        a.use_parallel_job(j)
        a.make()
        #modeller.info.time_mark()

        results = []
        for x in a.outputs:
            results.append((x['name'], x['molpdf'], x['DOPE score'], x['GA341 score'][1], x['Normalized DOPE score']))
        
        df = pd.DataFrame(results, columns=['name', 'molpdf', 'DOPE score', 'GA341 score', 'Normalized DOPE score'])
        return df

    def align3D(refe,r_chain, mobi, m_chain, output_directory=None):
        """Sequence-based structural alignment of two proteins."""

        from Bio.PDB import PDBParser, PDBIO, Superimposer
        from Bio.PDB.Polypeptide import is_aa

        from Bio import pairwise2
        from Bio.Align import substitution_matrices
        from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1


        def align_sequences(structA, structB):
            """
            Performs a global pairwise alignment between two sequences
            using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
            as implemented in Biopython. Returns the alignment, the sequence
            identity and the residue mapping between both original sequences.
            """

            def _get_pdb_sequence(structure):
                """
                Retrieves the AA sequence from a PDB structure.
                """

                _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, "X"))
                seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
                return seq

            resseq_A = _get_pdb_sequence(structA)
            resseq_B = _get_pdb_sequence(structB)

            sequence_A = "".join([i[1] for i in resseq_A])
            sequence_B = "".join([i[1] for i in resseq_B])
            alns = pairwise2.align.globalds(
                sequence_A,
                sequence_B,
                substitution_matrices.load("BLOSUM62"),
                one_alignment_only=True,
                open=-10.0,
                extend=-0.5,
                penalize_end_gaps=(False, False),
            )
            
            best_aln = alns[0]
            aligned_A, aligned_B, score, begin, end = best_aln

            # Equivalent residue numbering
            # Relative to reference
            mapping = {}
            aa_i_A, aa_i_B = 0, 0
            for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
                if aa_aln_A == "-":
                    if aa_aln_B != "-":
                        aa_i_B += 1
                elif aa_aln_B == "-":
                    if aa_aln_A != "-":
                        aa_i_A += 1
                else:
                    assert resseq_A[aa_i_A][1] == aa_aln_A
                    assert resseq_B[aa_i_B][1] == aa_aln_B
                    mapping[resseq_A[aa_i_A][0]] = resseq_B[aa_i_B][0]
                    aa_i_A += 1
                    aa_i_B += 1

            return mapping, score
    

        def parse_structure(filepath):
            """Parse a PDB/cif structure."""
            parser = PDBParser()
            return parser.get_structure(filepath, str(filepath))
    
        
        # Parse structures & take only the necessary chain
        s_reference = parse_structure(refe)
        try:
            reference = s_reference[0][r_chain]
        except KeyError:
            raise Exception(f"Chain {r_chain} not found in reference.")

        s_mobile = parse_structure(mobi)
        try:
            mobile = s_mobile[0][m_chain]
        except KeyError:
            raise Exception(f"Chain {m_chain} not found in mobile.")

        # Align sequences to get mapping between residues
        mapping, score = align_sequences(reference, mobile)

        refe_ca_list, mobi_ca_list = [], []
        for refe_res in mapping:
            refe_ca_list.append(reference[refe_res]["CA"])
            mobi_ca_list.append(mobile[mapping[refe_res]]["CA"])

        # Superimpose matching residues
        si = Superimposer()
        si.set_atoms(refe_ca_list, mobi_ca_list)
        si.apply(mobile.get_atoms())

        #print(f"RMSD between structures: {si.rms:4.2f}")
        
        #prepare output_filename
        output_filename = os.path.splitext(mobi)[0]+'_grid'+os.path.splitext(mobi)[1]
        if output_directory:
            
            # Write aligned mobile
            io = PDBIO()
            io.set_structure(mobile)
            io.save(os.path.join(output_directory, output_filename))
        else: 
            output_directory = homedirectory

        rot,tran = si.rotran
        rms = si.rms

        return os.path.join(output_directory, output_filename), rot, tran, rms, score

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
    
class Upload:
    from ipywidgets import FileUpload
    import codecs
    def __init__(self, **kwargs):
        
        return
    
    def receptor(self):
        print('Receptor pdb file:')
        upload = FileUpload(accept='.pdb', multiple=False)
        display(upload)
        self._ureceptor = upload
        return 
    
    
    def ligand(self):
        print('Ligands pdb files:')
        upload = FileUpload(accept='.pdb', multiple=True)
        display(upload)
        self._uligands = upload
        return
    
    
    def PARSEreceptor(self):
        if self._ureceptor:
            for i in self._ureceptor.value.keys(): 
                    filename = i
                    content = self._ureceptor.value[filename]['content']
                    with open(filename, 'w') as f:
                        f.write(codecs.decode(content))
            print('Receptor: ', filename, '\nStatus: OK\n')
        else: raise TypeError('Receptor not uploaded.')
        
        return filename
     
    def PARSEligands(self):  
        ligands_list=[]
        if self._uligands:
            for i in self._uligands.value.keys(): 
                filename = i
                content = self._uligands.value[filename]['content']
                with open(filename, 'w') as f:
                    f.write(codecs.decode(content))
                
                ligands_list.append(filename)
        else: raise TypeError('Ligands not uploaded.')
        print('Ligands: ', ligands_list, '\nStatus: OK\n')
        return ligands_list