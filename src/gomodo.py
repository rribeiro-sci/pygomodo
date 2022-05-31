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

import os, re, subprocess
homedirectory=os.path.dirname(__file__)

HHMdatabase=os.path.join(homedirectory,'databases/pyGOMODO_db/pyGOMODO')
templatesDatabase = os.path.join(homedirectory,'databases/costum_db/wo_isoforms/gpcr_db_wo_isoforms')
mypython='python'
hhsuitePath=os.environ['HHSUITEPATH']
hhsuitescripts=os.environ['HHSUITESCRIPTS']
processedPDB=os.path.join(homedirectory,'databases/ProcessedPdbs_wo_isoforms')
GPCR_Ref = os.path.join(homedirectory,'databases/GPCR_Ref.sqlite3')

class sharedFunctions:
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
                        newseq = sharedFunctions.charPerLine(newseq,80)
                        final = defline + newseq
                        if re.match(">[^\n]*\n([GALMFWKQESPVICYHRNDT]{1,80}\n)+$", final) != None:
                            with open(output, 'w') as f: f.write(final)
                            return final
        except: raise TypeError("""input value is not in FASTA""")

    def hhblits(seqfile, output, logfilePath, hhblitsDbsPath, ncpus, rounds):
        pos = output.find(".")
        if output[pos:]==".a3m":
            ext = " -oa3m "
        else:
            ext = " -o "
        f = open(output, "w+")
        command = hhsuitePath+"hhblits -cpu " + str(ncpus) + " -i " + seqfile + " -d "+ hhblitsDbsPath  + ext + output  + " -n " + str(rounds)
        subprocess.call(command, shell=True)
        f.close()
        sharedFunctions.update_logfile(logfilePath, command, "hhblits")
        return

    def update_logfile(logname, command, step):
        f = open(logname, "a+")
        f.write("\n" + step.upper() + ": " + command + "\n")
        f.close()
        return

    def getConformation(dir):
        import pandas as pd

        filename = os.path.join(dir, 'templates.hhr')
        

        col_names=sharedFunctions.hhr_name_parser(filename)
        df = pd.DataFrame(columns=col_names)
        df = sharedFunctions.hhr_summaryhit_parser(filename, df, col_names)
        # print(df)
        df['Sequence identity'] = sharedFunctions.extract_seq_identity(filename)
        df['State'] = sharedFunctions.retrieve_state(GPCR_Ref, df)

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
        import pandas as pd
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

    def hhmakemodel(logfilePath, hhsuitePath, input_file, mydir, seqname, templates_dir, nmodel = 1):
        """create a .ali file (pir format - necessary as Modeller input)"""
        command = mypython+" "+os.path.join(hhsuitescripts, "hhmakemodel.py ")+input_file+" "+templates_dir+" "+os.path.join(mydir,seqname + ".ali")+" ./ -m "+str(nmodel)
        subprocess.call(command, shell=True)
        sharedFunctions.update_logfile(logfilePath, command, "hhmakemodel")

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
            code, chain = sharedFunctions.divide_code(pdb_codes[i])
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
        uknp_name = sharedFunctions.extract_uknp_name(myseq + ".ali")
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
        
        import pandas as pd
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
            code, chain = sharedFunctions.divide_code(pdb_codes[i])
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
        uknp_name = sharedFunctions.extract_uknp_name(myseq + ".ali")
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
        
        import pandas as pd
        df = pd.DataFrame(results, columns=['name', 'molpdf', 'DOPE score', 'GA341 score', 'Normalized DOPE score'])
        return df


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
        if 'jobname' in kwargs:
            self._jobname = kwargs.pop('jobname')
        else: 
            import datetime
            self._jobname = 'pyGOMODO_'+datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
        
        self._jobdir = os.path.join(os.getcwd(), self._jobname)

        if 'ncpus' in kwargs:
            self._ncpus = kwargs.pop('ncpus')
        else: self._ncpus=1    
        
        self._uniprotID=None
        self._fasta=None
        self._filename='sequence.seq'
        self._cwd = os.getcwd()+"/"
        self._humanDB = True
        self._rounds=1
        
        #create the job folder  if not exists
        if not os.path.isdir(self._jobdir):
            subprocess.call('mkdir ' + self._jobdir, shell=True)
        
        #create log file is not exists
        self._name = self._filename[:self._filename.find('.')]
        
        self._logfilePath = os.path.join(self._jobdir, 'output.log')
        if not os.path.isfile(self._logfilePath):
            subprocess.call('touch '+ self._logfilePath, shell=True)
 
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
        else: self._HHMdatabase=HHMdatabase


        #is sequence in fasta format? Save in .fas file
        output_fas = os.path.join(self._jobdir, 'sequence.fas')
        self._fasta80 = sharedFunctions.searchFasta(self._fasta, output_fas)

        #Run HHBLITS
        print('Running...\n')
        sharedFunctions.hhblits(output_fas, os.path.join(self._jobdir, "query.a3m"), self._logfilePath, self._HHMdatabase, self._ncpus, self._rounds)
       

        def parsehhr():
            import pandas as pd

            filename = os.path.join(self._jobdir,'sequence.hhr')
            col_names=sharedFunctions.hhr_name_parser(filename)
            df = pd.DataFrame(columns=col_names)
            df = sharedFunctions.hhr_summaryhit_parser(filename, df, col_names)
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
        else: self._TemplatesDatabase=templatesDatabase

        #Run HHBLITS
        if self._hhmprofile:
            print('Running...\n')
            sharedFunctions.hhblits(self._hhmprofile, os.path.join(self._jobdir + "templates.hhr"), self._logfilePath, self._TemplatesDatabase, self._ncpus, self._rounds)
        else:
            print('Running...\n')
            sharedFunctions.hhblits(os.path.join(self._jobdir, "query.a3m"), os.path.join(self._jobdir, "templates.hhr"), self._logfilePath, self._TemplatesDatabase, self._ncpus, self._rounds)

        df =  sharedFunctions.getConformation(self._jobdir)
        self._templates=df
        print('Done!\n')
        return df

    def makeModels(self, **kwargs):
               
        if 'loop' in kwargs:
            self._loop = kwargs.pop('loop')
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
        else: self._hhsuitePath=hhsuitePath
        
        if 'processedPDB' in kwargs:
            self._processedPDB = kwargs.pop('processedPDB')
        else: self._processedPDB=processedPDB

        #create the folder "templates" if not exists
        templates_dir = os.path.join(self._jobdir, 'templates')
        if not os.path.isdir(templates_dir):
            subprocess.call('mkdir ' + templates_dir, shell=True)

        pir_name = os.path.join(templates_dir,  "mytemplate.pir")
     
        #myseq = sharedFunctions.extract_seqname(os.path.join(self._jobdir , "templates.hhr"))
        # write target sequence in pir format
        sharedFunctions.hhmakemodel(self._logfilePath, self._hhsuitePath, os.path.join(self._jobdir, "templates.hhr"), self._jobdir, 'sequence', templates_dir)


        #templates choosen by the user
        if 'templates' in kwargs:
            self._templates = kwargs.pop('templates')
        else: raise TypeError('Please select the templates')
        
    
        #copy the pdb (if does not exists yet) in the "templates" folder
        for i in range(0,len(self._templates)):
            pdb, chain = sharedFunctions.divide_code(self._templates[i])
            if not os.path.isfile(os.path.join(templates_dir, pdb + ".pdb")):
                if self._humanDB:
                    subprocess.call("cp "+os.path.join(self._processedPDB, pdb + "_proc.pdb ") + os.path.join(templates_dir, pdb + ".pdb"), shell=True)
                    #print("cp "+os.path.join(self._processedPDB, pdb + "_proc.pdb ") + os.path.join(templates_dir, pdb + ".pdb"))
                else:
                    sharedFunctions.download_pdb(pdb, self._cwd + "templates")

        #modeling (MODELLER)
        if self._loop:
            #print('*'*100 + '\n loop')
            models = sharedFunctions.modellerLoop('sequence', self._templates, pir_name, self._cwd, templates_dir, self._ncpus, self._nmodels, self._nloops)
        else:
            #print('*'*100 + '\n automodel')
            models = sharedFunctions.modeller('sequence', self._templates, pir_name, self._jobdir, templates_dir, self._ncpus, self._nmodels)

        return models

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
