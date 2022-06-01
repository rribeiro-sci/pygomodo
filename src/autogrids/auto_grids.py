import sys
import os
import array
import numpy
import Bio.SVDSuperimposer
from Bio                        import pairwise2
from Bio.pairwise2              import format_alignment
from Bio.SubsMat.MatrixInfo     import blosum62
from Bio.PDB.Superimposer       import Superimposer
from Bio.PDB.PDBIO              import PDBIO
from Bio.PDB.PDBParser          import PDBParser
from Bio.PDB.vectors            import Vector


def get_seq(res_list):

    letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
               'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
               'TYR':'Y','VAL':'V'}
    residues = []
    seq = ""

    for res in res_list:
        if res.id[0] == ' ':
            residues.append(letters[res.get_resname()])

    return seq.join(residues)

def get_chains_lists(pdb_filename, structure_id):

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(structure_id, pdb_filename)
    #Structure must be composed of exactly one model!
    model = structure[0]
    atoms_chains = []
    res_chains = []
    for chain in model :
        atoms_chains.append([])
        res_chains.append([])

    i = 0
    for chain in model:
        for residue in chain:
            res_chains[i].append(residue)

            for atom in residue:

                atoms_chains[i].append(atom)
        i += 1
    return atoms_chains, res_chains, structure

def get_chain_id(atoms_list):

    id = atoms_list[0].get_parent().get_parent().get_id()
    return id

def superimpose_get_rotranrms(fixed, moving):

    # 'fixed' and 'moving' are lists of Atom objects.
    # The moving atoms will be put on the fixed atoms.

    Superimposer

    sup = Superimposer()
    sup.set_atoms(fixed, moving)

    # calculate rot and tran matrices.

    rot,tran = sup.rotran
    rms = sup.rms


    return rms, rot, tran

def allign_3D(ref_list_atoms, ref_list_res, targ_list_res):

    seq1 = get_seq(targ_list_res)
    seq2 = get_seq(ref_list_res)

    alignment = pairwise2.align.globalds(seq1, seq2, blosum62, -10, -.5)[0]
    score = alignment[2]

    #set the two atoms lists to use to find the transformation matrix
    ref_list = []
    targ_list = []

    #create a list of indexes iterating over alignment, where each element contains a couple of value
    #(start and stop of the continuous match in the alignment)

    start_tar = 0
    stop_tar = 0
    start_ref = 0
    stop_ref = 0
    indexes_tar =[]
    indexes_ref = []
    continuous = False
    count_tar = 0
    count_ref = 0

    max = 0

    for i in range(len(alignment[1])):

        if (alignment[0])[i] != '-' :
            count_tar += 1

        if (alignment[1])[i] != '-' :
            count_ref += 1

        if (alignment[0])[i] != '-' and (alignment[1])[i] != '-' :
            if continuous == False:
                start_tar = count_tar
                start_ref = count_ref
                continuous = True
        else:
            if continuous == True:
                stop_tar = count_tar
                indexes_tar.append([start_tar,stop_tar])
                stop_ref = count_ref
                indexes_ref.append([start_ref,stop_ref])
                continuous = False
                if (stop_tar - start_tar) > max :
                    max = stop_tar-start_tar

    #check if the alignment produced a perfect match. if true perform the superimposition on the whole atom lists,
    #otherwise perform the superimposition on the atoms leading to continuous-matching subsequences of the alignment:
    #set k as the minimum length of the alignment's continuous subsequences you want to consider.
    #NB! max is used to fix the maximum value that k can assume, and refears to the maximum continuous subsequence length.

    if len(indexes_tar) == 0 and count_ref == len(ref_list_res):

        for res in targ_list_res:
            for atom in res:
                targ_list.append(atom)

        ref_list = ref_list_atoms

    else:
        k = max -10

        #extract from the target structure the atoms leading to continuous match
        #subsequences of the alignment with length >= keys

        for index_tar, index_ref in zip(indexes_tar, indexes_ref):

            if index_tar[1]-index_tar[0] >= k :
                for i in range(index_tar[0],index_tar[1]):

                    for atom in targ_list_res[i]:
                        targ_list.append(atom)

                for i in range(index_ref[0],index_ref[1]):

                    for atom in ref_list_res[i]:
                        ref_list.append(atom)

    #resize the two lists to perform superimposition

    if len(targ_list)>=len(ref_list):
        targ_list = targ_list[:len(ref_list)]

    else:

        ref_list = ref_list[:len(targ_list)]

    #try superimpose
    rms,rot,tran = superimpose_get_rotranrms(targ_list,ref_list)

    return rms, rot, tran, score

def auto_grids(receptor_filename, family, reference_files_dir):

    # NEW FAMILIES SHOULD BE ADDED HERE!

    families = {
                'IGR' : [reference_files_dir+'IGr_Ref.pdb', reference_files_dir+'IGr_Ref_grids.txt',],
                'GPCR' : [reference_files_dir+'GPCr_Ref.pdb', reference_files_dir+'GPCr_Ref_grids.txt',]
    }

    if family not in families:
        print('Argouments Error: family argument should be one of: ', families.keys())

    else:
        reference_structures_filename = (families[family])[0]
        reference_grids_filename = (families[family])[1]


    targ_atoms_lists, targ_res_lists, s1 = get_chains_lists(receptor_filename, "tar")
    ref_atoms_lists, ref_res_lists, s2 = get_chains_lists(reference_structures_filename, "ref")
    n_targ_chains = len(s1[0])
    n_ref_chains = len(s2[0])

    grids_file_ref = open(reference_grids_filename, "r")

    #create a dictionary containing the coordinates of the reference grids (as vectors) which can be
    #modified during the execution in order to not to modify the original reference grids file,
    #and a dictionary that will contain the transformed coordinates of the grids.

    ref_grids_dic = {}
    for line in grids_file_ref :

        toks = line.split()
        chain_id = toks[0]
        coor = Vector(toks[2],toks[3],toks[4])
        ref_grids_dic.update({chain_id : coor})

    targ_grids_dic = {}

    #Start the detection
    for i in range(n_targ_chains):

        #set a reasonable high value to minimize the rmsd
        score_opt = 0
        rms_opt = 1000
        targ_res_list = targ_res_lists[i]

        for j in range(n_ref_chains):

            #create copies in order to preserve the integrity of the original chains
            ref_atoms_list = ref_atoms_lists[j]
            ref_res_list = ref_res_lists[j]

            rms,rot,tran,score = allign_3D(ref_atoms_list,ref_res_list,targ_res_list)

            if rms < rms_opt:
                rms_opt = rms
                rot_opt, tran_opt = rot, tran
                score_opt = score
                opt = j

        targ_chain_id = get_chain_id(targ_atoms_lists[i])
        ref_chain_id = get_chain_id(ref_atoms_lists[opt])

        #set a threshold for the alignment score to judge the goodness of the calculated alignment
        if score_opt < 80 and rms_opt > 7:

            print('Error: no good structural alignment found for chain ',targ_chain_id)
            print(score_opt,rms_opt)
        else:

            #read the reference coordinates file and transform the relative coordinates currently stored in the dictionary,
            #then write them to the output target's grids coordinates file.
            #NB: this is important , because all the transformations performed of the reference chains change the
            #position of the chains themselves, but not the relative grids coordinates!

            ref_grid_coor = ref_grids_dic[ref_chain_id]
            targ_grid_coor = Vector.right_multiply(ref_grid_coor, rot_opt) + tran_opt

            #store the transformed coordinates in the target grids dictionary
            targ_grids_dic.update({targ_chain_id : targ_grid_coor})

            #print summary
            print("###############################################################")
            print("                                                               ")
            print("Target chain '"+targ_chain_id+"' Summary :                     ")
            print("                                                               ")
            print("reference chain used : '"+ref_chain_id+"'"                       )
            print("calculated rmsd :",rms_opt,                                     )
            print("calculated score alignment :",score_opt,                        )
            print("grid center coordinates : ",targ_grid_coor,                     )
            print("                                                               ")
            print("###############################################################")
            print("                                                               ")

    #return a dictionary containing the coordinates for each chain of the structure
    return targ_grids_dic
