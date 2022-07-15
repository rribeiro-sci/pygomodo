#!/usr/bin/env python



def align3D(refe,r_chain, mobi, m_chain, output_directory):
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

        return mapping
   

    def parse_structure(filepath):
        """Parse a PDB/cif structure."""
        parser = PDBParser()
        return parser.get_structure(filepath.name, str(filepath))
   
    
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
    mapping = align_sequences(reference, mobile)

    refe_ca_list, mobi_ca_list = [], []
    for refe_res in mapping:
        refe_ca_list.append(reference[refe_res]["CA"])
        mobi_ca_list.append(mobile[mapping[refe_res]]["CA"])

    # Superimpose matching residues
    si = Superimposer()
    si.set_atoms(refe_ca_list, mobi_ca_list)
    si.apply(mobile.get_atoms())

    #print(f"RMSD between structures: {si.rms:4.2f}")

    # Write aligned mobile
    io = PDBIO()
    io.set_structure(mobile)
    io.save(output_directory + filename)
