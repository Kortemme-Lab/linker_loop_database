#!/usr/bin/env python2.7

import pyrosetta
from pyrosetta import rosetta


def mutate_residues(pose, res_list, aa_list):
    '''Mutate a list of residues. The list of AAs could
    either be 1 letter code or 3 letter code.
    '''
    aa_name_map = {'A':'ALA', 'P':'PRO', 'V':'VAL', 'L':'LEU', 'I':'ILE', 'M':'MET',
                   'F':'PHE', 'Y':'TYR', 'W':'TRP', 'S':'SER', 'T':'THR', 'C':'CYS',
                   'K':'LYS', 'R':'ARG', 'H':'HIS', 'D':'ASP', 'E':'GLU', 'N':'ASN',
                   'Q':'GLN', 'G':'GLY'}


    mutater = rosetta.protocols.simple_moves.MutateResidue()
    for i in range(len(res_list)):
        name = aa_list[i] if len(aa_list[i]) == 3 else aa_name_map[aa_list[i]]
        mutater.set_res_name(name)
        mutater.set_target(res_list[i])
        mutater.apply(pose)

if __name__ == '__main__':
    pyrosetta.init()
    
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, 'inputs/structure_0.pdb')

    ft = rosetta.core.kinematics.FoldTree()
    ft.add_edge(1, pose.size(), 1)
    ft.add_edge(1, 22, -1)
    ft.add_edge(pose.size(), 23, -1)

    pose.fold_tree(ft)
    rosetta.core.pose.correctly_add_cutpoint_variants(pose)

    for i in range(8, 39):
        rosetta.core.conformation.idealize_position(i, pose.conformation())
        pose.set_phi(i, -120)
        pose.set_psi(i, 117)
        pose.set_omega(i, 180)

    mutate_residues(pose, list(range(1, pose.size() + 1)), ['ALA'] * pose.size())

    pose.dump_file('inputs/test_sheet.pdb')
