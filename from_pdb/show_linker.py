#!/usr/bin/env python2.7

import pyrosetta
from pyrosetta import rosetta

def show_linker(start_position, torsions):
    '''Insert a linker to the test sheet pose.'''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, 'inputs/test_sheet.pdb')
    
    ft = rosetta.core.kinematics.FoldTree()
    ft.add_edge(1, pose.size(), 1)
    ft.add_edge(1, 22, -1)
    ft.add_edge(pose.size(), 23, -1)

    pose.fold_tree(ft)
    rosetta.core.pose.correctly_add_cutpoint_variants(pose)


    for i, t in enumerate(torsions):
        pose.set_phi(start_position + i, t[0])
        pose.set_psi(start_position + i, t[1])
        pose.set_omega(start_position + i, t[2])

    pose.dump_file('test.pdb')


if __name__ == '__main__':
    pyrosetta.init()

    show_linker(7, [(-155.5666755419531, 161.44582156293308, 176.75106770275963), (-73.99647537198591, -24.068561633329338, -177.7598639587701), (-154.52193322619232, 168.75049876571543, -175.3235214934199), (-70.28629540101745, -32.87295199000924, -179.64596127188162), (-64.28346055429476, -39.50355723851481, 178.50692921320442), (-65.53359808192624, -44.032788075995526, 175.37957549084007), (-61.31512409913714, -43.10189690164581, -179.0552298144351), (-58.29029769682966, -43.916979339580045, 179.1890275600125)])
