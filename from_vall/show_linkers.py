#!/usr/bin/env python2

import os
import json

import pyrosetta
from pyrosetta import rosetta


def create_pose(size):
    '''Create a pose for a given size.'''
    pose = rosetta.core.pose.Pose()
    rosetta.core.pose.make_pose_from_sequence(pose, 'A' * size, 'fa_standard')
    return pose

def show_linkers(linker_db, length, num_output=100, output_path='debug'):
    '''Dump linkers into PDB files'''
    with open(linker_db, 'r') as f:
        linkers = json.load(f)

    print 'Load', len(linkers), 'linkers from the database', linker_db

    pose = create_pose(length)

    for i, l in enumerate(linkers):
        for j in range(length):
            pose.set_phi(j + 1, l['phis'][j])
            pose.set_psi(j + 1, l['psis'][j])
            pose.set_omega(j + 1, l['omegas'][j])

        #print i, l['pdb_id'], l['start_position'], l['sequence']
        if i >= num_output: break
        pose.dump_pdb(os.path.join(output_path, '{0}.pdb'.format(i)))


if __name__ == '__main__':
    pyrosetta.init()

    #show_linkers('linker_sheet_helix_3_BAB_with_padding.json', 11, 100)
    
    show_linkers('linker_helix_sheet_4_non_redundant.json', 6, 500)
