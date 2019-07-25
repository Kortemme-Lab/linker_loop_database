#!/usr/bin/env python2

import os
import json

import numpy as np

import pyrosetta
from pyrosetta import rosetta


def create_pose(size):
    '''Create a pose for a given size.'''
    pose = rosetta.core.pose.Pose()
    rosetta.core.pose.make_pose_from_sequence(pose, 'A' * size, 'fa_standard')
    return pose

def apply_linker_to_pose(linker, pose):
    '''Apply a linker to a pose.'''
    for i in range(len(linker['phis'])):
        pose.set_phi(i + 1, linker['phis'][i])
        pose.set_psi(i + 1, linker['psis'][i])
        pose.set_omega(i + 1, linker['omegas'][i])

def RMSD(P1, P2):
    '''Calculate the RMSD between two lists of xyz vectors'''
    d2s = [P1[i].distance_squared(P2[i]) for i in range(len(P1))]
    return np.sqrt(sum(d2s) / len(d2s)) 

def BB_RMSD(pose1, pose2):
    '''Calculate the backbone RMSD between two poses.'''
    bb_atoms = ['N', 'CA', 'C']
    P1 = [pose1.residue(i).xyz(a) for i in range(1, pose1.size() + 1) for a in bb_atoms]
    P2 = [pose2.residue(i).xyz(a) for i in range(1, pose2.size() + 1) for a in bb_atoms]

    return RMSD(P1, P2)

def c_term_RMSD(pose1, pose2):
    '''Calculate the C-terminal RMSD between two poses.'''
    bb_atoms = ['N', 'CA', 'C']
    P1 = [pose1.residue(pose1.size()).xyz(a) for a in bb_atoms]
    P2 = [pose2.residue(pose2.size()).xyz(a) for a in bb_atoms]

    return RMSD(P1, P2)

def filter_linkers_by_RMSD(linker_db, all_rmsd_cutoff=2, terminal_rmsd_cutoff=1.5):
    '''Filter the linkers by RMSD. Two linkers and considered
    redundant if the all backbone RMSD is within all_rmsd_cutoff
    and the C-term residue backbone RMSD is within terminal_rmsd_cutoff.
    '''
    # Load the database

    with open(linker_db, 'r') as f:
        linkers = json.load(f)
    
    linker_length = len(linkers[0]['phis'])

    print 'Load', len(linkers), 'linkers from the database', linker_db
    print 'The linker length is', linker_length 
  
    # Select the linkers

    pose_added = create_pose(linker_length)
    pose_new = create_pose(linker_length)

    selected_linkers = []
    cluster_sizes = []

    for i, l in enumerate(linkers):
        apply_linker_to_pose(l, pose_new)

        redundant = False
        for j, l_old in enumerate(selected_linkers):
            apply_linker_to_pose(l_old, pose_added)

            if BB_RMSD(pose_new, pose_added) < all_rmsd_cutoff \
                    and c_term_RMSD(pose_new, pose_added) < terminal_rmsd_cutoff:
                redundant = True
                cluster_sizes[j] += 1

                #print 'Found redundant linker:\n', l, '\n', l_old, '\n'###DEBUG
                #pose_new.dump_pdb('debug/pose_new.pdb')###DEBUG
                #pose_added.dump_pdb('debug/pose_added.pdb')###DEBUG
                #exit()###DEBUG
                
                break

        if not redundant:
            selected_linkers.append(l)
            cluster_sizes.append(1)

        if i % 100 == 0:
            print '{0}/{1} linkers tested. Found {2} non-redundant linkers.'.format(i, len(linkers), len(selected_linkers))
            
    print 'Selected', len(selected_linkers), 'non-redundant linkers'
    
    # Dump the selected linkers

    with open(linker_db[:-5] + '_non_redundant.json', 'w') as f:
        json.dump(selected_linkers, f)
   
    with open(linker_db[:-5] + '_cluster_sizes.json', 'w') as f:
        json.dump(cluster_sizes, f)

    return selected_linkers


if __name__ == '__main__':
    pyrosetta.init()
    
    filter_linkers_by_RMSD('linker_helix_sheet_2.json')
    filter_linkers_by_RMSD('linker_helix_sheet_3.json')
    filter_linkers_by_RMSD('linker_helix_sheet_4.json')
    filter_linkers_by_RMSD('linker_helix_sheet_5.json')
    
    filter_linkers_by_RMSD('linker_sheet_helix_2.json')
    filter_linkers_by_RMSD('linker_sheet_helix_3.json')
    filter_linkers_by_RMSD('linker_sheet_helix_4.json')
    filter_linkers_by_RMSD('linker_sheet_helix_5.json')
    
    filter_linkers_by_RMSD('linker_helix_helix_2.json')
    filter_linkers_by_RMSD('linker_helix_helix_3.json')
    filter_linkers_by_RMSD('linker_helix_helix_4.json')
    filter_linkers_by_RMSD('linker_helix_helix_5.json')

    filter_linkers_by_RMSD('linker_sheet_sheet_2.json')
    filter_linkers_by_RMSD('linker_sheet_sheet_3.json')
    filter_linkers_by_RMSD('linker_sheet_sheet_4.json')
    filter_linkers_by_RMSD('linker_sheet_sheet_5.json')
