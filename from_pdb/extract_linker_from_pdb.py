#!/usr/bin/env python2.7

import pyrosetta
from pyrosetta import rosetta


def extract_linker_from_pdb(pdb_file, chain, start, stop):
    '''Extract the torsions of a linker from a PDB file,
    the start and stop points are defined in pdb residue ids.
    '''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)

    start_seqpos = pose.pdb_info().pdb2pose(chain, start)
    stop_seqpos = pose.pdb_info().pdb2pose(chain, stop)

    torsions = []

    for i in range(start_seqpos, stop_seqpos + 1):
        torsions.append((pose.phi(i), pose.psi(i), pose.omega(i)))

    print torsions


if __name__ == '__main__':
    pyrosetta.init()

    extract_linker_from_pdb('inputs/2igd.pdb', 'A', 25, 32) 
    #extract_linker_from_pdb('inputs/2igd.pdb', 'A', 40, 47) 

