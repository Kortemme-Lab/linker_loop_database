#!/usr/bin/env python2.7

import pyrosetta
from pyrosetta import rosetta


def get_pose_abego(pose, abego_manager):
    '''Get the ABEGO sequence for a given pose.'''
    abego_list = []
    
    for i in range(1, pose.size()):
        abego_list.append(abego_manager.index2symbol(abego_manager.torsion2index_level1(
            pose.phi(i), pose.psi(i), pose.omega(i))))

    return abego_list


if __name__ == '__main__':
    
    pyrosetta.init()
    
    #vall_path = '/home/xingjie/Softwares/Rosetta/githubRepo/tools/fragment_tools/vall.jul19.2011.gz'
    vall_path = '/home/xingjie/Softwares/Rosetta/githubRepo/main/database/sampling/small.vall.gz'

    vall_provider = rosetta.protocols.frag_picker.VallProvider()
    vall_provider.vallChunksFromLibrary(vall_path)

    abego_manager = rosetta.core.sequence.ABEGOManager()

    for i in range(1, vall_provider.size() + 1):
        chunk = vall_provider.at(i)
        pose = chunk.get_pose()
        pose.dump_pdb('test.pdb')
        exit()
        
        #print get_pose_abego(chunk.get_pose(), abego_manager)

