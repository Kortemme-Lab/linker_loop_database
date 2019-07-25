#!/usr/bin/env python2.7
import re
import json

import pyrosetta
from pyrosetta import rosetta


def get_chunk_abego(chunk):
    '''Get the ABEGO sequence for a vall chunk.'''
    abego_manager = rosetta.core.sequence.ABEGOManager()
    abego_list = []

    for i in range(1, chunk.size() + 1):
        abego_list.append(abego_manager.index2symbol(abego_manager.torsion2index_level1(
            chunk.at(i).phi(), chunk.at(i).psi(), chunk.at(i).omega())))

    return ''.join(abego_list)

def get_chunk_ss(chunk):
    '''Get the secondary structure sequence of the chunk.'''
    return ''.join([chunk.at(i).ss() for i in range(1, chunk.size() + 1)])

def get_matched_positions(target, query):
    '''Get the positions starting from which the query regex matches 
    the target substring.
    '''
    positions = [m.start() for m in re.finditer(query, target)]
    return positions

def find_linker_loops(chunk, preSS, postSS, linker_pattern, abego_pattern):
    '''Find the sequence and torsions of the linker loop given the 
    definition of the preSS and postSS and linker length. 
    Return the list of linker loops. Note that the linker loop also
    contains the torsion of the stubs, which means each fragment is of
    length linker_length + 2'''
    sequence = chunk.get_sequence()
    ss = get_chunk_ss(chunk)
    if abego_pattern:
        abego = get_chunk_abego(chunk)
    linker_length = len(linker_pattern)

    regex = preSS + linker_pattern + postSS
    positions = get_matched_positions(ss, regex)

    linker_loops = []

    # Get the residue positions. Notice that the residues are stored in a vector1
    
    res_positions = [p + len(preSS) for p in positions]

    for p in res_positions:
        
        if abego_pattern:
            abego_seq = abego[p:p + len(abego_pattern)]
            if None is re.match(abego_pattern, abego_seq):
                continue

        phis = [chunk.at(i).phi() for i in range(p, p + linker_length + 2)]
        psis = [chunk.at(i).psi() for i in range(p, p + linker_length + 2)]
        omegas = [chunk.at(i).omega() for i in range(p, p + linker_length + 2)]
        l_seq = [chunk.at(i).aa() for i in range(p, p + linker_length + 2)]
        
        linker_loops.append({'phis' : phis, 'psis' : psis, 'omegas' : omegas, 'sequence' : l_seq,
            'pdb_id' : chunk.get_pdb_id(), 'start_position': p})

    return linker_loops

def make_linker_loop_database(linker_pattern, preSS_pattern, postSS_pattern, output_name, abego_pattern=None):
    '''Make a linker loop database.'''
    vall_path = '/home/xingjie/Softwares/Rosetta/githubRepo/tools/fragment_tools/vall.jul19.2011.gz'
    #vall_path = '/home/xingjie/Softwares/Rosetta/githubRepo/main/database/sampling/small.vall.gz'

    vall_provider = rosetta.protocols.frag_picker.VallProvider()
    vall_provider.vallChunksFromLibrary(vall_path)

    linker_loops = []

    for i in range(1, vall_provider.size() + 1):
        chunk = vall_provider.at(i)
       
        linker_loops += find_linker_loops(chunk, preSS_pattern, postSS_pattern,linker_pattern, abego_pattern)

   
    print 'Find', len(linker_loops), ' linker loops.'
    
    with open(output_name, 'w') as f:
         json.dump(linker_loops, f)

if __name__ == '__main__':
    
    pyrosetta.init()
   
    make_linker_loop_database('.' * 2, 'EEE', 'EEE', 'linker_sheet_sheet_2.json')
    make_linker_loop_database('.' * 3, 'EEE', 'EEE', 'linker_sheet_sheet_3.json')
    make_linker_loop_database('.' * 4, 'EEE', 'EEE', 'linker_sheet_sheet_4.json')
    make_linker_loop_database('.' * 5, 'EEE', 'EEE', 'linker_sheet_sheet_5.json')
    
    #make_linker_loop_database('.' * 2, 'EEE', 'HHHHHH', 'linker_sheet_helix_2.json')
    #make_linker_loop_database('.' * 3, 'EEE', 'HHHHHH', 'linker_sheet_helix_3.json')
    #make_linker_loop_database('.' * 4, 'EEE', 'HHHHHH', 'linker_sheet_helix_4.json')
    #make_linker_loop_database('.' * 5, 'EEE', 'HHHHHH', 'linker_sheet_helix_5.json')
    
    #make_linker_loop_database('.' * 2, 'HHHHHH', 'EEE', 'linker_helix_sheet_2.json')
    #make_linker_loop_database('.' * 3, 'HHHHHH', 'EEE', 'linker_helix_sheet_3.json')
    #make_linker_loop_database('.' * 4, 'HHHHHH', 'EEE', 'linker_helix_sheet_4.json')
    #make_linker_loop_database('.' * 5, 'HHHHHH', 'EEE', 'linker_helix_sheet_5.json')
    
    #make_linker_loop_database('.' * 2, 'HHHHHH', 'HHHHHH', 'linker_helix_helix_2.json')
    #make_linker_loop_database('.' * 3, 'HHHHHH', 'HHHHHH', 'linker_helix_helix_3.json')
    #make_linker_loop_database('.' * 4, 'HHHHHH', 'HHHHHH', 'linker_helix_helix_4.json')
    #make_linker_loop_database('.' * 5, 'HHHHHH', 'HHHHHH', 'linker_helix_helix_5.json')
    
    #make_linker_loop_database('EE....HHHHH', 'E', 'H', 'linker_sheet_helix_4_with_padding.json')
    #make_linker_loop_database('EE...HHHH', 'E', 'H', 'linker_sheet_helix_3_BAB_with_padding.json', abego_pattern='BBBABAAAA')
