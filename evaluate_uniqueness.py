#!/usr/bin/python

"""
Use the per-PDB-sequence domain locations (in addition to binding information found in the previously
generated FASTA file, create_fasta.py) to create alignments of each domain--ligand pair and evaluate 
per-sequence "uniqueness" as in:

Henikoff S and Henikoff JG (1994). "Position-based sequence weights." J Mol Biol 243(4): 574-578. 
  doi:10.1016/0022-2836(94)90032-9

Contact Shilpa Nadimpalli Kobren (snadimpa@princeton.edu) with questions
"""

import os
import sys
import gzip
import argparse
from subprocess import call, check_output, Popen, PIPE
from difflib import SequenceMatcher


########################################################################################################
# CONSTANTS
########################################################################################################

# path to where this script is currently located (and to where all data should be stored) -- this can
# be updated
DATAPATH = os.path.dirname(os.path.abspath(__file__))+'/'

DISTANCE_CUTOFF = 20.  # cutoff (in Angstroms) to consider values for the "mindist" score
PROXIMITY_CUTOFF = 3.6  # cutoff (in Angstroms) to consider an atom-to-atom distance as "close"

# full path to a tab-delimited file with columns PDB ID-PDB Chain, Domain Name (unique), and
#  comma-delimited list of (1-indexed domain match state : 0-indexed sequence position : amino acid value)
DOMAINS = DATAPATH+'processed_data/domains/BioLiP_2017-06-28-domains-pfam_v31.tsv.gz'


########################################################################################################
# DOMAIN-LIGAND ALIGNMENTS
########################################################################################################

def line_count_from_file(infile_name):
  """
  :param infile_name: full path to a file 
  :return: the number of lines in the specified file
  """

  if not os.path.isfile(infile_name):
    sys.stderr.write('No such file: '+infile_name+'\n')
    return 0

  if not infile_name.endswith('gz'):
    word_count = check_output(['wc', '-l', infile_name])
    return int(word_count.split()[0])

  decompress_file = Popen(('zcat', infile_name), stdout=PIPE)
  word_count = check_output(['wc', '-'], stdin=decompress_file.stdout)
  decompress_file.wait()
  return int(word_count.split()[0])


########################################################################################################

def ligand_groups(ligand_group_file=DATAPATH+'downloaded_data/ligand_groups.txt'):
  """
  :param ligand_group_file: full path to a tab-delimited file matching ligand super group identifiers
                            to original mmCIF ligand identifiers
  :return: mapping from original ligand identifier -> list of all super groups it is part of
  """

  if not os.path.isfile(ligand_group_file):
    sys.stderr.write('Could not open '+ligand_group_file+'\n')
    return {}

  ligand_to_group = {}  # ligand identifier -> [list of member super groups]

  ligand_handle = gzip.open(ligand_group_file) if ligand_group_file.endswith('gz') else open(ligand_group_file)
  for lig_line in ligand_handle:
    if lig_line.startswith('#'):
      continue

    ligand_group, original_id = lig_line[:-1].split('\t')

    if original_id not in ligand_to_group:
      ligand_to_group[original_id] = set()
    ligand_to_group[original_id].add(ligand_group)
  ligand_handle.close()

  # convert sets back to lists (without redundancies)
  for original_id in ligand_to_group.keys():

    # a single molecule must belong exclusively to the nucleic acid, ion, or small molecule groups
    if 'NUCACID_' in ligand_to_group[original_id]:
      for other_group in ['ION_', 'METABOLITE_', 'DRUGLIKE_', 'SM_']:
        ligand_to_group[original_id].discard(other_group)
    elif 'ION_' in ligand_to_group[original_id]:
      for other_group in ['METABOLITE_', 'DRUGLIKE_', 'SM_']:
        ligand_to_group[original_id].discard(other_group)

    ligand_to_group[original_id] = sorted(list(ligand_to_group[original_id]))

  return ligand_to_group


########################################################################################################

def translate_ligand(ligand_name):
  """
  :param ligand_name: we may want to rename certain ligands. For instance, we rename the nucleic acids groups
                      we generated while parsing PDB files
  :return: a "new" name of the ligand (or the same name, if no translation was needed
  """

  translations = {'NUCDNA': 'DNABASE_',
                  'NUCDNAB': 'DNABACKBONE_',
                  'NUCRNA': 'RNABASE_',
                  'NUCRNAB': 'RNABACKBONE_',
                  'III': 'PEPTIDE_'}

  return translations.get(ligand_name, ligand_name)


########################################################################################################

def union_ordered_lists(seq1, seq2):
  """
  :param seq1: ordered list of characters/integers to compare
  :param seq2: ordered list of characters/integers to compare
  :return: the union of unique characters/integers from the two sequences, in the same input order
  """

  sm = SequenceMatcher(a=seq1, b=seq2)
  union_entries = []  # final ordered list of unique entries between seq1 and seq2
  for (op, start1, end1, start2, end2) in sm.get_opcodes():
    if op == 'equal' or op == 'delete':
      # This range appears in both sequences, or only in the first one.
      union_entries += seq1[start1:end1]
    elif op == 'insert':
      # This range appears only in the second sequence.
      union_entries += seq2[start2:end2]
    elif op == 'replace':
      # There are different ranges in each sequence - add both.
      union_entries += seq1[start1:end1]
      union_entries += seq2[start2:end2]

  return union_entries


########################################################################################################

def combine_match_states(match_state_vectors):
  """
  :param match_state_vectors: list of lists of strings indicating match state positions
  :return: ordered list of all match states (including those present/absent in a subset of lists)
  """

  overall_match_states = []  # ordered list of complete match states

  for current_match_state_list in match_state_vectors:
    overall_match_states = union_ordered_lists(overall_match_states,
                                               current_match_state_list)
  return overall_match_states


########################################################################################################

def format_fasta_aln(matchstate_file, alignment_out):
  """
  :param matchstate_file: full path to a tab-delimited file containing a complete sequence identifier
                          and comma-delimited list of corresponding match states
  :param alignment_out: full path to a file to write the properly-formatted FASTA alignment
  :return: the name of the alignment output file (as input) upon successful completion, None otherwise
  """

  if not os.path.isfile(matchstate_file):
    sys.stderr.write('Could not open ' + matchstate_file + '\n')
    return None
  if not os.path.isdir('/'.join(alignment_out.split('/')[:-1])):
    sys.stderr.write('Could not write to ' + alignment_out + '\n')
    return None

  seqid_to_matchstate_val = {}  # sequence ID -> match state -> amino acid residue
  seqid_to_matchstate_pos = {}  # sequence ID -> match state -> 1-index AA position,...
  all_matchstates = []  # [[1,2,3,...], [2,3,4,...]]

  # process each match state line: seqID \t matchstate:aa_position:aa_value
  matchstate_handle = gzip.open(matchstate_file) if matchstate_file.endswith('gz') else open(matchstate_file)
  for mline in matchstate_handle:
    seqid, matchstate_orig = mline[:-1].split('\t')
    matchstates = [entry.split(':') for entry in matchstate_orig.split(',')]  # split apart the match state string

    seqid_to_matchstate_val[seqid] = {curr_mstate: aa_val for (curr_mstate, _, aa_val) in matchstates}
    seqid_to_matchstate_pos[seqid] = {curr_mstate: str(int(pos)+1) for (curr_mstate, pos, _) in matchstates}
    all_matchstates.append([curr_mstate for (curr_mstate, _, _) in matchstates])
  matchstate_handle.close()

  # get the overall list of ordered match states
  ordered_matchstates = combine_match_states(all_matchstates)

  # and write out the new FASTA-formatted alignment:
  aln_handle = gzip.open(alignment_out, 'w') if alignment_out.endswith('gz') else open(alignment_out, 'w')
  for seqid in sorted(seqid_to_matchstate_val.keys()):
    description = ','.join([mstate+':'+seqid_to_matchstate_pos[seqid].get(mstate, '-')
                            for mstate in ordered_matchstates])
    aln_handle.write('>' + seqid + '\t' + str(description) + '\n')
    aln_handle.write(''.join([seqid_to_matchstate_val[seqid].get(mstate, '-')
                              for mstate in ordered_matchstates]) + '\n\n')
  aln_handle.close()

  # return the intended filename to signify successful completion
  return alignment_out


########################################################################################################

def create_alignment_files(domain_file, fasta_dir, alignment_dir, distance='mindist'):
  """
  :param domain_file: full path to a tab-delimited list of domain matches, with columns corresponding to
                      [0] PDB identifier, [1] domain name, [2] ','-separated list of match state:index:value
  :param fasta_dir: full path to a directory containing FASTA files, where the header of each FASTA
                              file has been updated with the amino acid position, ligand type, and binding score
  :param alignment_dir: full path to a directory where the domain alignments should be stored
  :param distance: string corresponding to the distance score metric (for finding and updating filenames)
  :return: list of domain names, ligand types, and corresponding paths to alignment files
  """

  # confirm that the directories we are interested in exist
  for current_dir in [fasta_dir, alignment_dir]:
    if not os.path.isdir(alignment_dir):
      sys.stderr.write('Could not find directory ' + str(current_dir) + '\n')
      sys.exit(1)

  # create the output directory if need be
  if not os.path.isdir(alignment_dir + distance):
    call(['mkdir', alignment_dir + distance])

  ligand_to_group = ligand_groups()  # mapping of ligand ID -> super group

  # Process each domain separately:
  domain_ligand_pairs = set()  # keep track of which domain-ligand pairs we have started to process (no overwriting)

  # determine how many lines are in the domain_file (in order to print progress messages)
  domain_file_size = line_count_from_file(domain_file)
  progress_bars = [(str(rank*10)+'%', int(rank*(domain_file_size/10.))) for rank in range(1, 10)][::-1]

  domain_handle = gzip.open(domain_file) if domain_file.endswith('gz') else open(domain_file)

  process_index = 0
  for dom_line in domain_handle:
    if dom_line.startswith('#'):
      continue

    for progress_percent, progress_value in progress_bars:
      if process_index > progress_value:
        progress_bars = progress_bars[:-1]  # remove the last value to not reprint
        sys.stderr.write('Processed ' + progress_percent + ' (' + "{:,}".format(process_index) + '/' +
                         "{:,}".format(domain_file_size) + ') of the domain file entries.\n')
        break

    pdb_id_chain, domain_name, match_states = dom_line[:-1].split('\t')[:3]

    # (match state, 0-index pos, value)
    domain_positions = [tuple(entry.split(':')) for entry in match_states.split(',')]
    domain_start = int(min(map(int, [pos for (_, pos, _) in domain_positions])))
    domain_end = int(max(map(int, [pos for (_, pos, _) in domain_positions])))

    # location of the current FASTA file
    current_fasta_file = fasta_dir+pdb_id_chain[0]+'/'+pdb_id_chain[:2]+'/' + pdb_id_chain[:4]+'_'+distance+'.fa'
    if not os.path.isfile(current_fasta_file):
      continue

    # find the header that corresponds to our sequence of interest
    fasta_handle = open(current_fasta_file)
    for fasta_line in fasta_handle:
      if fasta_line.startswith('>' + pdb_id_chain):
        binding_positions = [entry.split('-') for entry in
                             fasta_line[fasta_line.find('bindingSiteRes=') + 15:fasta_line.rfind(';')].split(',')]

        # restrict to those binding positions that occur within the domain of interest with reasonable binding scores
        relevant_ligands = set([ligand_type for (pos, ligand_type, score) in binding_positions if
                                (domain_start <= int(pos) - 1 <= domain_end) and
                                ((distance != 'mindist' and float(score) > 0.) or
                                 (distance == 'mindist' and float(score) <= PROXIMITY_CUTOFF))])

        # for each ligand, write out a temporary alignment file (to be reformatted into FASTA format later)
        for ligand_type in relevant_ligands:

          # all super groups that ligand_type is part of:
          super_groups = ['ALL_', ligand_type] + ligand_to_group.get(ligand_type, [])

          # all entries that are not nucleic acids, ions, or peptides are considered small molecules, also
          if not ('NUCACID_' in super_groups or 'ION_' in super_groups or 'III' in super_groups):
            super_groups.append('SM_')

          # translate all names of the ligands (if need be):
          super_groups = [translate_ligand(orig_ligand_name) for orig_ligand_name in super_groups]

          for sub_ligand in super_groups:
            # name of file to write out the PDB ID and chain, domain start, domain end, and match states:
            tmp_aln_file = alignment_dir + distance + '/' + domain_name + '_' + sub_ligand + '_' + distance + '.tmpfa'

            # unique identifier for this particular domain
            sequence_id = pdb_id_chain + '_' + str(domain_start) + '_' + str(domain_end)

            if (domain_name, sub_ligand) not in domain_ligand_pairs:
              domain_ligand_pairs.add((domain_name, sub_ligand))
              aln_handle_out = open(tmp_aln_file, 'w')
            else:
              aln_handle_out = open(tmp_aln_file, 'a')  # do not overwrite previous values

            aln_handle_out.write(sequence_id + '\t' + match_states + '\n')
            aln_handle_out.close()
    fasta_handle.close()

    process_index += 1  # increment progress

  domain_handle.close()

  # reformat the temporary alignment files
  sys.stderr.write('Done!\n')
  progress_bars = [(str(rank * 10) + '%', int(rank * (len(domain_ligand_pairs) / 10.))) for rank in range(1, 10)][::-1]

  domain_name_list = set()
  ligand_type_list = set()
  successful_alignment_files = 0
  process_index = 0  # keep track of progress again
  for domain_name, ligand_type in domain_ligand_pairs:

    for progress_percent, progress_value in progress_bars:
      if process_index > progress_value:
        progress_bars = progress_bars[:-1]  # remove the last value to not reprint
        sys.stderr.write('FASTA-formatted ' + progress_percent + ' (' + "{:,}".format(process_index) + '/' +
                         "{:,}".format(len(domain_ligand_pairs)) + ') of the temporary alignment files.\n')
        break

    tmp_aln_file = alignment_dir + distance + '/' + domain_name + '_' + ligand_type + '_' + distance + '.tmpfa'

    # attempt the reformatting into a FASTA-formatted multiple alignment file
    if os.path.isfile(tmp_aln_file):
      new_aln_file = format_fasta_aln(tmp_aln_file,
                                      alignment_dir+distance+'/'+domain_name+'_'+ligand_type+'_'+distance+'.aln.fa')

      if new_aln_file:  # remove the temporary file if successful
        call(['rm', tmp_aln_file])
        domain_name_list.add(domain_name)
        ligand_type_list.add(ligand_type)
        successful_alignment_files += 1

    process_index += 1

  # print the success message
  sys.stderr.write('Done!\nSuccessfully created ' + "{:,}".format(successful_alignment_files) + ' alignment files ' +
                   'for ' + "{:,}".format(len(domain_name_list)) + ' domains and ' +
                   "{:,}".format(len(ligand_type_list)) + ' ligands.\n')
  sys.stderr.write('NEXT STEP to generate uniqueness weights: run\n' +
                   'python evaluate_uniqueness.py --distance '+distance+'\n')


########################################################################################################
# SEQUENCE UNIQUENESS SCORES
########################################################################################################

def henikoff_column_score(column):
  """
  :param column: ordered list of characters (corresponding to a column of a multiple alignment)
  :return: corresponding ordered list of scores according to the following:
           Henikoff S and Henikoff JG (1994). "Position-based sequence weights." J Mol Biol
           243(4): 574-578. doi:10.1016/0022-2836(94)90032-9
  """

  # unique characters appearing in column
  unique_chars = set(column)

  # number of times each unique character appears in the alignment column
  per_char_count = {curr_char: column.count(curr_char) for curr_char in unique_chars}

  # score corresponds to 1/(A*Bi), where A is the number of unique characters in the column, and
  # B is the number of times the character at position i appears in the column
  return [1./(len(unique_chars) * per_char_count[curr_char]) for curr_char in column]


########################################################################################################

def henikoff_alignment_score(seqid_to_sequence):
  """
  :param seqid_to_sequence: dictionary of sequence ID -> full sequence (where all sequences are the 
                            same length)
  :return: the raw (i.e., not yet normalized, doesn't necessarily sum to 1) uniqueness score for each
           sequence ID in the original alignment file
  """

  # remove the completely gapped columns (if any) from the score:
  new_seqid_to_sequence = remove_gapped_columns(seqid_to_sequence)

  # keep track of the per-column scores to generate overall per-sequence uniqueness scores
  seq_ids = sorted(new_seqid_to_sequence.keys())  # identifiers must be in the same order in each column
  sums_across_columns = [0.] * len(seq_ids)  # sum of per-column, per-sequence uniqueness scores

  for col_ind in xrange(len(new_seqid_to_sequence[seq_ids[0]])):
    current_column = [new_seqid_to_sequence[curr_seq_id][col_ind] for curr_seq_id in seq_ids]
    column_scores = henikoff_column_score(current_column)

    # update the running totals for each sequence's "uniqueness"
    for seq_ind in xrange(len(seq_ids)):
      sums_across_columns[seq_ind] += column_scores[seq_ind]

  return seq_ids, sums_across_columns


########################################################################################################

def remove_gapped_columns(seqid_to_sequence):
  """
  :param seqid_to_sequence: dictionary of sequence ID -> full sequence (where all sequences are the 
                            same length)
  :return: a dictionary of sequence ID -> sequence (where all completely gapped columns have been
           removed)
  """

  # check to see if there are any gapped columns
  gapped_columns = [index for index in xrange(len(seqid_to_sequence.values()[0]))
                    if set([sequence[index] for sequence in seqid_to_sequence.values()]) == {'-'}]

  if len(gapped_columns) < 1:
    return seqid_to_sequence

  # if there are, create new sequences excluding the fully-gapped columns
  seqid_to_ungapped_sequence = {sequence_id: [] for sequence_id in seqid_to_sequence.keys()}

  for index in xrange(len(seqid_to_sequence.values()[0])):
    if index in gapped_columns:
      continue

    for sequence_id, full_sequence in seqid_to_sequence.items():
      seqid_to_ungapped_sequence[sequence_id].append(full_sequence[index])

  # concatenate sequences and return
  return {sequence_id: ''.join(ungapped_seq) for sequence_id, ungapped_seq in seqid_to_ungapped_sequence.items()}


########################################################################################################

def find_closest_chain(pdb_id, pdb_chains, domain_location, ligand_type, distance):
  """
  :param pdb_id: string corresponding to the 4-digit PDB entry we are examining
  :param pdb_chains: the subset of structure chains (e.g., A, B, C...) that contain a domain instance
                     with an identical sequence in the same position
  :param domain_location: string of 0-index starting position '_' 0-index ending position
  :param ligand_type: string corresponding to the type of ligand we're calculating a distance to
  :return: one PDB chain that is the closest (of the subset provided) to the ligand of interest
  """

  ligand_to_group = ligand_groups()  # original ligand ID -> set of groups

  index_range = map(str, range(int(domain_location.split('_')[0])+1,
                               int(domain_location.split('_')[1])+2))  # all possible binding sites to consider

  distance_fasta = DATAPATH+'processed_data/fasta/'+pdb_id[0]+'/'+pdb_id[:2]+'/'+pdb_id+'_mindist.fa'
  if not os.path.isfile(distance_fasta):
    sys.stderr.write('Could not open '+distance_fasta+'\n')
    return sorted(pdb_chains)[0]

  chain_proximity = {}
  distance_handle = open(distance_fasta)
  for dist_line in distance_handle:
    if dist_line.startswith('>'+pdb_id) and dist_line[len('>'+pdb_id):len('>'+pdb_id)+1] in pdb_chains:

      current_chain = dist_line[len('>'+pdb_id):len('>'+pdb_id)+1]
      if current_chain not in chain_proximity:
        chain_proximity[current_chain] = {}  # aa position -> distance to ligand (within 5 angstroms)

      curr_bind_pos = dist_line[dist_line.find('bindingSiteRes=') + 15:dist_line.rfind(';')].split(',')
      curr_bind_pos = [entry.split('-') for entry in curr_bind_pos]

      for (aapos, current_ligand, binding_score) in curr_bind_pos:

        super_groups = ['ALL_', current_ligand] + ligand_to_group.get(current_ligand, [])
        if not ('NUCACID_' in super_groups or 'ION_' in super_groups or 'III' in super_groups):
          super_groups.append('SM_')

        # translate all names of the ligands (if need be):
        super_groups = [translate_ligand(orig_ligand_name) for orig_ligand_name in super_groups]

        # if we are looking at the correct ligand type and this position is within the domain range we want:
        if ligand_type in super_groups and aapos in index_range:
          if (distance == 'mindist' and float(binding_score) <= PROXIMITY_CUTOFF) or \
              (distance != 'mindist' and float(binding_score) >= PROXIMITY_CUTOFF):
            if aapos not in chain_proximity[current_chain]:
              chain_proximity[current_chain][aapos] = float(binding_score)

            if distance == 'mindist':
              chain_proximity[current_chain][aapos] = min(chain_proximity[current_chain][aapos],
                                                          float(binding_score))
            else:
              chain_proximity[current_chain][aapos] = max(chain_proximity[current_chain][aapos],
                                                          float(binding_score))
  distance_handle.close()

  total_proximity = []
  all_chains = sorted(list(chain_proximity.keys()))
  for chain_id, pos_to_dist in chain_proximity.items():
    total_proximity.append((len(pos_to_dist.keys()),  # total binding positions
                            (-1 if distance == 'mindist' else 1)*sum(pos_to_dist.values()),  # maximize this value
                            -1*all_chains.index(chain_id)))  # inverse of the chain ID ('A', 'B', 'C') -> (0, -1, -2)

  return all_chains[-1*sorted(total_proximity, reverse=True)[0][2]]


########################################################################################################

def clear_crystal_duplicates(alignment_file, ligand_type, distance):
  """
  :param alignment_file: full path to a fasta file with sequence ID (PDB ID, PDB chain, start index, end index)
  :return: a single PDB structure can have multiple (identical) protein chains, even if they belong
           to biological assemblies (there can be many different ones specified for the same structure).
           For each domain--ligand pair, we only consider sets of unique (PDB ID [no chain], ligand type, sequence)
           before evaluating uniqueness; we select the identical chain that is closest to the ligand.
  """

  unique_domhits = {}
  seqid_to_sequence = {}  # sequence ID (PDB ID, PDB Chain, start index, end index) -> string of complete sequence
  current_seq_id = ''  # keep track of the current sequence ID (while parsing the multiple alignment file)

  fasta_handle = gzip.open(alignment_file) if alignment_file.endswith('gz') else open(alignment_file)
  for aln_line in fasta_handle:
    if aln_line.startswith('>'):
      current_seq_id = aln_line[1:-1].split()[0]  # remove new line and starting '>'
      pdb_id = current_seq_id[:4]
      pdb_chain = current_seq_id[4]
      dom_loc = current_seq_id[6:]

      if pdb_id not in unique_domhits:
        unique_domhits[pdb_id] = {}
      if dom_loc not in unique_domhits[pdb_id]:
        unique_domhits[pdb_id][dom_loc] = set()
      unique_domhits[pdb_id][dom_loc].add(pdb_chain)

      seqid_to_sequence[current_seq_id] = []  # start to keep track of current sequence
    else:
      seqid_to_sequence[current_seq_id].append(aln_line.strip())  # sequence without newlines
  fasta_handle.close()

  for sequence_id in seqid_to_sequence.keys():
    seqid_to_sequence[sequence_id] = ''.join(seqid_to_sequence[sequence_id])

  # go through the sequences to store unique domains from chains:
  # pdbID -> (sequence, loc) -> [chains]  # identical chains will have had identical domain locations...
  unique_seq_ids = []
  for pdb_id in unique_domhits.keys():
    for dom_loc, chains in unique_domhits[pdb_id].items():

      if len(chains) > 1:

        # ONLY in this case do we bother checking the sequence itself:
        all_seqs = {}
        for chain in chains:
          current_seq = seqid_to_sequence[pdb_id+chain+'_'+dom_loc]
          if current_seq not in all_seqs:
            all_seqs[current_seq] = set()
          all_seqs[current_seq].add(chain)

        for current_seq, new_chains in all_seqs.items():
          # now, we should pick the closest:
          if len(new_chains) > 1:
            # we want to pick the CLOSEST chain when there are multiple chains in contact with the ligand
            unique_seq_ids.append(pdb_id + find_closest_chain(pdb_id, new_chains, dom_loc,
                                                              ligand_type, distance) + '_' + dom_loc)

          # otherwise, add each "unique" domain (i.e., unique sequence, but same location):
          else:
            unique_seq_ids.append(pdb_id + sorted(list(new_chains))[0] + '_' + dom_loc)

      else:
        unique_seq_ids.append(pdb_id + sorted(list(chains))[0] + '_' + dom_loc)

  return {seq_id: seq for seq_id, seq in seqid_to_sequence.items() if seq_id in unique_seq_ids}, \
         set([seq_id for seq_id in seqid_to_sequence.keys() if seq_id not in unique_seq_ids])


########################################################################################################

def overall_alignment_score(alignment_file, ligand_type, distance):
  """
  :param alignment_file: full path to a FASTA-formatted multiple alignment file
  :return: the raw (i.e., not yet normalized, doesn't necessarily sum to 1) uniqueness score for each
           sequence ID in the original alignment file
  """

  # read in the full sequence for each sequence ID
  seqid_to_sequence, missing_seqids = clear_crystal_duplicates(alignment_file, ligand_type, distance)

  # make sure that the sequences in this alignment are all the same length
  assert len(set([len(seq) for seq in seqid_to_sequence.values()])) == 1, "Varying sequence lengths in "+alignment_file

  return henikoff_alignment_score(seqid_to_sequence)


########################################################################################################

def normalize_scores(scores):
  """
  :param scores: list, set, array of positive numerical (i.e., int / float) values
  :return: the same list, set, array of values such that their sum is 1
  """

  current_total = sum(map(float, scores))

  return [float(current_value)/max(current_total, 1.) for current_value in scores]


########################################################################################################

def generate_uniqueness_scores(domain_names, ligand_types, alignment_files, output_file):
  """
  :param domain_names: list of domain names corresponding to each alignment file in alignment_files
  :param ligand_types: list of ligand binding types corresponding to each alignment file in alignment_files
  :param alignment_files: list of full paths to FASTA-formatted alignment files to calculate uniqueness scores for
  :param output_file: full path to an output file to write tab-delimited results. Each column will correspond to:
                      [0] domain name, [1] ligand binding type, [2] comma-delimited list of
                      pdbID_domain-start_domain-end: relative uniqueness
  :return: None, but print a success message
  """

  # make sure that we can actually open the file to write to; if so, do so
  if not os.path.isdir('/'.join(output_file.split('/')[:-1])):
    sys.stderr.write('Could not open '+output_file+' to write to. Exiting.\n')
    sys.exit(1)

  output_handle = gzip.open(output_file, 'w') if output_file.endswith('gz') else open(output_file, 'w')
  output_handle.write('\n'.join(['# Per domain sequence "uniqueness" scores for ' +
                                 "{:,}".format(len(set(domain_names))) + ' distinct domains contacting ' +
                                 "{:,}".format(len(set(ligand_types))) + ' distinct ligands, calculated as in',
                                 '# Henikoff S and Henikoff JG (1994). "Position-based sequence weights." ' +
                                 'J Mol Biol 243(4): 574-578. doi:10.1016/0022-2836(94)90032-9'])+'\n')
  output_handle.write('\t'.join(['#domain_name', 'ligand_binding_type', 'number_of_instances',
                                 'pdbID-pdbChain_domain-start_domain-end : relative_uniqueness, ...'])+'\n')

  progress_bars = [(str(rank * 10) + '%', int(rank * (len(alignment_files) / 10.))) for rank in range(1, 10)][::-1]

  all_uniqueness_scores = []

  # calculate scores for each alignment file
  for curr_ind, curr_aln_file in enumerate(alignment_files):

    for progress_percent, progress_value in progress_bars:
      if curr_ind > progress_value:
        progress_bars = progress_bars[:-1]  # remove the last value to not reprint
        sys.stderr.write('Calculated uniqueness scores for '+progress_percent+' ('+"{:,}".format(curr_ind)+'/' +
                         "{:,}".format(len(alignment_files)) + ') of the MSA files.\n')
        break

    distance = curr_aln_file[curr_aln_file.rfind('_')+1:curr_aln_file.rfind('.aln.fa')]
    trunc_aln_file = curr_aln_file.split('/')[-1].replace('_'+distance+'.aln.fa', '')
    ligand_type = trunc_aln_file[trunc_aln_file[:-1].rfind('_') + 1:]

    # PDB ID -> unnormalized uniqueness score
    ordered_seqids, seqid_to_score = overall_alignment_score(curr_aln_file, ligand_type, distance)
    relative_scores = normalize_scores(seqid_to_score)  # relative scores

    # write the new line out to file (domain name, ligand type, number of instances, sequence_id:relative_weight,...)
    all_uniqueness_scores.append('\t'.join([str(domain_names[curr_ind]), str(ligand_types[curr_ind]),
                                            str(len(ordered_seqids)),
                                            ','.join(str(ordered_seqids[curr_id])+':'+str(relative_scores[curr_id])
                                                     for curr_id in xrange(len(ordered_seqids)))])+'\n')

  all_uniqueness_scores.sort()
  for unique_line in all_uniqueness_scores:
    output_handle.write(unique_line)
  output_handle.close()
  sys.stderr.write('Done!\nSuccessfully wrote relative uniqueness scores to '+output_file+'\n')


########################################################################################################

if __name__ == "__main__":

  # Parse the command-line arguments
  parser = argparse.ArgumentParser(description='Generate multiple sequence alignments of each domain-' +
                                               'ligand pair to calculate per-sequence "uniqueness" scores.')

  parser.add_argument('--distance', type=str,
                      help='How to record the distance between receptor and ligand?',
                      default='mindist',
                      choices={'fracin4', 'mindist', 'meandist', 'maxstd', 'meanstd', 'sumstd',
                               'maxvdw', 'meanvdw', 'sumstd'})
  parser.add_argument('--create_alignments', dest='create_alignments', action='store_true',
                      help='Use the file containing all domain hits to create per-domain, per-ligand alignments?')
  parser.set_defaults(create_alignments=False)

  args = parser.parse_args()

  # ----------------------------------------------------------------------------------------------------
  if args.create_alignments:
    """
    For each domain, generate alignments of sequences for each type of ligand the domain contacts, to 
    downweight redundant sequences that are prevalent in structural data.
    """

    sys.stderr.write('Creating per-domain, per ligand type, multiple sequence alignments.\n')

    create_alignment_files(DOMAINS,
                           DATAPATH+'processed_data/fasta/',
                           DATAPATH+'processed_data/domains/alignments/',
                           args.distance)

  # ----------------------------------------------------------------------------------------------------
  else:
    """
    Once alignments have been made, use the alignments to generate a single "sequence uniqueness" file to 
    be used in the step to create binding potential scores
    """

    # check input directory existence:
    for subdir in ['', 'domains', 'domains/alignments', 'domains/alignments/'+args.distance]:
      if not os.path.isdir(DATAPATH+'processed_data/'+subdir):
        sys.stderr.write('ERROR, no such directory: '+DATAPATH+'processed_data/'+subdir+'\n')
        sys.stderr.write('Please run python '+os.getcwd()+'/evaluate_uniqueness.py --create_alignments --' +
                         args.distance+'\n')
        sys.exit(1)

    # get the list of processed alignment files
    current_domain_names = []
    current_ligand_types = []
    current_align_files = []

    # file format is [domain name]_[ligand type]_[distance].aln.fa
    for aln_file in os.listdir(DATAPATH+'processed_data/domains/alignments/'+args.distance):
      if not aln_file.endswith('_'+args.distance+'.aln.fa'):
        continue

      current_align_files.append(DATAPATH+'processed_data/domains/alignments/'+args.distance+'/'+aln_file)
      filename = aln_file.split('/')[-1].replace('_'+args.distance+'.aln.fa', '')

      # note, the ligand type may have a single *ending* underscore, and the domain name may include underscores
      curr_ltype = filename[filename[:-1].rfind('_')+1:]

      current_ligand_types.append(curr_ltype)
      current_domain_names.append(filename[:filename.rfind('_'+curr_ltype)])

    generate_uniqueness_scores(current_domain_names, current_ligand_types, current_align_files,
                               DATAPATH+'processed_data/domains/uniqueness-scores_'+args.distance+'.txt.gz')
