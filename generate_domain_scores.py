#!/usr/bin/python

"""
Using the previously computed per-sequence "uniqueness" scores (evaluate_uniqueness.py), along with the 
per-PDB-sequence binding information from the previously computed FASTA files (create_fasta.py), calculate
per domain-position binding potential weights with respect to each type of ligand that domain can bind to.

Contact Shilpa Nadimpalli Kobren (snadimpa@princeton.edu) with questions
"""

import os
import sys
import gzip
import argparse
from subprocess import call
from evaluate_uniqueness import ligand_groups, translate_ligand, normalize_scores


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
# DISTRIBUTION SUMMARY FUNCTIONS
########################################################################################################

def check_and_renormalize_distribution(distribution):
  """
  :param distribution: list of (value, weight) tuples such that the weights sum to 1
  :return: a renormalized distribution such that weights sum to 1
  """

  if not abs(1. - sum([weight for (_, weight) in distribution])) < 1e-6:
    normalized_weights = normalize_scores([weight for (_, weight) in distribution])
    return [(distribution[i][0], normalized_weights[i]) for i in xrange(len(distribution))]

  return distribution  # if the original weights were fine, return the original distribution


########################################################################################################

def weighted_sum(distribution):
  """
  :param distribution: list of (value, weight) tuples such that the weights sum to 1
  :return: the weighted sum of the distribution
  """

  fixed_dist = check_and_renormalize_distribution(distribution)

  return sum([weight * value for (value, weight) in fixed_dist])


########################################################################################################

def weighted_median(distribution):
  """
  :param distribution: list of (value, weight) tuples such that the weights sum to 1
  :return: the weighted median of the distribution
  """

  fixed_dist = check_and_renormalize_distribution(distribution)

  sorted_scores = sorted(fixed_dist)  # sorted by values first
  total_weight = 0.  # keep track of the total weight seen so far; the 0.5 mark will be the center

  for (value, weight) in sorted_scores:
    total_weight += weight
    if total_weight >= 0.5:
      return value  # value in the center

  # if we didn't exit the loop properly, simply return the unweighted median
  return sorted_scores[len(fixed_dist) / 2][0]


########################################################################################################

def weighted_fraction(distribution, max_cutoff=PROXIMITY_CUTOFF):
  """
  :param distribution: list of (value, weight) tuples such that the weights sum to 1
  :param max_cutoff: the maximum cutoff to consider an element in the distribution
  :return: the weighted fraction of the distribution that was <= the specified maximum cutoff
  """

  fixed_dist = check_and_renormalize_distribution(distribution)

  return sum([weight for (value, weight) in fixed_dist if value <= max_cutoff])


########################################################################################################

def choose_summary_function(distance):
  """
  :param distance: string corresponding to the distance metric to be used
  :return: a function that takes in a list of pair tuples with (value, relative weight) and returns a 
           single score AND the name of the function
  """

  # set the default way by which to "flatten" per-position score distributions into a single positional score
  if distance == 'mindist':
    return weighted_fraction, 'weighted_fraction_within_'+str(PROXIMITY_CUTOFF)+'A'

  elif distance in ['sumstd', 'sumvdw']:
    # the "sum" values range from 0 -> infinity, so we store the median
    return weighted_median, 'weighted_median'

  else:
    # the "mean", "max", and "fracin4" values range from 0->1, and, as the relative uniqueness scores sum to 1,
    # taking the weighted sum makes sense (resulting positional score must also range from 0 to 1)
    return weighted_sum, 'weighted_sum'


########################################################################################################
# PROCESS REQUIRED INPUT FILES
########################################################################################################

def process_uniqueness_file(uniqueness_file):
  """
  :param uniqueness_file: full path to a tab-delimited file with columns domain name, ligand type, and
                          comma-separated list of sequence ID: relative uniqueness weight
  :return: dictionary of domain_name -> ligand_type -> pddID-pdbChain_start_end -> relative_uniqueness
  """

  uniqueness = {}  # domain -> ligand type -> pdbID_start_end -> relative_uniqueness

  uniqueness_handle = gzip.open(uniqueness_file) if uniqueness_file.endswith('gz') else open(uniqueness_file)
  for weight_line in uniqueness_handle:
    if weight_line.startswith('#'):
      continue

    domain_name, ligand_type, number_instances, seq_to_uniqueness = weight_line[:-1].split('\t')[:4]

    if domain_name not in uniqueness:
      uniqueness[domain_name] = {}

    if ligand_type not in uniqueness[domain_name]:
      uniqueness[domain_name][ligand_type] = {}

    # format of these entries is 'pdbID_start_end:relative_uniqueness,...'
    seq_to_uniqueness = [entry.split(':') for entry in seq_to_uniqueness.split(',')]
    uniqueness[domain_name][ligand_type] = {seqid: float(rel_wt) for (seqid, rel_wt) in seq_to_uniqueness}
  uniqueness_handle.close()

  return uniqueness


########################################################################################################

def process_fasta_file(fasta_file, pdbid_pdbchain_subset, ligand_to_group, binding_positions):
  """
  :param fasta_file: full path to a FASTA-formatted file as generated by create_fasta.py
  :param pdbid_pdbchain_subset: subset of pdbID-pdbChain identifiers to consider
  :param ligand_to_group: get the mapping from ligand type to super group
  :param binding_positions: existing dictionary to be updated (if need be) of ligand_type -> pdbID-pdbChain ->
                            1-index AA position -> binding_score
  :return: None, but update the input binding_positions data structure
  """

  fasta_handle = gzip.open(fasta_file) if fasta_file.endswith('gz') else open(fasta_file)
  for fasta_line in fasta_handle:
    if fasta_line.startswith('>'):
      # format of these entries is pdbID-pdbChain bindingSiteRes=1-index AA position : ligand_type : score,...;
      pdbid_pdbchain = fasta_line[1:-1].split()[0]

      # only consider those sequences that had a domain match
      if pdbid_pdbchain not in pdbid_pdbchain_subset:
        continue

      curr_bind_pos = fasta_line[fasta_line.find('bindingSiteRes=') + 15:fasta_line.rfind(';')].split(',')
      curr_bind_pos = [entry.split('-') for entry in curr_bind_pos]

      for (aapos, current_ligand, binding_score) in curr_bind_pos:

        super_groups = ['ALL_', current_ligand] + ligand_to_group.get(current_ligand, [])
        if not ('NUCACID_' in super_groups or 'ION_' in super_groups or 'III' in super_groups):
          super_groups.append('SM_')

        # translate all names of the ligands (if need be):
        super_groups = [translate_ligand(orig_ligand_name) for orig_ligand_name in super_groups]

        for ligand_type in super_groups:

          if ligand_type not in binding_positions:
            binding_positions[ligand_type] = {}
          if pdbid_pdbchain not in binding_positions[ligand_type]:
            binding_positions[ligand_type][pdbid_pdbchain] = {}

          binding_positions[ligand_type][pdbid_pdbchain][aapos] = float(binding_score)

  fasta_handle.close()


########################################################################################################

def process_alignment_file(align_file, uniqueness_weights, lig_binding_positions, default_score):
  """
  :param align_file: full path to an alignment file, where we are interested in headers of the format 
                     'pdbID-pdbChain_start_end \t matchstate : 1-index AA position,...'
  :param uniqueness_weights: relative (i.e., sum to 1) weights for each pdbID-pdbChain_start_end domain identifier
  :param lig_binding_positions: dictionary of pdbid_pdbchain -> binding potential score
  :param default_score: if a binding potential score is unavailable, what is the default value? 
  :return: dictionary of match state -> [(value, weight), ]
  """

  local_mstate_values = {}  # match state -> pdbID-pdbChain -> (value, weight)

  aln_handle = gzip.open(align_file) if align_file.endswith('gz') else open(align_file)
  for aln_line in aln_handle:
    if aln_line.startswith('>'):

      # format of these entries is 'pdbID-pdbChain_start_end \t matchstate : 1-index AA position,...'
      seqid = aln_line[1:-1].split()[0]

      if seqid not in uniqueness_weights:  # do we have a uniqueness score for this sequence ID?
        continue
      relative_uniqueness = uniqueness_weights[seqid]

      pdbid_pdbchain = seqid.split('_')[0]  # pdbID-pdbChain (needed to index into binding positions)

      # store the binding potential weights for each match state in the domain
      mstates = [entry.split(':') for entry in aln_line[:-1].split('\t')[1].split(',')]
      for (matchstate, aapos) in mstates:

        if matchstate not in local_mstate_values:
          local_mstate_values[matchstate] = {}

        if pdbid_pdbchain not in lig_binding_positions or aapos not in lig_binding_positions[pdbid_pdbchain]:
          current_binding_weight = default_score
        else:
          current_binding_weight = lig_binding_positions[pdbid_pdbchain][aapos]

        local_mstate_values[matchstate][seqid] = (current_binding_weight, relative_uniqueness)
  aln_handle.close()

  return local_mstate_values


########################################################################################################
# POSITIONAL BINDING SCORES
########################################################################################################

def create_binding_scores(uniqueness_file, fasta_dir, alignment_dir, binding_score_dir, distance):
  """
  :param uniqueness_file: full path to a tab-delimited file with columns domain name, ligand type, and
                          comma-separated list of sequence ID: relative uniqueness weight
  :param fasta_dir: full path to a directory containing per-PDB-ID sequences and their positional scores
  :param alignment_dir: full path to a directory containing precomputed (domain, ligand) FASTA alignments
  :param binding_score_dir: full path to a directory to store output files
  :param distance: type of scoring metric to be used (e.g., mindist, fracin4, etc.) -- needed for naming
  :return: None, but print success message for the number of binding scores processed
  """

  # set the default way by which to "flatten" per-position score distributions into a single positional score
  summarize_position_func, column_name = choose_summary_function(distance)

  # all domain--ligand pairs to process are found in the uniqueness file:
  # domain -> ligand type -> pdbID_start_end -> relative_uniqueness
  uniqueness = process_uniqueness_file(uniqueness_file)

  total_processed_domains = 0
  progress_bars = [(str(rank * 10) + '%', int(rank * (len(uniqueness.keys()) / 10.))) for rank in range(1, 10)][::-1]

  # get the corresponding binding scores for all PDB files that matched each domain
  for domain_name in uniqueness.keys():

    for progress_percent, progress_value in progress_bars:
      if total_processed_domains > progress_value:
        progress_bars = progress_bars[:-1]  # remove the last value to not reprint
        sys.stderr.write('Processed ' + progress_percent + ' (' + "{:,}".format(total_processed_domains) + '/' +
                         "{:,}".format(len(uniqueness.keys())) + ') of interaction domains.\n')
        break

    # all unique pdbID-pdbChain identifiers with 1+ domain matches
    all_matching_pdbids = set([seqid.split('_')[0] for seqdict in uniqueness[domain_name].values()
                               for seqid in seqdict.keys()])

    # ligand_type -> pdbID-pdbChain -> 1-index AA position -> binding score
    binding_positions = {}
    ligand_mapping = ligand_groups()

    for pdbid in [pdbid_pdbchain[:4] for pdbid_pdbchain in all_matching_pdbids]:
      fasta_file = fasta_dir + pdbid[0] + '/' + pdbid[:2] + '/' + pdbid + '_' + distance + '.fa'

      if not os.path.isfile(fasta_file):
        sys.stderr.write('No such file: '+fasta_file+'\n')
        continue

      # update the binding_positions dictionary
      process_fasta_file(fasta_file, all_matching_pdbids, ligand_mapping, binding_positions)

    # store per-position, per ligand-binding type binding scores to an outfile:
    out_file = binding_score_dir + domain_name + '_binding-scores_' + distance + '.txt.gz'
    out_handle = gzip.open(out_file, 'w') if out_file.endswith('gz') else open(out_file, 'w')
    out_handle.write('# Continuous positional weights, calculated according to the '+distance+' statistic,' +
                     ' for '+domain_name+'\n')
    out_handle.write('\t'.join(['#ligand_type', 'match_state', column_name,
                                'distribution (pdbID-pdbChain_start_end : relative uniqueness weight : ' +
                                'positional score),...'])+'\n')

    # for each type of ligand-binding in this domain, obtain overall binding score distributions for each match state
    for ligand_type, seqid_to_uniqueness in uniqueness[domain_name].items():

      # read in the match state information:
      aln_file = alignment_dir + domain_name + '_' + ligand_type + '_' + distance + '.aln.fa'

      if not os.path.isfile(aln_file):
        sys.stderr.write('No such file: '+aln_file+'\n')
        continue

      # match state -> pdbID-pdbChain_start_end -> (binding score/value, relative uniqueness/weight)
      match_state_distributions = process_alignment_file(aln_file, seqid_to_uniqueness,
                                                         binding_positions[ligand_type],
                                                         DISTANCE_CUTOFF if distance == 'mindist' else 0.)

      # print the scores to file:
      for matchstate, distribution in match_state_distributions.items():

        flattened_score = summarize_position_func(distribution.values())

        if flattened_score > 0.:

          # write out the ligand type, match state, positional score, and complete distribution (worth recording)
          out_handle.write('\t'.join([ligand_type, str(matchstate), str(flattened_score),
                                      ','.join([seqid+':'+str(rel_wt)+':'+str(value) for
                                                seqid, (value, rel_wt) in sorted(distribution.items())
                                                if rel_wt > 0. and
                                                ((distance == 'mindist' and value < DISTANCE_CUTOFF) or
                                                 (distance != 'mindist' and value > 0.))])])+'\n')
    out_handle.close()
    total_processed_domains += 1

  sys.stderr.write('Successfully wrote positional scores for '+"{:,}".format(total_processed_domains)+' domains to ' +
                   binding_score_dir+'\n')


########################################################################################################

if __name__ == "__main__":

  # Parse the command-line arguments
  parser = argparse.ArgumentParser(description='Generate per-position binding potential scores for each ' +
                                               'domain with respect to each type of ligand it can bind to.')

  parser.add_argument('--distance', type=str,
                      help='How to record the distance between receptor and ligand?',
                      default='mindist',
                      choices={'fracin4', 'mindist', 'meandist', 'maxstd', 'meanstd', 'sumstd',
                               'maxvdw', 'meanvdw', 'sumstd'})

  args = parser.parse_args()

  # Make sure that all input files are present:
  uniqueness_scores_file = DATAPATH+'processed_data/domains/uniqueness-scores_'+args.distance+'.txt.gz'
  fasta_file_directory = DATAPATH+'processed_data/fasta/'  # subdirectory example: 2/2m/
  alignments_directory = DATAPATH+'processed_data/domains/alignments/'+args.distance+'/'

  if not os.path.isfile(uniqueness_scores_file):
    sys.stderr.write('Could not read uniqueness scores from '+uniqueness_scores_file+'\'n')
    sys.stderr.write('Please run python evaluate_uniqueness.py --distance '+args.distance+'\n')
    sys.exit(1)

  if not os.path.isdir(fasta_file_directory):
    sys.stderr.write('Could not find FASTA files in '+fasta_file_directory+'\n')
    sys.stderr.write('Please run python create_fasta.py --distance '+args.distance+'\n')
    sys.exit(1)

  if not os.path.isdir(alignments_directory):
    sys.stderr.write('Could not read per-domain alignments from '+alignments_directory+'\n')
    sys.stderr.write('Please run python evaluate_uniqueness.py --create_alignments --distance '+args.distance+'\n')
    sys.exit(1)

  # Create the new binding scores output directory if needed
  output_directory = DATAPATH+'processed_data/domains/binding_scores/'

  for subdir in ['', args.distance]:
    if not os.path.isdir(output_directory+subdir):
      call(['mkdir', output_directory+subdir])

  # Generate the binding scores:
  create_binding_scores(uniqueness_scores_file, fasta_file_directory, alignments_directory,
                        output_directory+args.distance+'/', args.distance)
