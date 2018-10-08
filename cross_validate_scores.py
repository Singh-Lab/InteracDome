#!/usr/bin/python

"""
Generate cross-validate the accuracy and consistency of domain-based binding propensities.
(1) compute precision, recall, true positive rate, false positive rate at varying binding
    propensities (10-fold CV) without any grouping of domain instances, when grouping
    identical instances, and when grouping instances with >90% sequence identity.
(2) randomly split all domain-ligand instances and compare the residue-to-ligand distances
    using PCC (to measure consistency)
(3) use empirical bootstrapping to estimate the standard error of each binding propensity

Contact Shilpa Nadimpalli Kobren (snadimpa@princeton.edu) with questions
"""

import os
import sys
import gzip
import argparse
import numpy as np
import scipy.stats as ss
from subprocess import call
from networkx import Graph, connected_components
from evaluate_uniqueness import henikoff_alignment_score
from generate_domain_scores import PROXIMITY_CUTOFF, choose_summary_function


########################################################################################################
# CONSTANTS
########################################################################################################

# path to where this script is currently located (and to where all data should be stored) -- this can
# be updated
DATAPATH = os.path.dirname(os.path.abspath(__file__))+'/'

# full path to where alignments (created using evaluate_uniqueness.py --create_alignments) can be found
ALN_PATH = DATAPATH+'processed_data/domains/alignments/'

# full path to all per-domain binding weights (created using generate_binding_scores.py)
SCORE_PATH = DATAPATH+'processed_data/domains/binding_scores/'

# full path to where output files containing cross-validation evaluations will be written
CV_PATH = DATAPATH+'processed_data/domains/cross_validation/'


########################################################################################################
# PROCESS INPUT MATCH STATE / SEQUENCE DATA
########################################################################################################

def process_domain_alignment(alignment_file):
  """
  :param alignment_file: full path to a FASTA-formatted alignment file
  :return: dictionary of sequence ID -> sequence (with gaps), where domains from identical crystal
           chains (from the same PDB ID) have already been excluded..
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
      unique_seq_ids.append(pdb_id + sorted(list(chains))[0] + '_' + dom_loc)

  return {seq_id: seq for seq_id, seq in seqid_to_sequence.items() if seq_id in unique_seq_ids}


########################################################################################################

def ordered_domain_matchstates(alignment_file):
  """
  :param alignment_file: full path to a FASTA-formatted alignment file
  :return: ordered list of all match states in the domain (including insertions, if the sequence identifers
           we are interested in contain binding positions in those sites)
  """

  if not os.path.isfile(alignment_file):
    sys.stderr.write('No alignment file ' + alignment_file + '\n')
    return []

  ordered_matchstates = []  # ordered list of strings corresponding to domain match states

  aln_handle = gzip.open(alignment_file) if alignment_file.endswith('gz') else open(alignment_file)
  for fasta_line in aln_handle:
    if fasta_line.startswith('>'):
      matchstate_to_index = fasta_line[1:-1].split()[1]  # e.g., 1:3,2:4,3:5,4:6,a4-0:7,...
      ordered_matchstates = [entry.split(':')[0] for entry in matchstate_to_index.split(',')]
      break
  aln_handle.close()

  return ordered_matchstates


########################################################################################################

def process_domain_distances(score_file, ligand_type, default_value):
  """
  :param score_file: full path to a tab-delimited file with columns ligand type, match state,
                     summary score statistic, comma-delimited list of sequence identifier, relative
                     uniqueness, and value
  :param ligand_type: type of ligand to be returning scores for
  :param default_value: default value to be used if no value specified
  :return: dictionary of match state -> sequence identifier -> "distance" value (only non-default values)
  """

  if not os.path.isfile(score_file):
    sys.stderr.write('No binding score file '+score_file+'\n')
    return {}

  matchstate_to_value = {}  # match state -> sequence ID -> "distance" value

  score_handle = gzip.open(score_file) if score_file.endswith('gz') else open(score_file)
  for score_line in score_handle:
    if score_line.startswith('#'):
      continue

    current_ligand_type, match_state, _, seqid_to_value = score_line[:-1].split('\t')[:4]

    # only grab scores for the ligand type we are interested in
    if current_ligand_type != ligand_type:
      continue

    for seq_id, _, value in [entry.split(':')[:3] for entry in seqid_to_value.split(',')]:
      if False and float(value) >= default_value:  # if the distance was always >= 20, for instance...
        continue
      if match_state not in matchstate_to_value:
        matchstate_to_value[match_state] = {}

      matchstate_to_value[match_state][seq_id] = value
  score_handle.close()

  return matchstate_to_value


########################################################################################################

def domain_ligand_types(score_file):
  """
  :param score_file: full path to a tab-delimited file containing positional scores for different
                     ligand types
  :return: a set of ligand types found in the input domain binding score file
  """

  ligand_types = set()

  score_handle = gzip.open(score_file) if score_file.endswith('gz') else open(score_file)
  for score_line in score_handle:
    if score_line.startswith('#'):
      continue
    ligand_types.add(score_line[:-1].split('\t')[0])
  score_handle.close()

  return ligand_types


########################################################################################################

def domain_distance_vectors(score_file, alignment_file, ligand_type, default_value=20.):
  """
  :param score_file: full path to a tab-delimited file containing positional scores for different
                     ligand types
  :param alignment_file: full path to a FASTA-formatted alignment file
  :param ligand_type: type of ligand to restrict results to (distance/values)
  :param default_value: 20. for minimum/average distance, 0. for all other scores
  :return: dictionary of sequence ID -> ordered tuple of "distance" values per match state 
  """

  # match state -> sequence ID -> value (only non-default values are returned by this function)
  distances_to_ligand = process_domain_distances(score_file, ligand_type, default_value)

  # all possible ordered match states (remove insertion states if we have no non-default values there)
  ordered_match_states = []
  for match_state in ordered_domain_matchstates(alignment_file):
    try:  # include all standard match state values
      ordered_match_states.append(str(int(match_state)))
    except ValueError:  # invalid literal for int() with base 10
      if match_state in distances_to_ligand:  # if this insertion state was involved in binding?
        ordered_match_states.append(match_state)

  # all sequence identifiers with non-default values and 1+ "binding residues"
  all_sequence_ids = set(unravel_list([seqid_to_value.keys() for seqid_to_value in distances_to_ligand.values()]))

  distance_vectors = {}  # sequence ID -> (values per position)
  for sequence_id in all_sequence_ids:

    current_values = []
    for match_state in ordered_match_states:
      current_values.append(float(distances_to_ligand.get(match_state, {}).get(sequence_id, default_value)))

    # store the tuple of values for this domain instance IFF there is at least one position within PROXIMITY_CUTOFF
    if True in [dist <= PROXIMITY_CUTOFF for dist in current_values]:
      distance_vectors[sequence_id] = tuple(current_values)

  return distance_vectors


########################################################################################################
# SPLIT DATA INTO INDEPENDENT FOLDS
########################################################################################################

def percent_identity(seq1, seq2):
  """
  :param seq1: string corresponding to a sequence
  :param seq2: string corresponding to a sequence *of the same length* as seq1
  :return: the percent of the two sequences that matched exactly (ranges from 0 to 1)
  """

  if len(seq1) != len(seq2):
    sys.stderr.write('Length of seq1 != length of seq2!\n')
    return 0.

  return sum([1. if seq1[i] == seq2[i] else 0. for i in xrange(len(seq1))]) / len(seq1)


########################################################################################################

def create_shared_sequence_groups(seqid_to_sequence=None, sequence_identity_cutoff=0.9):
  """
  :param seqid_to_sequence: dictionary of sequence ID -> sequence, where all sequences are the *same length*
  :param sequence_identity_cutoff: percent identity cutoff to assign sequences to the same "shared sequence" 
                                   group
  :return: list of lists, where each sub list contains a subset of the sequence IDs from seqid_to_sequence
  """

  # get the full list of sequence identifiers:
  sequence_ids = sorted(list(seqid_to_sequence.keys())) if seqid_to_sequence else []

  # assign "edges" in a graph between sequence IDs (nodes) with >=sequence_identity_cutoff identity
  edges = set()

  for index, seqid_1 in enumerate(sequence_ids[:-1]):
    for seqid_2 in sequence_ids[index+1:]:
      if percent_identity(seqid_to_sequence[seqid_1],
                          seqid_to_sequence[seqid_2]) >= sequence_identity_cutoff:
        edges.add((seqid_1, seqid_2))

  # if there are no two sequences with sequence identity >= sequence_identity_cutoff, there are no shared
  # sequence groups
  if len(edges) < 1:
    return []

  # otherwise, use Python's networkx package to extract connected components from the graph defined by the
  # edges computed (unfortunately chained instances will be one group)
  shared_group_graph = Graph()  # create an empty graph
  shared_group_graph.add_edges_from(edges)  # add edges

  # list of lists with 2+ elements each corresponding to shared sequence groups
  final_groups = []
  for group in connected_components(shared_group_graph):
    final_groups.append(group)
  return final_groups


########################################################################################################

def create_shared_sequence_folds(seqid_to_sequence=None, num_folds=10, sequence_identity_cutoff=0.9):
  """
  :param seqid_to_sequence: dictionary of sequence ID -> sequence, where all sequences are the *same length*
  :param num_folds: maximum number of distinct, non-empty folds containing a roughly equal number of sequence IDs
  :param sequence_identity_cutoff: percent identity cutoff to assign sequences to the same "shared sequence" group
  :return: list of sets, where each set contains a subset of the sequence IDs from seqid_to_sequence
  """

  # group the domain instances into "shared sequence" groups with >90% identity
  shared_sequence_groups = create_shared_sequence_groups(seqid_to_sequence, sequence_identity_cutoff)
  all_sequence_groups = range(len(shared_sequence_groups))  # assign a "key" per group corresponding to an int

  # include the set of sequence IDs that did not fall into one of the shared sequence groups:
  grouped_structures = set()  # all the structure IDs that fell into a group
  for group in shared_sequence_groups:
    for sequence_id in group:
      grouped_structures.add(sequence_id)

  all_sequence_groups += sorted([sequence_id for sequence_id in seqid_to_sequence.keys()
                                 if sequence_id not in grouped_structures])
  total_shared_groups = len(all_sequence_groups)

  # split the groups into folds (no more than the total number of groups) as specified in the command line
  sequence_folds = np.array_split(np.random.permutation(all_sequence_groups),
                                  min(len(all_sequence_groups), num_folds))

  # convert each group ID back into the appropriate list of sequence IDs:
  final_folds = []
  for current_fold in sequence_folds:

    fold_members = set()
    for item in current_fold:
      if item not in seqid_to_sequence:  # i.e., this is a group number
        for sequence_id in shared_sequence_groups[int(item)]:  # recast from numpy string back to int
          fold_members.add(sequence_id)  # e.g., 1un6B_58_84
      else:
        fold_members.add(item)

    if len(fold_members) > 0:
      final_folds.append(fold_members)

  # return the final folds
  return final_folds, total_shared_groups, len(set(seqid_to_sequence.values()))


########################################################################################################

def create_folds_from_sequences(seqid_to_sequence=None, num_folds=10, sequence_identity_cutoff=0.9):
  """
  :param seqid_to_sequence: dictionary of sequence ID -> sequence, where all sequences are the *same length*
  :param num_folds: maximum number of distinct, non-empty folds containing a roughly equal number of sequence IDs
  :param sequence_identity_cutoff: percent identity cutoff to assign sequences to the same "shared sequence" group
  :return: list of sets, where each set contains a subset of the sequence IDs from seqid_to_sequence
  """

  # first, get the folds of SHARED SEQUENCE groups:
  (shared_sequence_folds,
   total_shared_groups,
   total_unique_sequences) = create_shared_sequence_folds(seqid_to_sequence, num_folds, sequence_identity_cutoff)

  # now, split all sequences, regardless of sequence identity, into folds, too:
  sequence_folds = np.array_split(np.random.permutation(seqid_to_sequence.keys()),
                                  min(len(seqid_to_sequence.keys()), num_folds))
  final_folds = []
  for current_fold in sequence_folds:
    fold_members = set()
    for item in current_fold:
      fold_members.add(item)
    if len(fold_members) > 0:
      final_folds.append(fold_members)

  # return the folds generated by grouping sequences by sequence identity AND the folds without grouping:
  return shared_sequence_folds, total_shared_groups, total_unique_sequences, final_folds


########################################################################################################
# FUNCTIONS TO AVERAGE LISTS
########################################################################################################

def arithmetic_mean(list_of_values):
  """
  :param list_of_values: list/set/array or collection of numeric values to average
  :return: the arithmetic mean of the input list of values
  """

  # confirm that all values are numeric:
  try:
    values = map(float, list_of_values)
    return sum(values)/max(float(len(list_of_values)), 1.)

  except ValueError:
    sys.stderr.write('Non-numeric value in: ['+','.join(map(str, list_of_values))+']\n')
    return 0.


########################################################################################################

def unravel_list(list_of_lists):
  """
  :param list_of_lists: an ordered list of ordered lists
  :return: a "flattened" (unraveled) singular list containing all elements in their original order
  """

  return [entry for sub_list in list_of_lists for entry in sub_list]


########################################################################################################

def macro_average_curve(xs, ys, curve_type='pr', print_errors=False):
  """
  :param xs: dictionary of sequence ID -> [x-values]
  :param ys: dictionary of sequence ID -> [y-values]
  :param curve_type: either 'pr' for precision-recall or 'roc' for receiver-operator-characteristic;
                     needed to pad beginning/end of (x,y) coordinates to properly interpolate y-values
  :param print_errors: boolean indicating whether to print error messages (useful when inputs are
                      in the wrong format)
  :return: two lists corresponding to union of all x-values and interpolated y-values
  """

  # edit pairwise points so that there is only one y-value per x-value and the start/end positions
  # are appropriate for interpolation
  for sequence_id in xs.keys():
    current_points = zip(xs[sequence_id], ys[sequence_id])

    xvals, yvals = [], []  # declare empty lists to be padded on either end as necessary afterward

    # for each unique x-value, choose the *maximum* corresponding y-value (for better interpolation)
    for xval in sorted(list(set(xs[sequence_id]))):
      xvals.append(xval)
      yvals.append(max([p[1] for p in current_points if p[0] == xval]))

    # check the start and end:
    if xvals[0] > 0:  # add the starting point 0 and the corresponding "minimum" y-value
      xvals = [0.] + xvals
      yvals = [yvals[0]] + yvals

    # NOTE: for precision-threshold curves; we may never have gotten a binding propensity > 0.9 etc,
    # so we don't want to claim that precision was actually low at these thresholds at all.
    if xvals[-1] < 1 and curve_type in ['pr', 'roc']:
      xvals += [1.]
      yvals += [1. if curve_type == 'roc' else yvals[-1]]

    # reset the stored xs and ys for this sequence ID:
    xs[sequence_id] = xvals
    ys[sequence_id] = yvals

  if print_errors:
    sys.stderr.write('X-values (should be a LIST of ARRAYS):\n')
    sys.stderr.write(str([xs[i] for i in xs.keys()]) + '\n')
    sys.stderr.write('Concatenate those distinct arrays (should be ONE array):\n')
    sys.stderr.write(str(unravel_list(xs.values())) + '\n')
    sys.stderr.write('Extract the unique elements (should still be ONE array):\n')
    sys.stderr.write(str(sorted(list(set(unravel_list(xs.values()))))) + '\n\n\n')

  # aggregate all unique x-values
  all_xs = np.unique(np.concatenate([xs[i] for i in xs.keys()]))

  # Then interpolate all corresponding y-values at these x coordinates
  all_ys = [0.]*len(all_xs)  # start with a y-value sum of 0 at each x-coordinate
  all_ys_cnt = [0.]*len(all_xs)  # ALSO keep track of how many non-negative values you even had..
  for sequence_id in xs.keys():
    for k, x in enumerate(all_xs):
      y = piecewise(x, xs[sequence_id], ys[sequence_id])
      if y < 0:
        continue
      all_ys[k] += y
      all_ys_cnt[k] += 1.

    # all_ys += np.piecewise(current_ys, [current_ys < 0, current_ys >= 0], [lambda xv: 0, lambda xv: xv])
    # all_ys_cnt += np.piecewise(current_ys, [current_ys < 0, current_ys >= 0], [lambda xv: 0, lambda xv: 1])
    # all_ys += np.interp(all_xs, np.array(xs[i], dtype="float"), np.array(ys[i], dtype="float"))

  # Finally average those y-values at each individual x-value
  all_ys = [y/all_ys_cnt[k] for k, y in enumerate(all_ys)]
  # all_ys /= all_ys_cnt

  # return lists of overall x-values and corresponding interpolated y-values
  return list(all_xs), list(all_ys)


########################################################################################################

def piecewise(x, xvals, yvals):
  """
  :param x: a particular x-value
  :param xvals: the x-values over which the piecewise function is defined
  :param yvals: the corresponding y-values at each x-value
  :return: evaluate the piecewise function
  """

  if x <= xvals[0]:
    return yvals[0]
  elif x > xvals[-1]:
    return -1 if yvals[-1] < 1 else yvals[-1]
  else:
    for i in xrange(len(xvals)-1):
      if xvals[i] < x <= xvals[i+1]:
        return yvals[i]


########################################################################################################
# ACCURACY FUNCTIONS
########################################################################################################

def auc_trapz(xvals, yvals, curve_type):
  """
  :param xvals: list of either false positive rates or recalls (floats), corresponding to x-values
  :param yvals: list of corresponding true positive rates or precisions (floats), corresponding to y-values
  :param curve_type: string, either 'pr' for a precision-recall curve, or 'roc' for a receiver-operator-
                     characteristic curve
  :return: the area under the specified curve, calculated using the trapezoidal function
  """

  assert len(xvals) == len(yvals), "X- and Y-values are not 1-to-1"

  # cannot calculate any area on empty lists of points:
  if len(xvals) == len(yvals) == 0:
    return 0.

  pts = sorted(zip(xvals, yvals))  # sort points in the X-direction.

  xs, ys = [], []
  for xval in sorted(list(set(xvals))):
    yval = max(map(float, [p[1] for p in pts if p[0] == xval]))
    xs.append(xval)
    ys.append(yval)

  if xs[0] > 0:
    xs = [0.] + xs
    ys = [ys[0]] + ys
  if xs[-1] < 1:
    xs += [1.]
    ys += [1. if curve_type == 'roc' else ys[-1]]  # all ROC curves must end at (1, 1)

  return np.trapz(ys, xs)  # calculate the return the area under the ROC curve


########################################################################################################

def receiver_operator_curve(binary_test_vector, continuous_training_vector):
  """
  :param binary_test_vector: list of booleans (True/False) or ints (1/0) corresponding to positive and 
                             negative examples
  :param continuous_training_vector: list of corresponding continuous values, where the higher values
                                     should indicate more positive examples
  :return: a sorted (low->high) list of pairs corresponding to false positive rate (x-value) and 
           true positive rate (y-value)
  """

  fprs, tprs, thresholds = [], [], []  # keep track of the paired false positive and true positive rates

  # total true positives and true negatives in this binary test vector
  total_positives = sum([1 if entry else 0 for entry in binary_test_vector])
  total_negatives = sum([1 if not entry else 0 for entry in binary_test_vector])

  if not (total_positives > 0 and total_negatives > 0):
    sys.stderr.write('Cannot compute TPR and FPR on all positive or all negative test vector: \n' +
                     ', '.join(map(str, binary_test_vector)) + '\n')
    sys.exit(1)

  # threshold specifically at the increments we are given in the continuous training vector
  sorted_thresholds = sorted(list(set(continuous_training_vector)), reverse=True) + [0.]  # high -> low

  # keep track of all binary positions at each *lower* threshold
  indices_by_threshold = {threshold: [index for index, entry in enumerate(continuous_training_vector)
                                      if entry >= threshold] for threshold in sorted_thresholds}

  for threshold in sorted_thresholds:
    current_matches = [binary_test_vector[index] for index in indices_by_threshold[threshold]]
    true_positives = sum([1 if entry else 0 for entry in current_matches])
    false_positives = sum([1 if not entry else 0 for entry in current_matches])

    fpr = float(false_positives) / total_negatives
    tpr = float(true_positives) / total_positives

    if fpr > 0 or tpr > 0:
      fprs.append(fpr)
      tprs.append(tpr)
      thresholds.append(threshold)

    if false_positives == total_negatives:
      break  # no need to explore further, we've already reached our maximum x-value

  # make sure our results make sense:
  if fprs[-1] < 1:
    sys.stderr.write('False positive rate does not finish at 1: ' + ', '.join(map(str, fprs)) + '\n')
    sys.exit(1)

  return fprs, tprs, thresholds


########################################################################################################

def calculate_receiver_operator(test_values, training_values):
  """
  :param test_values: list of binary values, where 1 indicates a positive example, and 0 indicates negative
  :param training_values: list of continuous values, where high values indicate possible positive examples,
                          and low values indicate possible negative examples
  :return: a sorted list of false positive rates (x), true positive rates (y), and the area under the receiver
           operator characteristic curve, calculated using the trapezoidal function
  """

  false_positive_rate, true_positive_rate, threshold = receiver_operator_curve(test_values, training_values)
  area_under_roc_curve = auc_trapz(false_positive_rate, true_positive_rate, 'roc')

  sorted_points = sorted(zip(false_positive_rate, true_positive_rate, threshold))  # sort the output from low->high
  return ([fpr for (fpr, _, _) in sorted_points],
          [tpr for (_, tpr, _) in sorted_points],
          area_under_roc_curve)


########################################################################################################

def precision_recall_curve(binary_test_vector, continuous_training_vector):
  """
  :param binary_test_vector: list of booleans (True/False) or ints (1/0) corresponding to positive and 
                             negative examples
  :param continuous_training_vector: list of corresponding continuous values, where the higher values
                                     should indicate more positive examples
  :return: a sorted (low->high) list of pairs corresponding to recall (x-value) and precision (y-value)
  """

  recalls, precisions, thresholds = [], [], []  # keep track of the paired recall and precision values

  # total true positives in this binary test vector
  total_positives = sum([1 if entry else 0 for entry in binary_test_vector])

  if not (total_positives > 0):
    sys.stderr.write('Cannot compute precision and recall on all positive or all negative test vector: \n' +
                     ', '.join(map(str, binary_test_vector)) + '\n')
    sys.exit(1)

  # threshold specifically at the increments we are given in the continuous training vector
  sorted_thresholds = sorted(list(set(continuous_training_vector)), reverse=True)  # high -> low
  if sorted_thresholds[-1] > 0:
    sorted_thresholds += [0.]

  # keep track of all binary positions at each *lower* threshold
  indices_by_threshold = {threshold: [index for index, entry in enumerate(continuous_training_vector)
                                      if entry >= threshold] for threshold in sorted_thresholds}

  for threshold in sorted_thresholds:
    current_matches = [binary_test_vector[index] for index in indices_by_threshold[threshold]]
    current_positives = sum([1 if entry else 0 for entry in current_matches])

    recall = float(current_positives) / total_positives
    precision = float(current_positives) / len(current_matches)

    if current_positives > 0:
      recalls.append(recall)
      precisions.append(precision)
      thresholds.append(threshold)

    if current_positives == total_positives:
      break  # no need to explore further, we've already reached our maximum x-value (recall)

  # make sure our results make sense:
  if not (recalls[0] > 0) or recalls[-1] < 1:
    sys.stderr.write('Recall starts at 0 and/or does not finish at 1: ' + ', '.join(map(str, recalls)) + '\n')
    sys.exit(1)

  return recalls, precisions, thresholds


########################################################################################################

def calculate_precision_recall(test_values, training_values):
  """
  :param test_values: list of binary values, where 1 indicates a positive example, and 0 indicates negative
  :param training_values: list of continuous values, where high values indicate possible positive examples,
                          and low values indicate possible negative examples
  :return: a sorted list of recall values (x), precision values (y), and the area under the precision
           recall curve, calculated using the trapezoidal function
  """

  # NOTE: we must compute the precision-recall ourselves, because we should never start with a recall of
  # zero, nor should we end with a precision of 0.

  recall, precision, _ = precision_recall_curve(test_values, training_values)
  area_under_pr_curve = auc_trapz(recall, precision, 'pr')

  sorted_points = sorted(zip(recall, precision))  # sort the output from low->high
  return ([rc for (rc, _) in sorted_points],
          [pr for (_, pr) in sorted_points],
          area_under_pr_curve)


########################################################################################################

def calculate_precision_threshold(test_values, training_values):
  """
  :param test_values: list of binary values, where 1 indicates a positive example, and 0 indicates negative
  :param training_values: list of continuous values, where high values indicate possible positive examples,
                          and low values indicate possible negative examples
  :return: a sorted list of recall values (x), precision values (y), and the area under the precision
           recall curve, calculated using the trapezoidal function
  """

  # NOTE: we must compute the precision-recall ourselves, because we should never start with a recall of
  # zero, nor should we end with a precision of 0.

  _, precision, threshold = precision_recall_curve(test_values, training_values)

  sorted_points = sorted(zip(threshold, precision))  # sort the output from high->low
  return ([th for (th, _) in sorted_points],
          [pr for (_, pr) in sorted_points],
          0.)


########################################################################################################

def consistency_by_splits(domain_name, num_splits=2, sequence_identity_cutoff=1.0,
                          distance='mindist', default_value=20.):
  """
  :param domain_name: full name of the domain being processed (e.g., PF00096_zf-C2H2)
  :param num_splits: how to split instances to calculate consistency between residue/ligand distances
  :param sequence_identity_cutoff: maximum allowed sequence identity for instances to be allowed in
                                   different splits
  :param distance: metric used to assess distance for binding potential scores
  :param default_value: "default" value stored in case of no information, depends on the distance
                        metric being used
  :return: list of output lines for consistency
  """

  average_consistency = []  # list of lines to write out to file

  score_file = SCORE_PATH+distance+'/'+domain_name+'_binding-scores_'+distance+'.txt.gz'

  # get the set of ligands that this domain binds to (e.g., DNABASE_, ATP, METABOLITE_):
  ligand_types = domain_ligand_types(score_file)
  total_processed_ligands = 0
  progress_bars = [(str(rank * 10) + '%', int(rank * (len(ligand_types) / 10.))) for rank in range(1, 10)][::-1]

  for ligand_type in sorted(list(ligand_types)):
    for progress_percent, progress_value in progress_bars:
      if total_processed_ligands > progress_value:
        progress_bars = progress_bars[:-1]  # remove the last value to not reprint
        sys.stderr.write('Processed ' + progress_percent + ' (' + "{:,}".format(total_processed_ligands) + '/' +
                         "{:,}".format(len(ligand_types)) + ') of ligands for '+domain_name+'.\n')
        break

    # full sequences for each domain instance (in contact with the current ligand type)
    alignment_file = ALN_PATH + distance + '/' + domain_name + '_' + ligand_type + '_' + distance + '.aln.fa'
    all_sequences = process_domain_alignment(alignment_file)  # sequence ID -> fully aligned sequence

    # create a per-sequence "distance-to-ligand" vector, i.e. sequence ID -> (values per position)
    distance_vectors = domain_distance_vectors(score_file, alignment_file, ligand_type, default_value)
    all_sequences = {seq_id: sequence for seq_id, sequence in all_sequences.items() if seq_id in distance_vectors}
    total_structures = len(set([seq_id[:4] for seq_id in all_sequences.keys()]))  # unique PDB IDs (without chains)

    # keep track of the fraction of positive (i.e. binding) positions to measure the "easiness" of the task
    average_fraction_positives = []
    for current_distance_vector in distance_vectors.values():
      binding_values = [True if (distance in ['mindist', 'meandist'] and current_distance <= PROXIMITY_CUTOFF) or
                                (distance not in ['mindist', 'meandist'] and current_distance > 0) else False
                        for current_distance in current_distance_vector]
      if True in binding_values and False in binding_values:  # some binding/non-binding positions
        average_fraction_positives.append(binding_values.count(True) / float(len(binding_values)))

    total_match_states = len(distance_vectors.values()[0])  # we don't need the actual match state names

    (shared_sequence_folds,
     total_shared_groups,
     total_unique_sequences,
     all_sequence_folds) = create_folds_from_sequences(all_sequences, num_splits, sequence_identity_cutoff)

    # update information for different types of curves that we are keeping track of:
    final_results_lines = [domain_name, ligand_type, str(total_match_states),
                           str(arithmetic_mean(average_fraction_positives)),
                           '|'.join([str(len(all_sequences.keys())),
                                     str(total_unique_sequences),
                                     str(total_shared_groups),
                                     str(total_structures)])]

    # group the domain instances into "shared sequence" groups with > specified sequence identity,
    # then return folds containing sequence IDs (ordered list of non-empty sets)
    for domain_instance_folds in [shared_sequence_folds, all_sequence_folds]:
      if len(domain_instance_folds) < 2:
        consistencies = '--'
        consistency = '--'
      else:
        consistency = []
        for _ in xrange(10):  # 10 random splits of the shared sequence groups
          folds = np.array_split(np.random.permutation(range(total_shared_groups)), num_splits)
          average_distances = []
          for half in folds:
            half_seq_ids = set()
            for grp_index in half:
              for seq_id in domain_instance_folds[grp_index]:
                half_seq_ids.add(seq_id)
            average_distances.append([arithmetic_mean([distance_vectors[seq_id][position]
                                                       for seq_id in half_seq_ids])
                                      for position in xrange(total_match_states)])

          for split_index in xrange(num_splits-1):
            for next_split_index in xrange(split_index+1, num_splits):
              pcc, _ = ss.pearsonr(average_distances[split_index], average_distances[next_split_index])
              consistency.append(pcc)
        consistencies = ','.join(map(str, consistency))
        consistency = str(arithmetic_mean(consistency))
      final_results_lines += [consistencies, consistency]

    # Write out results to file:
    average_consistency.append('\t'.join(final_results_lines) + '\n')
  total_processed_ligands += 1

  return average_consistency


########################################################################################################

def accuracy_by_cross_validation(domain_name, num_folds=10, sequence_identity_cutoff=0.9,
                                 distance='mindist', default_value=20.):
  """
  :param domain_name: full name of the domain being processed (e.g., PF00096_zf-C2H2)
  :param num_folds: number of folds to split the input data (domain instances) into for cross validation
  :param sequence_identity_cutoff: maximum allowed sequence identity between two domains for them to be considered
                                   in separate folds
  :param distance: metric used to assess distance for binding potential scores
  :param default_value: "default" value stored in case of no information, depends on the distance
                                 metric being used
  :return: list of output lines for precision-recall, list of output lines for receiver-operator, and 
           the final number of folds that was actually used
  """

  final_average_accuracy = {'pr': [], 'roc': [], 'pt': []}  # accuracy type -> list of lines to write out to file
  concatenate_vectors = True

  score_file = SCORE_PATH+distance+'/'+domain_name+'_binding-scores_'+distance+'.txt.gz'

  # get the set of ligands that this domain binds to (e.g., DNABASE_, ATP, METABOLITE_):
  ligand_types = domain_ligand_types(score_file)
  total_processed_ligands = 0
  progress_bars = [(str(rank * 10) + '%', int(rank * (len(ligand_types) / 10.))) for rank in range(1, 10)][::-1]

  for ligand_type in sorted(list(ligand_types)):

    for progress_percent, progress_value in progress_bars:
      if total_processed_ligands > progress_value:
        progress_bars = progress_bars[:-1]  # remove the last value to not reprint
        sys.stderr.write('Processed ' + progress_percent + ' (' + "{:,}".format(total_processed_ligands) + '/' +
                         "{:,}".format(len(ligand_types)) + ') of ligands for '+domain_name+'.\n')
        break

    # full sequences for each domain instance (in contact with the current ligand type)
    alignment_file = ALN_PATH + distance + '/' + domain_name + '_' + ligand_type + '_' + distance + '.aln.fa'
    all_sequences = process_domain_alignment(alignment_file)  # sequence ID -> fully aligned sequence
    total_structures = len(set([seq_id[:4] for seq_id in all_sequences.keys()]))  # unique PDB IDs (without chains)

    # create a per-sequence "distance-to-ligand" vector, i.e. sequence ID -> (values per position)
    distance_vectors = domain_distance_vectors(score_file, alignment_file, ligand_type, default_value)

    # keep track of the fraction of positive (i.e. binding) positions to measure the "easiness" of the task
    average_fraction_positives = []
    for current_distance_vector in distance_vectors.values():
      binding_values = [True if (distance in ['mindist', 'meandist'] and current_distance <= PROXIMITY_CUTOFF) or
                                (distance not in ['mindist', 'meandist'] and current_distance > 0) else False
                        for current_distance in current_distance_vector]
      if True in binding_values and False in binding_values:
        average_fraction_positives.append(binding_values.count(True) / float(len(binding_values)))

    total_match_states = len(distance_vectors.values()[0])  # we don't need the actual match state names

    # group the domain instances into "shared sequence" groups with > specified sequence identity,
    # then return folds containing sequence IDs (ordered list of non-empty sets)
    (shared_sequence_folds,
     total_shared_groups,
     total_unique_sequences,
     all_sequence_folds) = create_folds_from_sequences(all_sequences, num_folds, sequence_identity_cutoff)

    # update information for different types of curves that we are keeping track of:
    final_results_lines = {curve_type: [domain_name, ligand_type, str(total_match_states),
                                        str(arithmetic_mean(average_fraction_positives)),
                                        '|'.join([str(len(all_sequences.keys())),
                                                  str(total_unique_sequences),
                                                  str(total_shared_groups),
                                                  str(total_structures)])] for curve_type in ['pr', 'pt', 'roc']}

    # now, actually perform the cross validation !
    for domain_instance_folds in [shared_sequence_folds, all_sequence_folds]:

      # the following variables will be assigned values in the next steps:
      recall = {'Real': {}, 'Random': {}}  # xs
      precision = {'Real': {}, 'Random': {}}  # ys
      area_under_pr_curve = {'Real': {}, 'Random': {}}

      prop_thresholds = {'Real': {}, 'Random': {}}  # xs
      prop_precisions = {'Real': {}, 'Random': {}}  # ys
      area_under_prop_curve = {'Real': {}, 'Random': {}}

      false_positive_rate = {'Real': {}, 'Random': {}}  # xs
      true_positive_rate = {'Real': {}, 'Random': {}}  # ys
      area_under_roc_curve = {'Real': {}, 'Random': {}}

      # ---------------------------------------------------------------------------------------------------------
      # hold out each fold in turn, compute new scores on the training set, and test on the held-out set
      total_held_out_instances = 0
      for test_index, test_set in enumerate(domain_instance_folds):

        # set of domain instances to generate our positional weights from (training):
        training_set = set()
        for training_group in [other_group for train_index, other_group in enumerate(domain_instance_folds)
                               if train_index != test_index]:
          for training_item in training_group:
            training_set.add(training_item)

        # the sizes of the testing and training sets must be non-zero (here is where single folds normally fail)
        if len(test_set) < 1 or len(training_set) < 1:
          continue

        # get the "uniqueness" for each domain sequence in the training set
        sorted_seq_ids, relative_uniqueness = henikoff_alignment_score({sequence_id: all_sequences[sequence_id]
                                                                        for sequence_id in training_set})

        # method by which to "flatten" a distribution of values and relative weights (e.g., mindist):
        summary_value_function, _ = choose_summary_function(distance)

        # for each match state, get the "training" binding potential score
        training_values = []
        for match_state in xrange(total_match_states):
          current_distribution = [(distance_vectors[sequence_id][match_state],
                                   relative_uniqueness[seq_index]) for seq_index, sequence_id in
                                  enumerate(sorted_seq_ids) if sequence_id in distance_vectors]
          flattened_score = summary_value_function(current_distribution)
          training_values.append(flattened_score)

        if len(set(training_values)) < 2:  # Don't include training vector if it is all 0s
          continue

        # concatenate the test vectors and duplicate the corresponding training vectors
        concatenation_count = 0  # repetitions of our binding potential scores
        test_vector = []  # back-to-back concatenations of the binary binding test vectors
        for test_values in [values_tuple for sequence_id, values_tuple in distance_vectors.items()
                            if sequence_id in test_set]:
          current_test_vector = [True if (distance in ['mindist', 'meandist'] and
                                          current_distance <= PROXIMITY_CUTOFF) or
                                         (distance not in ['mindist', 'meandist'] and
                                          current_distance > 0) else False
                                 for current_distance in test_values]

          if len(set(current_test_vector)) < 2:  # Don't include test vector if it is all 0s or all 1s
            continue
          test_vector.append(current_test_vector)
          concatenation_count += 1

        # Don't bother testing on this "group" if it contained no individual vectors with some binding positions
        if len(test_vector) < 1:
          continue

        if concatenate_vectors:
          test_vector = [unravel_list(test_vector)]  # concatenate all the lists (should be just one for LOOCV)
          training_vector = [concatenation_count * training_values]  # and duplicate the training list accordingly
        else:
          training_vector = [training_values]*concatenation_count

        # randomize the training vectors 10x
        for hold_out_index in xrange(len(test_vector)):
          randomize_count = 10
          random_test_vector = randomize_count * test_vector[hold_out_index]  # repetitions of the true test labels
          random_training_vector = []  # back-to-back concatenations of randomized binding potential weights
          for _ in xrange(randomize_count):  # randomize 10 times!
            random_training_vector.append(list(np.random.permutation(training_vector[hold_out_index])))
          random_training_vector = unravel_list(random_training_vector)

          # calculate the "actual" & "random" precision-recall, precision-threshold, receiver-operator curves and area:
          for (xs, ys, auc, accuracy_func) in [(recall, precision, area_under_pr_curve, calculate_precision_recall),
                                               (prop_thresholds, prop_precisions, area_under_prop_curve,
                                                calculate_precision_threshold),
                                               (false_positive_rate, true_positive_rate,
                                                area_under_roc_curve, calculate_receiver_operator)]:

            for (value_type, test_values, training_values) in [('Real', test_vector[hold_out_index],
                                                                training_vector[hold_out_index]),
                                                               ('Random', random_test_vector, random_training_vector)]:
              (xs[value_type][total_held_out_instances],
               ys[value_type][total_held_out_instances],
               auc[value_type][total_held_out_instances]) = accuracy_func(test_values, training_values)
          total_held_out_instances += 1

      # average all values to print out
      for (xs, ys, auc, curve_type) in [(recall, precision, area_under_pr_curve, 'pr'),
                                        (prop_thresholds, prop_precisions, area_under_prop_curve, 'pt'),
                                        (false_positive_rate, true_positive_rate, area_under_roc_curve, 'roc')]:

        final_scores = []

        for value_type in ['Real', 'Random']:

          if len(xs[value_type]) > 1:  # multiple curves to be averaged; this is the usual case
            average_xs, average_ys = macro_average_curve(xs[value_type], ys[value_type], curve_type)
            average_auc = arithmetic_mean(auc[value_type].values())

          elif len(xs[value_type]) == 1:
            average_xs, average_ys = list(xs[value_type].values()[0]), list(ys[value_type].values()[0])
            average_auc = auc[value_type].values()[0]

          else:  # we were unable to cross validate this set (too few shared sequence groups, for instance)
            average_xs = average_ys = ['-']
            average_auc = 0

          final_scores.append(','.join(map(str, average_xs)))
          final_scores.append(','.join(map(str, average_ys)))
          final_scores.append(str(average_auc))

        # keep track of the final lines to print:
        real_xs, real_ys, real_auc, random_xs, random_ys, random_auc = final_scores[:6]
        final_results_lines[curve_type] += [real_xs, real_ys, random_xs, random_ys, real_auc, random_auc]

    # and now, store the final results line:
    for curve_type in ['pr', 'pt', 'roc']:
      final_average_accuracy[curve_type].append('\t'.join(final_results_lines[curve_type])+'\n')
    total_processed_ligands += 1

  # return all output lines and the average number of folds in the end
  return final_average_accuracy


########################################################################################################

def bootstrapped_stderr(domain_name, sequence_identity_cutoff, distance, default_value):
  """
  :param domain_name: full name of the domain being processed (e.g., PF00096_zf-C2H2)
  :param sequence_identity_cutoff: maximum allowed sequence identity between two domains for them to be considered
                                   in separate folds
  :param distance: metric used to assess distance for binding potential scores
  :param default_value: "default" value stored in case of no information, depends on the distance
                        metric being used
  :return: list of output lines for bootstrapped standard deviation
  """

  standard_errors = []  # list of lines to write out to file

  score_file = SCORE_PATH+distance+'/'+domain_name+'_binding-scores_'+distance+'.txt.gz'

  # method by which to "flatten" a distribution of values and relative weights:
  summary_value_function, _ = choose_summary_function(distance)

  # get the set of ligands that this domain binds to (e.g., DNABASE_, ATP, METABOLITE_):
  ligand_types = domain_ligand_types(score_file)
  total_processed_ligands = 0
  progress_bars = [(str(rank * 10) + '%', int(rank * (len(ligand_types) / 10.))) for rank in range(1, 10)][::-1]

  for ligand_type in sorted(list(ligand_types)):

    for progress_percent, progress_value in progress_bars:
      if total_processed_ligands > progress_value:
        progress_bars = progress_bars[:-1]  # remove the last value to not reprint
        sys.stderr.write('Processed ' + progress_percent + ' (' + "{:,}".format(total_processed_ligands) + '/' +
                         "{:,}".format(len(ligand_types)) + ') of ligands for '+domain_name+'.\n')
        break

    # full sequences for each domain instance (in contact with the current ligand type)
    alignment_file = ALN_PATH + distance + '/' + domain_name + '_' + ligand_type + '_' + distance + '.aln.fa'

    # NOTE: many of these sequences may be "cancelled out" if they do not contain a position within 3.6A of a ligand
    all_sequences = process_domain_alignment(alignment_file)  # sequence ID -> fully aligned sequence

    # create a per-sequence "distance-to-ligand" vector, i.e. sequence ID -> (values per position)
    distance_vectors = domain_distance_vectors(score_file, alignment_file, ligand_type, default_value)
    all_sequences = {seq_id: sequence for seq_id, sequence in all_sequences.items() if seq_id in distance_vectors}
    total_structures = len(set([seq_id[:4] for seq_id in all_sequences.keys()]))  # unique PDB IDs (without chains)
    distance_vectors = {seq_id: distance_vect for seq_id, distance_vect in distance_vectors.items()
                        if seq_id in all_sequences}

    # it's possible that we can't actually bootstrap these values
    if len(distance_vectors.keys()) < 1 or len(all_sequences.keys()) < 1:
      continue

    # keep track of the fraction of positive (i.e. binding) positions to measure the "easiness" of the task
    average_fraction_positives = []
    for current_distance_vector in distance_vectors.values():
      binding_values = [True if (distance in ['mindist', 'meandist'] and current_distance <= PROXIMITY_CUTOFF) or
                                (distance not in ['mindist', 'meandist'] and current_distance > 0) else False
                        for current_distance in current_distance_vector]
      if True in binding_values and False in binding_values:
        average_fraction_positives.append(binding_values.count(True) / float(len(binding_values)))

    total_match_states = len(distance_vectors.values()[0])  # we don't need the actual match state names

    # group the domain instances into "shared sequence" groups with > specified sequence identity,
    # then return folds containing sequence IDs (ordered list of non-empty sets)
    (shared_sequence_folds,
     total_shared_groups,
     total_unique_sequences,
     all_sequence_folds) = create_folds_from_sequences(all_sequences, len(all_sequences), sequence_identity_cutoff)

    # update information for different types of curves that we are keeping track of:
    final_results_lines = [domain_name, ligand_type, str(total_match_states),
                           str(arithmetic_mean(average_fraction_positives)),
                           '|'.join([str(len(all_sequences.keys())),
                                     str(total_unique_sequences),
                                     str(total_shared_groups),
                                     str(total_structures)])]

    all_seq_ids = sorted(distance_vectors.keys())
    bootstrapped_values = []
    for _ in xrange(1001):  # 1000 bootstraps
      if len(bootstrapped_values) < 1:  # store the actual binding propensities first
        current_seq_ids = all_seq_ids
      else:
        current_seq_ids = [all_seq_ids[i] for i in np.random.choice(len(all_seq_ids), len(all_seq_ids), replace=True)]

      # get the "uniqueness" for each domain sequence in the bootstrapped set
      sorted_seq_ids, relative_uniqueness = henikoff_alignment_score({sequence_id: all_sequences[sequence_id]
                                                                      for sequence_id in current_seq_ids})

      # for each match state, get the "training" binding propensity
      current_scores = []
      for match_state in xrange(total_match_states):
        current_distribution = [(distance_vectors[sequence_id][match_state],
                                 relative_uniqueness[seq_index]) for seq_index, sequence_id in
                                enumerate(sorted_seq_ids) if sequence_id in distance_vectors]
        flattened_score = summary_value_function(current_distribution)
        current_scores.append(flattened_score)
      bootstrapped_values.append(current_scores)

    # now, calculate the real binding propensities:
    real_binding_scores = ','.join(map(str, bootstrapped_values[0]))

    # the standard errors and 95% confidence interval sizes:
    std_errs, confidence_intervals = [], []
    for i in xrange(total_match_states):
      current_distribution = [bootstrapped_values[j][i] for j in range(1, 1001)]
      std_errs.append(np.std(current_distribution))
      current_distribution.sort()
      confidence_intervals.append(current_distribution[974]-current_distribution[25])
    std_errs = ','.join(map(str, std_errs))
    confidence_intervals = ','.join(map(str, confidence_intervals))

    # Write out results to file:
    final_results_lines.extend([real_binding_scores, std_errs, confidence_intervals])
    standard_errors.append('\t'.join(final_results_lines)+'\n')
  total_processed_ligands += 1

  return standard_errors


########################################################################################################
# OUTPUT
########################################################################################################

def format_outfile_header(num_folds, sequence_identity_cutoff, distance, curve_type):
  """
  :param num_folds: integer indicating the average number of folds across ligand types found in this file
  :param sequence_identity_cutoff: maximum allowed sequence identity between two domains for them to be considered
                                   in separate folds
  :param distance: distance metric used to calculate binding potential weights (e.g., mindist)
  :param curve_type: string indicating precision-recall (pr) or receiver-operator-characteristic (roc)
  :return: the complete header for the tab-delimited output file
  """

  # all output files have the same corresponding starting columns:
  starting_values = ['#domain_name',
                     'ligand_type',
                     'number_of_matchstates',
                     'fraction_binding_positions',
                     '|'.join(['instances',
                               'instances_unique',
                               'instances_<'+str(int(sequence_identity_cutoff*100))+'seq_identity',
                               'structures'])]

  if curve_type in ['stderr']:
    header = '\n'.join(['# Binding propensity standard errors and 95% confidence intervals computed from ' +
                        '1000 bootstrapped samples',
                        '# All domain sequences (aligned by ligand type) are in ' + ALN_PATH + distance + '/',
                        '# All binding propensities are found in ' + SCORE_PATH + distance + '/',
                        '\t'.join(starting_values +
                                  ['binding_propensities',
                                   'standard_errors',
                                   '95%_confidence_intervals'])]) + '\n'
    return header

  elif curve_type in ['cons']:
    header = '\n'.join(['# Consistency of residue-to-ligand distances across 10 randomized '+str(num_folds)+'-way ' +
                        '#   splits of domain instances computed using Pearson\'s Correlation Coefficient(PCC)',
                        '# All domain sequences (aligned by ligand type) are in ' + ALN_PATH + distance + '/',
                        '# All binding propensities are found in ' + SCORE_PATH + distance + '/',
                        '\t'.join(starting_values +
                                  ['individual_consistencies',
                                   'average_consistency'])]) + '\n'
    return header

  if curve_type == 'pr':
    x_name = 'recall'
    y_name = 'precision'
    au_name = 'auPRC'
  elif curve_type == 'roc':
    x_name = 'false_positive_rate'
    y_name = 'true_positive_rate'
    au_name = 'auROC'
  elif curve_type == 'pt':
    x_name = 'binding_propensity'
    y_name = 'precision'
    au_name = 'auPBC'
  else:
    x_name = 'x-value'
    y_name = 'y-value'
    au_name = 'AUC'

  header = '\n'.join(['# ~' + str(num_folds) + '-fold cross-validated '+x_name+'s and '+y_name+'s',
                      '# '+au_name+'s are calculated *before macro averaging* '+curve_type.upper()+' curves, ' +
                      'average of these values are given here',
                      '# All domain sequences (aligned by ligand type) are in ' + ALN_PATH + distance + '/',
                      '# All binding propensities are found in ' + SCORE_PATH + distance + '/',
                      '\t'.join(starting_values +
                                ['real_'+x_name,
                                 'real_'+y_name,
                                 'random_'+x_name,
                                 'random_'+y_name,
                                 'averaged_real_'+au_name,
                                 'averaged_random_'+au_name,
                                 'nogrouping_real_'+x_name,
                                 'nogrouping_real_' + y_name,
                                 'nogrouping_random_' + x_name,
                                 'nogrouping_random_' + y_name,
                                 'nogrouping_averaged_real_' + au_name,
                                 'nogrouping_averaged_random_' + au_name])]) + '\n'

  return header


########################################################################################################
# MAIN
########################################################################################################

if __name__ == "__main__":

  if not os.path.isdir(SCORE_PATH+'mindist/'):
    sys.stderr.write('No such directory: ' + SCORE_PATH + 'mindist/\n' +
                     'Please run: python generate_domain_scores.py\n')
    sys.exit(1)

  total_interaction_domains = len([domfile for domfile in
                                   os.listdir(SCORE_PATH+'mindist/')
                                   if domfile.endswith('_binding-scores_mindist.txt.gz')])

  if total_interaction_domains < 1:
    sys.stderr.write('No binding score files in ' + SCORE_PATH + 'mindist/\n' +
                     'Please run: python generate_domain_scores.py\n')
    sys.exit(1)

  # --------------------------------------------------------------------------------------------------
  # parse the command-line arguments
  parser = argparse.ArgumentParser(description='Assess the "accuracy" of our binding propensities by measuring ' +
                                               'how well scores generated from a subset of structures are ' +
                                               'indicative of binding in other structural instances.')

  parser.add_argument('--distance', type=str,
                      help='How to record the distance between receptor and ligand?',
                      default='mindist',
                      choices={'fracin4', 'mindist', 'meandist', 'maxstd', 'meanstd', 'sumstd',
                               'maxvdw', 'meanvdw', 'sumstd'})
  parser.add_argument('--num_folds', type=int,
                      help='Number of folds for cross-validation (to split training/test sets)',
                      default=10)  # for leave one out, specify a very high number (we take the min with # groups)
  parser.add_argument('--seq_identity', type=float,
                      help='Maximum (noninclusive) allowed sequence identity between domain instances to be ' +
                           'considered in the same cross-validated fold',
                      default=0.9)
  parser.add_argument('--consistency', dest='calculate_consistency', action='store_true', default=False,
                      help='Calculate distance-to-ligand consistencies across 50-50 splits of domain-ligand instances')
  parser.add_argument('--stderr', dest='calculate_stderr', action='store_true', default=False,
                      help='Calculate bootstrapped standard errors of all binding propensities')
  parser.add_argument('--start', type=int,
                      help='Starting index of sorted list of HMMs to cross validate scores for',
                      default=0,
                      choices=range(0, total_interaction_domains))
  parser.add_argument('--end', type=int,
                      help='Ending index of sorted list of HMMs to cross validate scores for',
                      default=total_interaction_domains,
                      choices=range(1, total_interaction_domains+1))

  args = parser.parse_args()
  sys.stderr.write('python cross_validate_scores.py --start '+str(args.start)+' --end '+str(args.end)+'\n')

  # --------------------------------------------------------------------------------------------------
  # get list of domains to process:
  for subdir in ['', args.distance]:
    if not os.path.isdir(SCORE_PATH+subdir):
      sys.stderr.write('No such directory: '+SCORE_PATH+subdir+'\n' +
                       'Please run: python generate_domain_scores.py --distance '+args.distance+'\n')
      sys.exit(1)

  domain_names = sorted([binding_score_file[:binding_score_file.find('_binding-scores_'+args.distance+'.txt.gz')]
                         for binding_score_file in os.listdir(SCORE_PATH+args.distance)])[args.start:args.end]

  # set the default "distance" value based on the distance metric indicated
  default_distance_value = 20. if args.distance in ['mindist', 'meandist'] else 0.

  # --------------------------------------------------------------------------------------------------
  # create directories to store output in:
  required_subdirs = ['cross_validation', 'cross_validation/'+args.distance]
  if args.calculate_consistency:
    required_subdirs.append('cross_validation/'+args.distance+'/consistency')
  elif args.calculate_stderr:
    required_subdirs.append('cross_validation/'+args.distance+'/standard_error')
  else:
    required_subdirs.extend(['cross_validation/'+args.distance+'/precision_recall',
                             'cross_validation/'+args.distance+'/receiver_operator',
                             'cross_validation/'+args.distance+'/precision_threshold'])
  for subdir in required_subdirs:
    if not os.path.isdir(DATAPATH+'processed_data/domains/'+subdir):
      call(['mkdir', DATAPATH+'processed_data/domains/'+subdir])

  # --------------------------------------------------------------------------------------------------
  # process each interaction domain, one at a time
  for current_domain_name in domain_names:
    sys.stderr.write('Processing '+current_domain_name+'...\n')

    accuracy_by_ligand = {}  # keep track of the accuracy for each domain-ligand pair
    outfiles = {}  # full paths to corresponding output files

    # distance-to-ligand consistencies across 50-50 splits of domain-ligand instances
    if args.calculate_consistency:
      accuracy_by_ligand['con'] = consistency_by_splits(current_domain_name,
                                                        2,  # number of folds = two
                                                        0.9,  # sequence identity = 90%
                                                        args.distance,
                                                        default_distance_value)
      outfiles['con'] = CV_PATH + args.distance + '/consistency/' + current_domain_name + '-con-' + \
                        args.distance+'.txt.gz'

    # bootstrapped standard errors
    elif args.calculate_stderr:
      accuracy_by_ligand['stderr'] = bootstrapped_stderr(current_domain_name,
                                                         1.0,  # sequence identity = 100%
                                                         args.distance,
                                                         default_distance_value)
      outfiles['stderr'] = CV_PATH + args.distance + '/standard_error/' + current_domain_name + '-stderr-' + \
                           args.distance + '.txt.gz'

    # precision-recall, receiver-operator characteristic, precision-threshold (default)
    else:
      accuracy_by_ligand = accuracy_by_cross_validation(current_domain_name,
                                                        args.num_folds,
                                                        args.seq_identity,
                                                        args.distance,
                                                        default_distance_value)
      outfiles = {score_type: CV_PATH + args.distance + '/' +
                              ('precision_recall' if score_type == 'pr' else
                               ('receiver_operator' if score_type == 'roc' else
                                'precision_threshold')) + '/' +
                              current_domain_name + '-' + score_type + '-' + args.distance + '.txt.gz'
                  for score_type in ['pr', 'roc', 'pt']}

    # --------------------------------------------------------------------------------------------------
    # write out results
    for score_type in accuracy_by_ligand.keys():
      accuracy_handle = gzip.open(outfiles[score_type], 'w') if outfiles[score_type].endswith('gz') else \
                        open(outfiles[score_type], 'w')
      accuracy_handle.write(format_outfile_header(str(args.num_folds), args.seq_identity, args.distance, score_type))
      for result_line in sorted(accuracy_by_ligand[score_type]):
        accuracy_handle.write(result_line)
      accuracy_handle.close()
      sys.stderr.write('See ' + outfiles[score_type] + '\n')
