#!/usr/bin/python

"""
Use the receptor atom -> ligand atom distances computed using calculate_distances.py to create per-PDB fasta files
where the sequence header line contains the following comma-delineated binding site information: 1-index residue 
position:type of ligand:score.

Contact Shilpa Nadimpalli Kobren (snadimpa@princeton.edu) with questions
"""

import os
import sys
import gzip
import argparse
from subprocess import call


########################################################################################################
# CONSTANTS
########################################################################################################

DATAPATH = os.getcwd()  # path to where all data should be stored; this can be updated
DISTANCE_CUTOFF = 20.  # maximum distance (in Angstroms) between a ligand and receptor atom to record


########################################################################################################

def mean(numbers):
  """
  :param numbers: list, set, or collection of numbers (positive, negative floats/ints)
  :return: the average of the values in the input number list
  """
  return float(sum(numbers)) / max(len(numbers), 1)


########################################################################################################

def create_biolip_fasta_files(distance_file, fasta_file, current_pdb_id, distance_metric):
  """
  :param distance_file: full path to a "_distances.txt.gz" file as created by calculate_distances.py
  :param fasta_file: full path to an output "_distancetype.fa" where output should be written to
  :param current_pdb_id: PDB ID that the distance file should correspond to
  :param distance_metric: type of per-residue score measuring the binding status of an amino acid residue.
                          (e.g., mindist, meandist, fracin4, maxstd, meanstd, sumstd, maxvdw, meanvdw, sumvdw)
  :return: full path to a recently created fasta file where the sequence ID header line contains binding info
  """

  if os.path.getsize(distance_file) <= 180:
    sys.stderr.write('Incomplete file: ' + distance_file + '\n')
    return

  receptor_chain_sequences = {}  # pdbID-pdbChain -> complete corresponding AA sequence
  binding_residues = {}  # pdbID-pdbChain -> AA residue index -> [(atom ID, ligand_type, "distance")]
  rename_ligand_id = {}  # original ligand ID -> newly named ligand ID (in this case, unnecessary)

  distance_handle = gzip.open(distance_file) if distance_file.endswith('gz') else open(distance_file)
  for distline in distance_handle:
    if distline.startswith('#'):
      continue

    (pdb_chain_id, residue_aa_index, residue_aa_value, residue_atom_id, residue_atom_value,
     ligand_id, ligand_atom_value, euc_dist, vdw_dist, vdw_error, std_dist, std_error) = distline[:-1].split('\t')[:12]

    if current_pdb_id != '' and pdb_chain_id[:-1] != current_pdb_id:
      continue

    # get the receptor sequence (if this line contains that information)
    if len(distline[:-1].split('\t')) > 12 and distline[:-1].split('\t')[12] != '':
      receptor_sequence = distline[:-1].split('\t')[12]

    # make sure that the residue at the position is what we expect, and that we haven't exceeded the size of
    #  the receptor change. This should never happen (unless there is a bug to investigate)
    if int(residue_aa_index) - 1 >= len(receptor_sequence) or \
            receptor_sequence[int(residue_aa_index) - 1] != residue_aa_value:
      sys.stderr.write('ERROR: [' + residue_aa_index + '] != ' + residue_aa_value + ' in ' + receptor_sequence + '\n')
      continue

    # store the sequence (the sequence *will* appear on the first line of a new pdbID-chain
    if pdb_chain_id not in receptor_chain_sequences:
      receptor_chain_sequences[pdb_chain_id] = receptor_sequence

    # start storing atom distances for this particular receptor chain residue:
    if pdb_chain_id not in binding_residues:
      binding_residues[pdb_chain_id] = {}
    if int(residue_aa_index) not in binding_residues[pdb_chain_id]:
      binding_residues[pdb_chain_id][int(residue_aa_index)] = []

    # set the current_distance accordingly:
    current_distance = None

    # Euclidean distance to a particular ligand atom
    if distance_metric in ['mindist', 'meandist', 'fracin4']:
      current_distance = float(euc_dist)

    # Overlap area (integrated) between two normal distributions with 1.5 standard deviation across the board
    elif distance_metric in ['maxstd', 'meanstd', 'sumstd']:
      if std_dist and std_dist != 'N/A' and std_error and float(std_error) < float(std_dist):
        current_distance = float(std_dist)
      else:
        current_distance = 'N/A'

    # Overlap area (integrated) between two normal distributions with vdw interaction radii as standard deviation
    elif distance_metric in ['maxvdw', 'maxstd', 'sumvdw']:
      if vdw_dist and vdw_dist != 'N/A' and vdw_error and float(vdw_error) < float(vdw_dist):
        current_distance = float(vdw_dist)
      else:
        current_distance = 'N/A'

    try:
      new_atom_info = (residue_atom_id,
                       rename_ligand_id.get(ligand_id, 'SM' if len(rename_ligand_id.keys()) > 0 else ligand_id),
                       float(current_distance))
      binding_residues[pdb_chain_id][int(residue_aa_index)].append(new_atom_info)
    except ValueError:  # if distance cannot be converted to a float, then we do not keep track of it
      continue
  distance_handle.close()

  # create the fasta file as specified:
  fasta_outhandle = open(fasta_file, 'w')

  for pdb_chain_id in sorted(receptor_chain_sequences.keys()):
    if pdb_chain_id not in binding_residues:
      continue

    fasta_outhandle.write('>' + pdb_chain_id + ' bindingSiteRes=')

    final_binding_residues = []
    for residue_aa_index, current_binding_residues in binding_residues[pdb_chain_id].items():
      for ligand_id in set([lig_id for (_, lig_id, _) in current_binding_residues]):

        v = None  # final score for this particular residue_aa_index to the current ligand type

        if distance_metric == 'fracin4':  # FRACTION of atoms that are within 4AA of THIS ligand.
          total_atoms = [float(res_atom_id.split('/')[1])
                         for (res_atom_id, lig_id, _) in current_binding_residues if lig_id == ligand_id][0]
          v = len(set([res_atom_id for (res_atom_id, lig_id, curr_dist) in current_binding_residues
                       if lig_id == ligand_id and float(curr_dist) <= 4])) / total_atoms

        elif distance_metric == 'mindist':  # MINIMUM distance of ANY atom to this ligand...
          v = min([curr_dist for (_, lig_id, curr_dist) in current_binding_residues if lig_id == ligand_id])

        elif distance_metric == 'meandist':  # AVERAGE distance of ANY atom to a ligand.
          scores = []
          for residue_atom_id in set([res_atom_id for (res_atom_id, _, _) in current_binding_residues]):
            current_score = [curr_dist for (res_atom_id, lig_id, curr_dist) in current_binding_residues
                             if lig_id == ligand_id and res_atom_id == residue_atom_id]
            scores.append(min(current_score) if len(current_score) > 0 else DISTANCE_CUTOFF)
          v = mean(scores)

        elif distance_metric in ['maxstd', 'maxvdw']:  # MAXIMUM "closeness" score of ANY atom
          v = 0
          for residue_atom_id in set([res_atom_id for (res_atom_id, _, _) in current_binding_residues]):
            current_score = sum([curr_dist for (res_atom_id, lig_id, curr_dist) in current_binding_residues
                                 if lig_id == ligand_id and res_atom_id == residue_atom_id])
            if current_score > v:
              v = current_score

        elif distance_metric in ['sumstd', 'sumvdw']:  # SUM of "closeness" scores of all atoms!
          # Note: sum of an empty list returns 0
          v = sum([curr_dist for (_, lig_id, curr_dist) in current_binding_residues if lig_id == ligand_id])

        final_binding_residues.append((residue_aa_index, ligand_id, v))
    fasta_outhandle.write(','.join([str(res_aa_index) + '-' + lig_id + '-' + str(final_score)
                                    for (res_aa_index, lig_id, final_score) in sorted(final_binding_residues)]) +
                          ';\n' + receptor_chain_sequences[pdb_chain_id] + '\n\n')
  fasta_outhandle.close()
  sys.stderr.write('Wrote to ' + fasta_file + '\n')


########################################################################################################

if __name__ == "__main__":
  # Parse the command-line arguments
  parser = argparse.ArgumentParser(description='Process the receptor-ligand distance files to create FASTA ' +
                                               'files where each AA position is associated with a binding score ' +
                                               'to 1+ ligand types.')

  parser.add_argument('--prefix', type=str,
                      help='Two letter prefix to subset of PDB IDs to process',
                      default='10')
  parser.add_argument('--distance_metric', type=str,
                      help='How to record the distance between receptor and ligand?',
                      default='mindist',
                      choices={'mindist', 'fracin4', 'meandist', 'maxstd', 'meanstd', 'sumstd',
                               'maxvdw', 'meanvdw', 'sumstd'})

  args = parser.parse_args()

  # ----------------------------------------------------------------------------------------------------
  # confirm that a concatenated list of all current BioLiP annotations exists:
  current_annotation_file = DATAPATH+'processed_data/annotations/current_annotations.txt'

  if not os.path.isfile(current_annotation_file):
    sys.stderr.write('Could not open '+current_annotation_file+'\n')
    sys.stderr.write('Please run python '+os.getcwd()+'/download_biolip.py\n')
    sys.exit(1)

  # find the subset of PDB IDs to run on:
  subset_pdb_ids = set()
  pdb_id_handle = open(current_annotation_file)
  for annot_line in pdb_id_handle:
    this_pdb_id = annot_line.strip().split('\t')[0]
    if this_pdb_id.startswith(args.prefix):
      subset_pdb_ids.add(this_pdb_id)
  pdb_id_handle.close()

  if len(subset_pdb_ids) < 1:
    sys.stderr.write('No PDB IDs to process with prefix ' + args.prefix + '!\n')
    sys.exit(1)

  # ----------------------------------------------------------------------------------------------------
  # create a FASTA file for each receptor protein sequence that contains the distance to a particular
  #  ligand

  # check distance output directory (make sure it exists)
  distance_out_directory = DATAPATH + 'processed_data/distances/' + args.prefix[0] + '/' + args.prefix + '/'
  if not os.path.isdir(distance_out_directory):
    sys.stderr.write('No such directory: '+distance_out_directory+'\n')
    sys.exit(1)

  # create fasta output directory (if need be)
  for subdir in ['', args.prefix[0], args.prefix[0] + '/' + args.prefix]:
    if not os.path.isdir(DATAPATH + 'processed_data/fasta/' + subdir):
      call(['mkdir', DATAPATH + 'processed_data/fasta/' + subdir])
  fasta_out_directory = DATAPATH + 'processed_data/fasta/' + args.prefix[0] + '/' + args.prefix + '/'

  distlist_files = sorted([biolip_file for biolip_file in os.listdir(distance_out_directory)
                           if biolip_file.endswith('_distances.txt.gz') and
                           biolip_file.split('_')[0] in subset_pdb_ids])

  for infile in distlist_files:
    distance_infile = distance_out_directory + infile
    fasta_outfile = fasta_out_directory+infile.replace('_distances.txt', '_'+args.distance_metric+'.fa')
    create_biolip_fasta_files(distance_infile, fasta_outfile, infile[27:31], args.distance_metric)
