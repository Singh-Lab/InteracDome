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

# path to where this script is currently located (and to where all data should be stored) -- this can
# be updated
DATAPATH = os.path.dirname(os.path.abspath(__file__))+'/'
DISTANCE_CUTOFF = 20.  # maximum distance (in Angstroms) between a ligand and receptor atom to record


########################################################################################################

def mean(numbers):
  """
  :param numbers: list, set, or collection of numbers (positive, negative floats/ints)
  :return: the average of the values in the input number list
  """
  return float(sum(numbers)) / max(len(numbers), 1)


########################################################################################################

def receptor_ligand_distance(euc_dist, vdw_dist, vdw_error, std_dist, std_error, distance):
  """
  :param euc_dist: Euclidean distance between two atoms specified by their (x, y, z) coordinates
  :param vdw_dist: overlap area between two Gaussians centered at teh two (x, y, z) atom locations, with 
                   standard deviations set the respective van der Waals' interaction radii of the two atoms
  :param vdw_error: computation error from calculating the overlap area (integral from -Inf to Inf)
  :param std_dist: as in vdw_dist, but where the standard deviations of the Gaussians are set to 1.5
  :param std_error: as in vdw_error, except between two Gaussians with standard deviations set to 1.5
  :param distance: metric used to determine the distance between the two atoms
  :return: the "distance" between the two atoms according to the specified distance metric
  """

  # Euclidean distance to a particular ligand atom
  if distance in ['mindist', 'meandist', 'fracin4']:
    current_distance = float(euc_dist)

  # Overlap area (integrated) between two normal distributions with 1.5 standard deviation across the board
  elif distance in ['maxstd', 'meanstd', 'sumstd']:
    if std_dist and std_dist != 'N/A' and std_error and float(std_error) < float(std_dist):
      current_distance = float(std_dist)
    else:
      current_distance = 'N/A'

  # Overlap area (integrated) between two normal distributions with vdw interaction radii as standard deviation
  elif distance in ['maxvdw', 'maxstd', 'sumvdw']:
    if vdw_dist and vdw_dist != 'N/A' and vdw_error and float(vdw_error) < float(vdw_dist):
      current_distance = float(vdw_dist)
    else:
      current_distance = 'N/A'

  else:
    current_distance = 'N/A'

  return current_distance


########################################################################################################

def residue_binding_score(residue, ligand_type, distance):
  """
  :param residue: list of tuples (atom ID, ligand atom type, distance)
  :param ligand_type: type of ligand that we are currently considering
  :param distance: scoring metric for condensing information about all heavy side-chain residue atoms 
                   into a single score
  :return: score for the residue according to the specified distance metric
  """

  # final score for this particular residue_aa_index to the current ligand type
  residue_score = 0. if distance != 'mindist' else DISTANCE_CUTOFF

  if distance == 'fracin4':  # FRACTION of atoms that are within 4AA of THIS ligand.
    total_atoms = [float(res_atom_id.split('/')[1])
                   for (res_atom_id, lig_id, _) in residue if lig_id == ligand_type][0]
    residue_score = len(set([res_atom_id for (res_atom_id, lig_id, curr_dist) in residue
                             if lig_id == ligand_type and float(curr_dist) <= 4])) / total_atoms

  elif distance == 'mindist':  # MINIMUM distance of ANY atom to this ligand...
    residue_score = min([curr_dist for (_, lig_id, curr_dist) in residue if lig_id == ligand_type])

  elif distance == 'meandist':  # AVERAGE distance of ANY atom to a ligand.
    scores = []
    for residue_atom_id in set([res_atom_id for (res_atom_id, _, _) in residue]):
      current_score = [curr_dist for (res_atom_id, lig_id, curr_dist) in residue
                       if lig_id == ligand_type and res_atom_id == residue_atom_id]
      scores.append(min(current_score) if len(current_score) > 0 else DISTANCE_CUTOFF)
    residue_score = mean(scores)

  elif distance in ['maxstd', 'maxvdw']:  # MAXIMUM "closeness" score of ANY atom
    residue_score = 0.
    for residue_atom_id in set([res_atom_id for (res_atom_id, _, _) in residue]):
      current_score = sum([curr_dist for (res_atom_id, lig_id, curr_dist) in residue
                           if lig_id == ligand_type and res_atom_id == residue_atom_id])
      if current_score > residue_score:
        residue_score = current_score

  elif distance in ['meanstd', 'meanvdw']:  # AVERAGE "closeness" score across all atoms
    scores = []
    for residue_atom_id in set([res_atom_id for (res_atom_id, _, _) in residue]):
      current_score = sum([curr_dist for (res_atom_id, lig_id, curr_dist) in residue
                           if lig_id == ligand_type and res_atom_id == residue_atom_id])
      scores.append(current_score)
    residue_score = mean(scores)

  elif distance in ['sumstd', 'sumvdw']:  # SUM of "closeness" scores of all atoms!
    # Note: sum of an empty list returns 0
    residue_score = sum([curr_dist for (_, lig_id, curr_dist) in residue if lig_id == ligand_type])

  return residue_score


########################################################################################################

def create_biolip_fasta_files(distance_file, fasta_file, current_pdb_id, distance):
  """
  :param distance_file: full path to a "_distances.txt.gz" file as created by calculate_distances.py
  :param fasta_file: full path to an output "_distancetype.fa" where output should be written to
  :param current_pdb_id: PDB ID that the distance file should correspond to
  :param distance: type of per-residue score measuring the binding status of an amino acid residue.
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
  receptor_sequence = ''
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
    current_distance = receptor_ligand_distance(euc_dist, vdw_dist, vdw_error, std_dist, std_error, distance)

    try:
      new_atom_info = (residue_atom_id,
                       rename_ligand_id.get(ligand_id, 'SM' if len(rename_ligand_id.keys()) > 0 else ligand_id),
                       float(current_distance))
      binding_residues[pdb_chain_id][int(residue_aa_index)].append(new_atom_info)
    except ValueError:  # if distance cannot be converted to a float, then we do not keep track of it
      continue

  distance_handle.close()

  # create the fasta file as specified:
  fasta_outhandle = gzip.open(fasta_file, 'w') if fasta_file.endswith('gz') else open(fasta_file, 'w')

  for pdb_chain_id in sorted(receptor_chain_sequences.keys()):
    if pdb_chain_id not in binding_residues:
      continue

    fasta_outhandle.write('>' + pdb_chain_id + ' bindingSiteRes=')

    final_binding_residues = []
    for residue_aa_index, current_binding_residue in binding_residues[pdb_chain_id].items():
      for ligand_id in set([lig_id for (_, lig_id, _) in current_binding_residue]):

        current_residue_score = residue_binding_score(current_binding_residue, ligand_id, distance)

        final_binding_residues.append((residue_aa_index, ligand_id, current_residue_score))
    fasta_outhandle.write(','.join([str(res_aa_index) + '-' + lig_id + '-' + str(final_score)
                                    for (res_aa_index, lig_id, final_score) in sorted(final_binding_residues)]) +
                          ';\n' + receptor_chain_sequences[pdb_chain_id] + '\n\n')
  fasta_outhandle.close()
  sys.stderr.write('Wrote to ' + fasta_file + '\n')


########################################################################################################

def create_hmmer_fasta_file(annotation_dir=DATAPATH+'downloaded_data/annotations/',
                            annotation_file=DATAPATH+'processed_data/annotations/current_annotations.txt'):
  """
  :param annotation_dir: full path to a directory containing annotation files downloaded from BioLiP
  :param annotation_file: full path to a file containing tab-delimited information about PDB identifiers,
                          chains, and sequences
  :return: None, but print success message upon successful write of FASTA-formatted file (to find domains)
  """

  # first, get the date of the last downloaded annotation file:
  year, month, day = sorted([map(int, fname[fname.find('BioLiP_')+7:fname.rfind('.txt')].split('-'))
                             for fname in os.listdir(annotation_dir)
                             if fname.startswith('BioLiP_') and fname.endswith('.txt')])[-1]
  latest_date = str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)

  out_fasta_file = DATAPATH+'processed_data/annotations/BioLiP_'+latest_date+'_nonredundant.fa'
  recorded_chains = {}  # keep track of which chains have already been written to the file

  annot_handle = open(annotation_file)
  for annot_line in annot_handle:
    seq_id = ''.join(annot_line.strip().split()[:2])
    seq = annot_line.strip().split()[-1]

    if seq not in recorded_chains:
      recorded_chains[seq] = set()
    recorded_chains[seq].add(seq_id)
  annot_handle.close()

  # write nonredundant results to outfile
  fasta_handle = open(out_fasta_file, 'w')
  for seq_index, (seq, seq_ids) in sorted(recorded_chains):
    fasta_handle.write('>seq'+str(seq_index)+' '+','.join(sorted(list(seq_ids)))+'\n'+seq+'\n')
  fasta_handle.close()

  sys.stderr.write('All protein receptor chains (as of '+latest_date+') are found in '+out_fasta_file+'\n')


########################################################################################################

if __name__ == "__main__":

  # Parse the command-line arguments
  parser = argparse.ArgumentParser(description='Process the receptor-ligand distance files to create FASTA ' +
                                               'files where each AA position is associated with a binding score ' +
                                               'to 1+ ligand types.')

  parser.add_argument('--prefix', type=str,
                      help='Two letter prefix to subset of PDB IDs to process',
                      default='10')
  parser.add_argument('--force', dest='force', action='store_true',
                      help='Forcibly overwrite fasta files that have already been written; otherwise skip')
  parser.add_argument('--hmmer_input', dest='hmmer_input', action='store_true',
                      help='Create single FASTA file to run domain-finding algorithm on')
  parser.add_argument('--distance', type=str,
                      help='How to record the distance between receptor and ligand?',
                      default='mindist',
                      choices={'mindist', 'fracin4', 'meandist', 'maxstd', 'meanstd', 'sumstd',
                               'maxvdw', 'meanvdw', 'sumvdw'})
  parser.set_defaults(force=False)
  parser.set_defaults(hmmer_input=False)

  args = parser.parse_args()

  # ----------------------------------------------------------------------------------------------------
  # confirm that a concatenated list of all current BioLiP annotations exists:
  current_annotation_file = DATAPATH+'processed_data/annotations/current_annotations.txt'

  if not os.path.isfile(current_annotation_file):
    sys.stderr.write('Could not open '+current_annotation_file+'\n')
    sys.stderr.write('Please run python '+DATAPATH+'/download_biolip.py\n')
    sys.exit(1)

  if args.hmmer_input:
    create_hmmer_fasta_file(DATAPATH+'downloaded_data/annotations/', current_annotation_file)
    sys.exit(0)

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
    fasta_outfile = fasta_out_directory+infile.replace('_distances.txt.gz', '_'+args.distance+'.fa.gz')

    # if this file already exists, and we're not forcibly overwriting, then skip!
    if os.path.isfile(fasta_outfile) and not args.force:
      continue
      
    create_biolip_fasta_files(distance_infile, fasta_outfile, infile[27:31], args.distance)

  sys.stderr.write('Successfully wrote '+"{:,}".format(len(distlist_files))+' '+args.distance+' fasta files to ' +
                   fasta_out_directory+'\n')
