#!/usr/bin/python

"""
Create "distance to ligand" files for each BioLiP PDB that we have -- functions have been
updated to be run on the cluster

Contact Shilpa Nadimpalli Kobren (snadimpa@princeton.edu) with questions
"""

import os
import sys
import math
import gzip
import argparse
from subprocess import call
from string import ascii_letters, digits
from random import choice


########################################################################################################
# CONSTANTS
########################################################################################################

DATAPATH = os.getcwd()  # path to where all data should be stored; this can be updated
DISTANCE_CUTOFF = 20.  # maximum distance (in Angstroms) between a ligand and receptor atom to record


########################################################################################################
# STEP 1: DISTANCES BETWEEN LIGANDS
########################################################################################################

def random_filename(size=10, chars=ascii_letters+digits):
  """
  :param size: length of random string to return
  :param chars: set of characters from which to randomly select (with replacement)
  :return: a random string of length size comprised of chars (lowercase ASCII letters and digits 0-9)
  """

  return ''.join(choice(chars) for _ in range(size))


########################################################################################################

def get_vdw_interaction_radii(vdw_filename=DATAPATH+'downloaded_data/vdw.txt'):
  """
  :param vdw_filename: full path to a tab-delineated file with the columns:
  :return: dictionary of element -> van der Waals' interaction radii (or ionic radii, if the interaction
           radius is unavailable)
  """

  vdw = {}  # element -> van der Waals interaction radii

  vdw_inhandle = gzip.open(vdw_filename) if vdw_filename.endswith('gz') else open(vdw_filename)
  for vdw_line in vdw_inhandle:
    if vdw_line.startswith('#'):
      continue
    element = vdw_line.strip().split('\t')[1].upper()
    ionic = vdw_line.strip().split('\t')[3] if len(vdw_line.strip().split('\t')) > 3 else ''
    vandw = vdw_line.strip().split('\t')[5] if len(vdw_line.strip().split('\t')) > 5 else ''

    vdw[element] = vandw if vandw != '' else (ionic if ionic != '' else '1.5')
  vdw_inhandle.close()

  return vdw


########################################################################################################

def get_euclidean_distance(pt_1, pt_2):
  """
  :param pt_1: n-length tuple corresponding to the first point
  :param pt_2: n-length tuple corresponding to the second point
  :return: Euclidean distance between the two points pt_1 and pt_2
  """

  assert len(pt_1) == len(pt_2), "Points (" + str(pt_1) + ") and (" + str(pt_2) + ") have different dimensions."

  return math.sqrt(sum([(pt_1[i] - pt_2[i]) * (pt_1[i] - pt_2[i]) for i in xrange(len(pt_1))]))


########################################################################################################

def subtype_nuc_ligand(ligand_structure_file):
  """
  :param ligand_structure_file: full path to a ligand structural PDB file (downloaded from BioLiP)
  :return: (string) full path to a new, temporary copy of the original ligand structure file with an 
           extra column specifying whether the atom line corresponds to RNA, DNA, backbone, or base
  """

  # check that the file exists and that it is of type "nucleic acid"
  if not os.path.isfile(ligand_structure_file) or '_NUC_' not in ligand_structure_file:
    return ligand_structure_file

  # name of the temporary, updated file that we will be writing to
  tmp_ligand_file = DATAPATH + ligand_structure_file.split('/')[-1] + '-tmpnuc'
  tmp_nuchandle = open(tmp_ligand_file, 'w')

  curr_base_id = ''  # keep track of the current base ID
  curr_base_lines = []  # and all lines from the original ligand file corresponding to that base

  orig_ligand_handle = open(ligand_structure_file)
  for lig_line in orig_ligand_handle:

    # inappropriately short / non-ATOM entry line:
    if len(lig_line) < 81:
      tmp_nuchandle.write(lig_line)
      continue

    # reached a new base, so process and reset:
    if curr_base_id != lig_line[22:26].strip():

      if len(curr_base_lines) > 0:

        # "deoxy-ribonucleic acid" (DNA) will not have any O2 backbone atoms
        atoms = [base_line[12:16].strip() for base_line in curr_base_lines]
        nucleic_acid_type = 'R' if "O2'" in atoms else 'D'

        # HETATM entires, positions 67-72 (1-indexed) are unused, so we fill our information in here
        for base_line in curr_base_lines:
          new_base_lines = list(base_line)
          new_base_lines[67] = nucleic_acid_type
          new_base_lines[68] = 'N'
          new_base_lines[69] = 'A'

          # if we have a phosphate, the atom corresponds to a backbone, so add 'B' to designate this:
          if "'" in base_line[12:16].strip() or 'P' in base_line[12:16].strip():
            new_base_lines[70] = 'B'

          # write out all the updated lines to the new temporary file
          tmp_nuchandle.write(''.join(new_base_lines))

      # reset the current list of lines for the bext base to process:
      curr_base_lines = []
      curr_base_id = lig_line[22:26].strip()

    # no matter what, add the current line (to be eventually processed)
    curr_base_lines.append(lig_line)
  orig_ligand_handle.close()

  # repeat the processing step for the very last set of entries:
  if len(curr_base_lines) > 0:
    atoms = [base_line[12:16].strip() for base_line in curr_base_lines]
    nucleic_acid_type = 'R' if "O2'" in atoms else 'D'

    for base_line in curr_base_lines:
      new_base_lines = list(base_line)
      new_base_lines[67] = nucleic_acid_type
      new_base_lines[68] = 'N'
      new_base_lines[69] = 'A'
      if "'" in base_line[12:16].strip() or 'P' in base_line[12:16].strip():
        new_base_lines[70] = 'B'
      tmp_nuchandle.write(''.join(new_base_lines))

  tmp_nuchandle.close()
  return tmp_ligand_file  # return the full path to the newly-updated NUC ligand file


########################################################################################################

def create_distlist_files(annotation_file, pdb_ids, receptor_pdb_dir, ligand_pdb_dir,
                          distlist_dir, annot_dir, include_backbone=False):
  """
  :param annotation_file: full path to a tab-delineated annotation file (as downloaded from BioLiP),
                          e.g., DBPATH+'biolip/annotations/BioLiP_allfiles-2015-09-26.txt'
  :param pdb_ids: subset of PDB IDs to create distlist files for
  :param receptor_pdb_dir: full path to a directory containing all receptor protein PDB structures
  :param ligand_pdb_dir: full path to a directory containing all ligand PDB structures (corresponding to proteins)
  :param distlist_dir: full path to a directory containing the processed "_distances.txt.gz" output files
  :param annot_dir: full path to a directory containing the processed "_annotation.txt.gz" output files
  :param include_backbone: boolean whether to include protein backbone->ligand proximities in calculations or not
  :return: None, but create a "distlist" file in the appropriate directory corresponding to each of the PDB 
           IDs specified, and update the annotation file accordingly.
  """

  # 3-letter to 1-letter abbreviation mapping for amino acids
  aaabbrv = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
             'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
             'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
             'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

  # van der waals radii to calculate what's "close enough"
  vdw = get_vdw_interaction_radii()

  # For each PDB ID that we are running on:
  pdb_id_index = 0  # keep track of our progress
  total_pdbs = str(len(pdb_ids))
  for current_pdb_id in pdb_ids:

    # print status message:
    pdb_id_index += 1
    sys.stderr.write('[' + str(pdb_id_index) + '/' + total_pdbs + '] ' + current_pdb_id + '\n')

    # find full directories to store output files (create it if necessary)
    prefix = current_pdb_id[:1] + '/' + current_pdb_id[:2] + '/'

    for outerdir in [distlist_dir, annot_dir]:
      for innerdir in ['', current_pdb_id[:1], current_pdb_id[:1] + '/' + current_pdb_id[:2]]:
        if not os.path.isdir(str(outerdir) + str(innerdir)):
          call(['mkdir', str(outerdir) + str(innerdir)])

    # two output files: (1) an updated annotation file, and (2) a "distlist" file containing pairwise Euc distances
    new_annot_file = annot_dir + prefix + current_pdb_id + '_annotation.txt'
    distlist_file = distlist_dir + prefix + current_pdb_id + '_distances.txt.gz'

    annot_outhandle = gzip.open(new_annot_file, 'w') if new_annot_file.endswith('gz') else open(new_annot_file, 'w')
    distlist_outhandle = gzip.open(distlist_file, 'w') if distlist_file.endswith('gz') else open(distlist_file, 'w')
    distlist_outhandle.write('\n'.join(['# All pairwise distances between receptor protein chain residue atoms ' +
                                        'and ligand atoms for ' + current_pdb_id,
                                        '# NOTE: columns 8-11 contain the overlap area between Gaussian ' +
                                        'distributions centered at each atom with ',
                                        '#   standard deviations set to either the van der Waals radii of ' +
                                        'the two atoms or to 1.5']) + '\n')
    distlist_outhandle.write('\t'.join(['#pdbID-pdbChain', 'receptor_aa_1-index', 'receptor_aa_value',
                                        'receptor_atom_id', 'receptor_atom_value', 'ligand_id', 'ligand_atom_value',
                                        'euclidean_distance',
                                        'overlap_vdw_radii', 'integral_error_vdw_radii',
                                        'overlap_1.5', 'integral_error_1.5',
                                        'full_receptor_sequence']) + '\n')

    receptor_sequences_output = set()  # keep track of which pdbID-pdbChain receptor sequences we have processed

    # process appropriate entries from the annotation file one-by-one
    annotation_inhandle = gzip.open(annotation_file) if annotation_file.endswith('gz') else open(annotation_file)
    for annotation in annotation_inhandle:

      # restrict to the subset of PDB IDs that we are trying to process:
      annot_pdb_id = annotation[:-1].split('\t')[0]
      if annot_pdb_id != current_pdb_id:
        continue

      # get the rest of the information from the annotation line:
      (pdb_chain,
       resolution,
       binding_site_id,
       ligand_id,
       ligand_chain,
       ligand_serial_number,
       binding_site_residue__pdb_index,
       binding_site_residue__one_index,
       catalytic_site_residue__pdb_index,
       catalytic_site_residue__one_index,
       enzyme_commission_number,
       go_terms,
       binding_affinity__manual,
       binding_affinity__binding_moad,
       binding_affinity__pdb_bind_cn,
       binding_affinity__binding_db,
       uniprot_id,
       pubmed_id,
       full_receptor_sequence) = annotation[:-1].split('\t')[1:]

      # attempt to find the receptor structure file and ligand structure file:
      receptor_structure_file = receptor_pdb_dir + annot_pdb_id + pdb_chain + '.pdb'
      ligand_structure_file = ligand_pdb_dir+annot_pdb_id+'_'+ligand_id+'_'+ligand_chain+'_'+ligand_serial_number+'.pdb'

      if not os.path.isfile(receptor_structure_file):
        sys.stderr.write('Could not open ' + receptor_structure_file + '\n')
        continue
      if not os.path.isfile(ligand_structure_file):
        sys.stderr.write('Could not open ' + ligand_structure_file + '\n')
        continue

      # --------------------------------------------------------------------------------
      # PROCESS RECEPTOR FILE

      positions = {}  # PDB index -> residue AA
      atoms = {}  # PDB index -> set((atom, X, Y, Z))

      current_aa_position = -100000  # keep track of which residue we are processing (multiple atoms per AA residue)

      receptor_handle = open(receptor_structure_file)
      for pdb_line in receptor_handle:

        # skip ahead to the atom entries (ignore hydrogens)
        if not pdb_line.startswith('ATOM') or pdb_line[76:78].strip() == 'H':
          continue

        residue_aa_value = aaabbrv.get(pdb_line[17:20].strip(), 'X')
        residue_aa_position = int(pdb_line[22:26].strip())
        loc_x = float(pdb_line[30:38].strip())
        loc_y = float(pdb_line[38:46].strip())
        loc_z = float(pdb_line[46:54].strip())
        ligand_atom_value = pdb_line[76:78].strip()

        # moved on to a new amino acid in the protein chain
        if residue_aa_position != current_aa_position:
          current_aa_position = residue_aa_position
          positions[current_aa_position] = residue_aa_value  # keep track of the amino acid value at each index
          atoms[current_aa_position] = set()

          # skip AA backbone atoms if specified (i.e., include just side chain proximity)
          if not include_backbone:
            for _ in xrange(3):
              try:
                next(receptor_handle)
              except StopIteration:
                continue
          continue

        # at each intermediate atom per residue, store all the location information
        atoms[current_aa_position].add((loc_x, loc_y, loc_z, ligand_atom_value))
      receptor_handle.close()

      # generate, for each residue ID, a sorted list of atoms (x, y, x, atom_value, atom_index):
      for residue_aa_index in atoms.keys():
        atoms[residue_aa_index] = sorted(list(atoms[residue_aa_index]))  # sort the atoms by x,y,z location...
        for atom_index in xrange(len(atoms[residue_aa_index])):  # add an index to each atom
          atoms[residue_aa_index][atom_index] = tuple(list(atoms[residue_aa_index][atom_index]) + [atom_index + 1])

      # --------------------------------------------------------------------------------
      # PROCESS LIGAND FILE

      # if the ligand is a nucleic acid, subtype it to RNA, DNA, backbone, or base:
      if ligand_id == 'NUC':
        ligand_structure_file = subtype_nuc_ligand(ligand_structure_file)

        # keep track of separate sets of binding residues corresponding to the four new types of ligands
        new_binding_residues = {'DNA': set(), 'DNAB': set(), 'RNA': set(), 'RNAB': set()}
      else:
        new_binding_residues = {'': set()}

      # start to process at the first receptor amino acid position that we have information for:
      start_index = int(min(map(int, positions.keys())))

      ligand_handle = open(ligand_structure_file)
      for pdb_line in ligand_handle:

        # again, skip forward to the atom entries and ignore hydrogens in distance calculations
        if not pdb_line.startswith('HETATM') or pdb_line[76:78].strip() == 'H':
          continue

        loc_x = float(pdb_line[30:38].strip())
        loc_y = float(pdb_line[38:46].strip())
        loc_z = float(pdb_line[46:54].strip())
        ligand_sub_type = pdb_line[66:72].strip()  # this entry was added into the temporary file processed
        ligand_atom_value = pdb_line[76:78].strip()

        # recall that receptor_atoms is now a sorted list of tuples of (x, y, z, atom_value, atom_index)
        for receptor_aa_position, receptor_atoms in atoms.items():
          for (ratom_x, ratom_y, ratom_z, ratom_value, ratom_index) in receptor_atoms:

            # calculate the distance between this ligand atom and each receptor atom
            current_euclidean_distance = get_euclidean_distance((loc_x, loc_y, loc_z), (ratom_x, ratom_y, ratom_z))

            # keep track of the distance (within reason, DISTANCE_CUTOFF specified above)
            if current_euclidean_distance <= DISTANCE_CUTOFF:
              distlist_outhandle.write('\t'.join(map(str, [annot_pdb_id + pdb_chain,
                                                           receptor_aa_position - start_index + 1,
                                                           positions[receptor_aa_position],
                                                           str(ratom_index) + '/' + str(len(receptor_atoms)),
                                                           ratom_value,
                                                           ligand_id + ligand_sub_type,
                                                           ligand_atom_value,
                                                           current_euclidean_distance,
                                                           'N/A', 'N/A', 'N/A', 'N/A'])))

              # to save harddrive space in the resulting files, only write out the full receptor sequence ONCE
              if annot_pdb_id + pdb_chain not in receptor_sequences_output:
                receptor_sequences_output.add(annot_pdb_id + pdb_chain)
                distlist_outhandle.write('\t' + ''.join([positions[aa_index] if aa_index in positions else 'X'
                                                         for aa_index in
                                                         xrange(start_index,
                                                                int(max(map(int, positions.keys()))) + 1)]) +
                                         '\n')
              else:
                distlist_outhandle.write('\t\n')

            # call this receptor_aa_position a binding residue IFF the distance is close enough
            if current_euclidean_distance <= 2 + vdw.get(ligand_atom_value, 1.5) + vdw.get(ratom_value, 1.5):
              new_binding_residues[ligand_sub_type].add(receptor_aa_position)

      ligand_handle.close()

      # remove the temporary file we had created, just in case
      if '-tmpnuc' in ligand_structure_file:
        call(['rm', ligand_structure_file])

      # the very last step is to update the annotation file for this particular PDB ID:
      for ligand_sub_type in new_binding_residues.keys():
        if len(new_binding_residues[ligand_sub_type]) > 0:
          annot_outhandle.write('\t'.join([annot_pdb_id,
                                           pdb_chain,
                                           resolution,
                                           binding_site_id,
                                           ligand_id + ligand_sub_type,
                                           ligand_chain,
                                           ligand_serial_number,
                                           ' '.join([positions[aa_index] + str(aa_index)
                                                     for aa_index in
                                                     sorted(list(new_binding_residues[ligand_sub_type]))]),
                                           ' '.join([positions[aa_index] + str(aa_index - start_index + 1)
                                                     for aa_index in
                                                     sorted(list(new_binding_residues[ligand_sub_type]))]),
                                           catalytic_site_residue__pdb_index,
                                           catalytic_site_residue__one_index,
                                           enzyme_commission_number,
                                           go_terms,
                                           binding_affinity__manual,
                                           binding_affinity__binding_moad,
                                           binding_affinity__pdb_bind_cn,
                                           binding_affinity__binding_db,
                                           uniprot_id,
                                           pubmed_id,
                                           ''.join([positions[aa_index] if aa_index in positions else 'X'
                                                    for aa_index in
                                                    xrange(start_index, int(max(map(int, positions.keys()))) + 1)])]) +
                                '\n')
    annotation_inhandle.close()  # close the input annotation file (done processing each line)

    distlist_outhandle.close()  # close the two new output files that were created
    annot_outhandle.close()


########################################################################################################

def update_distlist_files(distlist_infiles):
  """
  :param distlist_infiles: ordered list of full paths to "distlist" files to add extra info to:
  :return: None, create a new tab-delineated file with the headings:
           [0] pdbID-pdbChain, 
           [1] residue number (1-indexed), 
           [2] residue value, 
           [3] atomID (#/total), 
           [4] ligandID, 
           [5] ligandAtomValue, 
           [6] euclidean distance, 
           [7] VDW INTEGRAL AREA, 
           [8] VDW INTEGRAL ERROR, 
           [9] STD INTEGRAL AREA, 
           [10] STD INTEGRAL ERROR, 
           [11] receptor sequence
  """

  # input all the van der Waals radii that we know about
  vdw = get_vdw_interaction_radii()

  # Keep track of all (distance, set of radii) that we calculated the integral and error for,
  # so that we don't redo it ! MEMORY INTENSIVE but less time intensive (hopefully)
  for distlist_filename in distlist_infiles:
    sys.stderr.write('Adding integrals to ' + distlist_filename + '\n')

    # open a new R file to calculate the integrals (fastest option)
    rscript_integral_results = DATAPATH + random_filename() + '.txt'

    rscript_file = DATAPATH + random_filename() + '.R'
    rscript_handle = open(rscript_file, 'w')
    rscript_handle.write('sink("' + rscript_integral_results + '")\n')
    rscript_handle.write("""min.f1f2 <- function(x, mu1, mu2, sd1, sd2) {
                              f1 <- dnorm(x, mean=mu1, sd=sd1)
                              f2 <- dnorm(x, mean=mu2, sd=sd2)
                              pmin(f1,f2)
                         }\n""")

    lines = []  # keep track of the original lines from the input file (to be updated)
    integrals = {}  # we only need to keep track of the integrals between specific distances (avoid redundancies)

    distlist_inhandle = gzip.open(distlist_filename) if distlist_filename.endswith('gz') else open(distlist_filename)
    for dline in distlist_inhandle:

      if dline.startswith('#'):
        lines.append([dline])
        continue

      lines.append(dline[:-1].split('\t'))  # keep track of the original line
      receptor_atom, _, ligand_atom, euclidean_distance = dline[:-1].split('\t')[4:8]

      if float(euclidean_distance) <= DISTANCE_CUTOFF:
        sd1 = vdw.get(receptor_atom, '1.5')
        sd2 = vdw.get(ligand_atom, '1.5')

        # calculate the integrals if we haven't already
        for rsd1, rsd2 in sorted([[sd1, sd2], ['1.5', '1.5']]):
          if (euclidean_distance, (rsd1, rsd2)) not in integrals:
            rscript_handle.write('writeLines(c("' + euclidean_distance + '\t' + rsd1 + '\t' + rsd2 + '"))\n')
            rscript_handle.write("""integrate(min.f1f2, -Inf, Inf, mu1 = 0, mu2 = """ + str(euclidean_distance) + """, 
                                         sd1 = """ + str(rsd1) + """, sd2 = """ + str(rsd2) + """)\n""")
    distlist_inhandle.close()

    # Run for a single complete file at a time (fewer calls)
    rscript_handle.write('sink()\n')
    rscript_handle.close()
    call(['Rscript', rscript_file])
    call(['rm', rscript_file])

    if rscript_integral_results.endswith('gz'):
      integrals_inhandle = gzip.open(rscript_integral_results)
    else:
      integrals_inhandle = open(rscript_integral_results)

    for int_line in integrals_inhandle:
      if 'with absolute error' not in int_line:
        euc_dist, sd1, sd2 = int_line.strip().split('\t')
        overlap = integrals_inhandle.next().strip()
        integrals[(euc_dist, tuple(sorted([sd1, sd2])))] = (overlap.split()[0], overlap.split()[-1])
    integrals_inhandle.close()
    call(['rm', rscript_integral_results])

    # finally, write out the new results (i.e., original files with added columns)
    processed_filename = distlist_filename.replace('_distances', '_distances-woverlap')

    if processed_filename.endswith('gz'):
      processed_outhandle = gzip.open(processed_filename, 'w')
    else:
      processed_outhandle = open(processed_filename, 'w')

    # Write out each of the original lines, in order, with new features:
    for origline in lines:
      if origline[0].startswith('#'):
        processed_outhandle.write('\t'.join(origline))
        continue

      processed_outhandle.write('\t'.join(origline[:8]) + '\t')
      if float(origline[7]) <= DISTANCE_CUTOFF:
        processed_outhandle.write('\t'.join([integrals[(origline[7],
                                                        tuple(sorted([vdw.get(origline[4], '1.5'),
                                                                      vdw.get(origline[6], '1.5')])))][0],
                                             integrals[(origline[7],
                                                        tuple(sorted([vdw.get(origline[4], '1.5'),
                                                                      vdw.get(origline[6], '1.5')])))][1],
                                             integrals[(origline[7], tuple(sorted(['1.5', '1.5'])))][0],
                                             integrals[(origline[7], tuple(sorted(['1.5', '1.5'])))][1]]))
      else:
        processed_outhandle.write('\t'.join(['N/A', 'N/A', 'N/A', 'N/A']))  # we didn't calculate integrals here
      processed_outhandle.write('\t' + '\t'.join(origline[12:]) + '\n')  # sequence is written, regardless
    processed_outhandle.close()

    # overwrite the original file IFF we finished without an error
    call(['mv', processed_filename, distlist_filename])


########################################################################################################

if __name__ == "__main__":

  # Parse the command-line arguments
  parser = argparse.ArgumentParser(description='Process the BioLiP receptor/ligand PDB files to calculate ' +
                                               'distances between receptor and ligand atoms.')

  parser.add_argument('--prefix', type=str,
                      help='Two letter prefix to subset of PDB IDs to process',
                      default='10')
  parser.add_argument('--function', type=str,
                      help='Step in the pipeline to run: e.g., create_distlist, update_distlist, create_fasta',
                      default='create_distances',
                      choices={'create_distances', 'update_distances'})

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
  if args.function == 'create_distances':
    """
    Calculate the distances between the receptor protein chains and their ligands to create files
    named "_distances.txt" in the biolip/distances/ directory
    """

    receptor_directory = DATAPATH + 'biolip/downloaded_data/receptor/'
    ligand_directory = DATAPATH + 'biolip/downloaded_data/ligand/'
    distance_out_directory = DATAPATH + 'biolip/processed_data/distances/'
    annotations_out_directory = DATAPATH + 'biolip/processed_data/annotations/'

    # create the directories we need (if we need them)
    for maindir in [distance_out_directory, annotations_out_directory]:
      for subdir in ['', args.prefix[0], args.prefix[0] + '/' + args.prefix]:
        if not os.path.isdir(maindir + subdir):
          call(['mkdir', maindir + subdir])

    if not os.path.isdir(receptor_directory) or not os.path.isdir(ligand_directory):
      sys.stderr.write('No such directories:\n' + receptor_directory + '\n' + ligand_directory + '\n')
      sys.exit(1)

    # create initial distance files:
    create_distlist_files(current_annotation_file, subset_pdb_ids, receptor_directory, ligand_directory,
                          distance_out_directory, annotations_out_directory)

  # ----------------------------------------------------------------------------------------------------
  elif args.function == 'update_distances':
    """
    Use the previously calculated distances between the receptor protein chains and their ligands
    to create files with new distance metrics that require Gaussian overlap calculations
    """

    distance_out_directory = DATAPATH + 'biolip/processed_data/distances/' + args.prefix[0] + '/' + args.prefix + '/'
    if not os.path.isdir(distance_out_directory):
      sys.stderr.write('No such directory: ' + distance_out_directory + '\n')
      sys.exit(1)

    distlist_files = sorted([distance_out_directory + biolip_file for biolip_file in os.listdir(distance_out_directory)
                             if biolip_file.endswith('_distances.txt.gz') and
                             biolip_file.split('_')[0] in subset_pdb_ids])

    update_distlist_files(distlist_files)
