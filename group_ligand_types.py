#!/usr/bin/python

"""
Download ligand information from the Chemical Components Dictionary, DrugBank, and the Human Metabolome
Database and compare molecule SMILES strings using the Tanimoto coefficient to group ligands into 
supergroups "ION_", "DRUGLIKE_", and "METABOLITE_"

Contact Shilpa Nadimpalli Kobren (snadimpa@princeton.edu) with questions
"""

import os
import sys
import gzip
import re
import argparse
import xml.etree.ElementTree as element_tree
from subprocess import call
from string import ascii_letters, digits
from random import choice


########################################################################################################

# path to where this script is currently located (and to where all data should be stored) -- this can
# be updated
DATAPATH = os.path.dirname(os.path.abspath(__file__))+'/downloaded_data/'

TANIMOTO_CUTOFF = 0.9  # SMILES strings at or above this cutoff are considered the same molecule


########################################################################################################
# DOWNLOAD AND PARSE ORIGINAL INFORMATION
########################################################################################################

def download_and_parse_ccd(ccd_path=DATAPATH+'ccd/',
                           outfile=DATAPATH+'ccd/chemical_components_dictionary-parsed.tsv'):
  """
  :param ccd_path: full path to directory where downloaded files and parsed file should be stored
  :param outfile: full path to a file to write output to
  :return: None, but print a success message upon completion; parse and combine files from
           wwPDB's Chemical Components Dictionary
  """

  # download required file
  if not os.path.isdir(ccd_path):
    call(['mkdir', ccd_path])

  components_file = ccd_path + 'components.cif.gz'

  if not os.path.isfile(components_file):
    call(['wget', 'ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz', '-O', components_file])

  # parse the downloaded file
  mmcif_to_molecular_info = {}  # ligand_id -> "Information Name" -> "Information Value"
  mmcif_id = 'xxxxxxx'  # starting identifier (to be replaced)

  components_handle = gzip.open(components_file) if components_file.endswith('.gz') else open(components_file)
  for comp_line in components_handle:

    # extract the mmCIF ligand ID
    if comp_line.startswith('_chem_comp.id '):
      mmcif_id = comp_line.strip().split()[-1]
      mmcif_to_molecular_info[mmcif_id] = {}

    # extract the primary name of the ligand
    elif comp_line.startswith('_chem_comp.name '):
      if '"' in comp_line:
        ligand_name = comp_line[comp_line.find('"') + 1:].strip()[:-1]
      elif len(comp_line.strip().split()) > 1:
        ligand_name = comp_line.strip().split()[-1]
      else:
        ligand_name = '?'
      mmcif_to_molecular_info[mmcif_id]['ligand_name'] = ligand_name

    # store all synonym names for the ligand
    elif comp_line.startswith('_chem_comp.pdbx_synonyms '):
      if '"' in comp_line:
        synonyms = comp_line[comp_line.find('"') + 1:].strip()[:-1]
      elif len(comp_line.strip().split()) > 1:
        synonyms = comp_line.strip().split()[-1]
      else:
        synonyms = '?'
      mmcif_to_molecular_info[mmcif_id]['synonyms'] = synonyms if synonyms != 'None' else '?'

    # extract the primary SMILES string (this likely matches the canonical SMILES string)
    elif comp_line.startswith(mmcif_id + ' SMILES '):
      if 'SMILES' not in mmcif_to_molecular_info[mmcif_id]:
        mmcif_to_molecular_info[mmcif_id]['SMILES'] = set()

      smiles_value = (' '.join(comp_line.strip().split()[2:-2]).replace('"', ''),
                      comp_line.strip().split()[-2],
                      comp_line.strip().split()[-1].replace('"', ''))

      mmcif_to_molecular_info[mmcif_id]['SMILES'].add(smiles_value)

    # extract the canonical SMILES string
    elif comp_line.startswith(mmcif_id + ' SMILES_CANONICAL '):
      if 'canonical_SMILES' not in mmcif_to_molecular_info[mmcif_id]:
        mmcif_to_molecular_info[mmcif_id]['canonical_SMILES'] = set()

      canonical_smiles_value = (' '.join(comp_line.strip().split()[2:-2]).replace('"', ''),
                                comp_line.strip().split()[-2],
                                comp_line.strip().split()[-1].replace('"', ''))

      mmcif_to_molecular_info[mmcif_id]['canonical_SMILES'].add(canonical_smiles_value)

    # extract the InChI string
    elif comp_line.startswith(mmcif_id + ' InChI '):
      no_truncation = len(comp_line.strip().split()) > 4  # keep track of whether the line has been truncated

      inchi_value = (' '.join(comp_line.strip().split()[2:-2 if no_truncation else -1]).replace('"', ''),
                     comp_line.strip().split()[-2 if no_truncation else -1],
                     comp_line.strip().split()[-1].replace('"', '').replace('InChI=', '') if no_truncation else '')

      mmcif_to_molecular_info[mmcif_id]['InChI'] = inchi_value

    # extract the InChI Key
    elif comp_line.startswith(mmcif_id + ' InChIKey '):
      no_truncation = len(comp_line.strip().split()) > 4

      inchi_key_value = (' '.join(comp_line.strip().split()[2:-2 if no_truncation else -1]).replace('"', ''),
                         comp_line.strip().split()[-2 if no_truncation else -1],
                         comp_line.strip().split()[-1].replace('"', '') if no_truncation else '')

      mmcif_to_molecular_info[mmcif_id]['InChIkey'] = inchi_key_value
  components_handle.close()

  # write all the parsed information to the specified tab-delimited file
  out_handle = gzip.open(outfile, 'w') if outfile.endswith('gz') else open(outfile, 'w')
  out_handle.write('\n'.join(['# All ligand information was downloaded from wwPDB\'s Chemical Components Dictionary:',
                              '# ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz',
                              '\t'.join(['#ligand_mmcif_id', 'ligand_name', 'synonyms', 'InChI_key', 'InChI',
                                         'SMILES', 'canonical_SMILES'])])+'\n')

  for ligand_mmcif_id in sorted(mmcif_to_molecular_info.keys()):
    ligand_info = mmcif_to_molecular_info[ligand_mmcif_id]
    out_handle.write('\t'.join([ligand_mmcif_id, ligand_info.get('ligand_name', ''),
                                ligand_info.get('synonyms', '?'),
                                ':'.join(ligand_info.get('InChIkey', set())),
                                ':'.join(ligand_info.get('InChI', set())),
                                ','.join([':'.join(smiles_value) for smiles_value in
                                          sorted(list(ligand_info.get('SMILES', set())))]),
                                ','.join([':'.join(smiles_value) for smiles_value in
                                          sorted(list(ligand_info.get('canonical_SMILES', set())))])]) + '\n')
  out_handle.close()

  sys.stderr.write('Parsed all molecular information for ' + "{:,}".format(len(mmcif_to_molecular_info.keys())) +
                   ' ligands into '+outfile+'\n')


########################################################################################################

def download_and_parse_hmdb(hmdb_path=DATAPATH+'hmdb/',
                            outfile=DATAPATH+'hmdb/human_metabolome_database-parsed.tsv'):
  """
  :param hmdb_path: full path to directory where downloaded files and parsed file should be stored
  :param outfile: full path to a file to write output to
  :return: None, but print a success message upon completion; parse files from the Human Metabolome Database (HMDB)
  """

  # download required file
  if not os.path.isdir(hmdb_path):
    call(['mkdir', hmdb_path])

  metabolites_file = hmdb_path + 'hmdb_metabolites.zip'

  if not os.path.isfile(metabolites_file):
    call(['wget', 'http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip', '-O', metabolites_file])
    call(['unzip', metabolites_file, '-d', hmdb_path[:-1]])

  # write to the output handle as we parse the files, 1 by 1:
  out_handle = gzip.open(outfile, 'w') if outfile.endswith('gz') else open(outfile, 'w')
  out_handle.write('\n'.join(['# All information downloaded from the Human Metabolome Database (HMDB), version 4.0:',
                              '# http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip',
                              '\t'.join(['#hmdb_id', 'metabolite_name', 'category', 'synonyms', 'iupac_name',
                                         'SMILES', 'InChI', 'InChIKey'])])+'\n')

  files_processed = 0

  # create an element tree parser:
  prefix = ''
  for event, e in element_tree.iterparse(hmdb_path + 'hmdb_metabolites.xml',
                                         events=('start', 'end', 'start-ns', 'end-ns')):
    if event == 'start-ns':
      prefix = '{'+e[1]+'}'
    if event == 'end':
      if e.tag == prefix+'metabolite':

        # NOTE: we have to try/except because some entries have Unicode characters that are not included in
        # basic ASCII (and cannot be written out) -- THIS is why we write directly instead of storing

        try:  # hmdb_id
          hmdb_id = ','.join(str(child.text).strip() for child in e.findall(prefix+'accession')[:1])
          out_handle.write(hmdb_id)
        except:
          out_handle.write('')

        name = ''
        try:  # metabolite_name
          name = ','.join(str(child.text).strip() for child in e.findall(prefix+'name')[:1])
          out_handle.write('\t' + name)
        except:
          out_handle.write('\t')

        categories = []
        try:  # category
          for c2 in e:
            if c2.tag == prefix+'taxonomy':
              for child in c2:
                if 'class' in child.tag:
                  categories.append(child.text)

          for child in e:
            if child.tag == prefix + 'ontology':
              for c2 in child:
                if c2.tag == prefix + 'root':
                  found_disposition = False
                  found_source = False
                  for c3 in c2:
                    if c3.tag == prefix + 'term' and c3.text == 'Disposition':
                      found_disposition = True
                    if found_disposition and c3.tag == prefix + 'descendants':
                      for c4 in c3:  # these are all "descendant"
                        for c5 in c4:
                          if c5.tag == prefix + 'term' and c5.text == 'Source':
                            found_source = True
                          if found_source and c5.tag == prefix + 'descendants':
                            for c6 in c5:  # these are all "descendant"
                              for c7 in c6:
                                if c7.tag == prefix + 'term':
                                  categories.append(c7.text)
          out_handle.write('\t' + ','.join(sorted(list(set([str(a).strip().lower() for a in categories])))))
        except:
          out_handle.write('\t')
          i = 0
          for l in sorted(list(set([a for a in categories]))):
            try:
              out_handle.write((',' if i > 0 else '') + str(l).strip().lower())
              i += 1
            except:
              pass

        synonyms = []
        try:  # synonyms
          for c2 in e:
            if c2.tag == prefix + 'synonyms':
              synonyms += [str(child.text).strip() for child in c2]
          out_handle.write('\t' + ','.join(sorted(list(set([a for a in synonyms if a != name])))))
        except:
          out_handle.write('\t')
          i = 0
          for l in sorted(list(set([a for a in synonyms if a != name]))):
            try:
              out_handle.write((',' if i > 0 else '') + l)
              i += 1
            except:
              pass

        try:  # iupac name
          iupac = ','.join(str(child.text).strip() for child in e.findall(prefix+'iupac_name')[:1])
          out_handle.write('\t' + iupac)
        except:
          out_handle.write('\t')

        try:  # smiles string
          smiles = ','.join(str(child.text).strip() for child in e.findall(prefix+'smiles')[:1])
          out_handle.write('\t' + smiles)
        except:
          out_handle.write('\t')

        try:  # inchi
          inchi = ','.join(str(child.text).strip() for child in e.findall(prefix+'inchi')[:1]).replace('InChI=', '')
          out_handle.write('\t' + inchi)
        except:
          out_handle.write('\t')

        try:  # inchi key
          inchikey = ','.join(str(child.text).strip() for child in
                              e.findall(prefix+'inchikey')[:1]).replace('InChIKey=', '')
          out_handle.write('\t' + inchikey)
        except:
          out_handle.write('\t')

        out_handle.write('\n')
        files_processed += 1

  out_handle.close()
  sys.stderr.write('Processed information for '+"{:,}".format(files_processed)+' metabolites into '+outfile+'\n')
  sys.exit(1)


########################################################################################################

def download_and_parse_drugbank(drugbank_path=DATAPATH + 'drugbank/',
                                outfile=DATAPATH + 'drugbank/drugbank-parsed.tsv'):
  """
  :param drugbank_path: full path to directory where downloaded files and parsed file should be stored
  :param outfile: full path to a file to write output to
  :return: None, but print a success message upon completion; parse files from DrugBank
  """

  # download required file
  if not os.path.isdir(drugbank_path):
    call(['mkdir', drugbank_path])

  drug_file = drugbank_path + 'drugbank_all_full_database.xml'

  if not os.path.isfile(drug_file):
    sys.stderr.write('Please download the most recent DrugBank file from https://www.drugbank.ca/releases/latest\n')
    sys.stderr.write('and store it locally as '+drug_file+'\n')
    sys.exit(1)

  # write out header for eventual parsed DrugBank file
  out_handle = gzip.open(outfile, 'w') if outfile.endswith('gz') else open(outfile, 'w')
  out_handle.write('\n'.join(['# All information downloaded from DrugBank, Version 5.0:',
                              '# https://www.drugbank.ca/releases/latest',
                              '\t'.join(['#drugbank_id', 'drug_name', 'synonyms', 'SMILES', 'InChI', 'InChIKey',
                                         'other_ids'])])+'\n')

  # start to parse the drugbank XML file:
  e = element_tree.parse(drug_file).getroot()
  prefix = e.tag[e.tag.find('{'):e.tag.find('}') + 1]

  processed_drugs = 0
  for drug in e.findall(prefix + 'drug'):

    smiles, inchi, inchikey, = '', '', ''
    ids = set()

    if drug.find(prefix + 'calculated-properties') is not None:
      for prop in drug.find(prefix + 'calculated-properties').findall(prefix + 'property'):
        if prop.find(prefix + 'kind').text == 'SMILES':
          smiles = re.sub(' +', ' ', prop.find(prefix + 'value').text.strip())
        elif prop.find(prefix + 'kind').text == 'InChI':
          inchi = re.sub(' +', ' ', prop.find(prefix + 'value').text.replace('InChI=', '').strip())
        elif prop.find(prefix + 'kind').text == 'InChIKey':
          inchikey = re.sub(' +', ' ', prop.find(prefix + 'value').text.strip())

    if drug.find(prefix + 'external-identifiers') is not None:

      for prop in drug.find(prefix + 'external-identifiers').findall(prefix + 'external-identifier'):
        ids.add(
          re.sub(' +', ' ', prop.find(prefix + 'resource').text + ':' + prop.find(prefix + 'identifier').text.strip()))

    if drug.find(prefix + 'products') is not None:
      products = [re.sub(' +', ' ', products.find(prefix + 'name').text.strip()) for products in
                  drug.find(prefix + 'products').findall(prefix + 'product') if
                  products is not None and products.find(prefix + 'name') is not None]
    else:
      products = []

    if drug.find(prefix + 'synonyms') is not None:
      synonyms = [re.sub(' +', ' ', syn.text.strip()) for syn in drug.find(prefix + 'synonyms').findall('synonym') if
                  syn is not None]
    else:
      synonyms = []

    drugbankid = re.sub(' +', ' ', drug.find(prefix + 'drugbank-id').text.strip()) if drug.find(
      prefix + 'drugbank-id') is not None else ''
    drugname = re.sub(' +', ' ', drug.find(prefix + 'name').text.strip()) if drug.find(
      prefix + 'name') is not None else ''

    # We have to do this because some drug names have unicode characters that are not in
    # basic ASCII and therefore cannot be written out. Oops. Should be irrelevant for us, really...
    try:
      out_handle.write(drugbankid)
    except:
      out_handle.write('')

    try:
      out_handle.write('\t' + drugname)
    except:
      out_handle.write('\t')

    try:
      out_handle.write('\t' + ','.join(sorted(list(set(synonyms + products)))))
    except:
      out_handle.write('\t')
      i = 0
      for l in sorted(list(set(synonyms + products))):
        try:
          out_handle.write((',' if i > 0 else '') + l)
          i += 1
        except:
          pass

    try:
      out_handle.write('\t' + smiles)
    except:
      out_handle.write('\t')

    try:
      out_handle.write('\t' + inchi)
    except:
      out_handle.write('\t')

    try:
      out_handle.write('\t' + inchikey)
    except:
      out_handle.write('\t')

    try:
      out_handle.write('\t' + ','.join(sorted(list(ids))))
    except:
      out_handle.write('\t')
      i = 0
      for l in sorted(list(ids)):
        try:
          out_handle.write((',' if i > 0 else '') + l)
          i += 1
        except:
          pass
    out_handle.write('\n')
    processed_drugs += 1

  out_handle.close()
  sys.stderr.write('Processed information for '+"{:,}".format(processed_drugs)+' drugs into '+outfile+'\n')


########################################################################################################
# CALCULATE TANIMOTO SMILES->SMILES COMPARISONS
########################################################################################################

def random_filename(size=10, chars=ascii_letters + digits):
  """
  :param size: length of random string to return
  :param chars: set of characters from which to randomly select (with replacement)
  :return: a random string of length size comprised of chars (lowercase ASCII letters and digits 0-9)
  """

  return ''.join(choice(chars) for _ in range(size))


########################################################################################################

def process_alternate_ligands(alt_ligand_file, smiles_string_index, all_alt_smiles_file, all_alt_names_file):
  """
  :param alt_ligand_file: full path to a tab-delimited file where one of the column values is a
                                molecular structure SMILES string
  :param smiles_string_index: the 0-indexed column of the SMILES string in alternate_ligand_file
  :param all_alt_smiles_file: full path to an intermediate file required by "babel" containing just SMILES strings
  :param all_alt_names_file: full path to an intermediate file containing an ordered list of ligands in the smiles file
  :return: none, but print success message upon successful write of files
  """

  # Go through the alternate ligand->SMILES file, and keep an ordered list of molecular IDs, names, and SMILES strings
  alt_ligands = []

  alt_ligand_handle = gzip.open(alt_ligand_file) if alt_ligand_file.endswith('gz') else open(alt_ligand_file)
  for ligand_line in alt_ligand_handle:

    # make sure that this line may contain a SMILES string
    if ligand_line.startswith('#') or len(ligand_line[:-1].split('\t')) < smiles_string_index + 1:
      continue

    alt_smiles = ligand_line[:-1].split('\t')[smiles_string_index]

    if len(alt_smiles) > 0:
      alt_ligand_id, alt_ligand_name = ligand_line[:-1].split('\t')[:2]
      alt_ligands.append((alt_ligand_id, alt_ligand_name, str(alt_smiles)))

  alt_ligand_handle.close()

  # Write out *just the SMILES strings* to the temporary file (in the same order as we stored the ligands)
  if all_alt_smiles_file.endswith('gz'):
    temp_smiles_handle = gzip.open(all_alt_smiles_file, 'w')
  else:
    temp_smiles_handle = open(all_alt_smiles_file, 'w')
  if all_alt_names_file.endswith('gz'):
    temp_names_handle = gzip.open(all_alt_names_file, 'w')
  else:
    temp_names_handle = open(all_alt_names_file, 'w')

  for (alt_ligand_id, alt_ligand_name, alt_smiles) in alt_ligands:
    temp_smiles_handle.write(alt_smiles + '\n')
    temp_names_handle.write('\t'.join([alt_ligand_id, alt_ligand_name, alt_smiles])+'\n')
  temp_smiles_handle.close()
  temp_names_handle.close()

  sys.stderr.write('Wrote all SMILES strings to '+all_alt_smiles_file+'\n')
  sys.stderr.write('Wrote all ligand info (in the same order) to '+all_alt_names_file+'\n')


########################################################################################################

def calculate_tanimoto(alt_ligand_file, all_alt_smiles_file, all_alt_names_file, output_file,
                       orig_ligand_file=DATAPATH+'ccd/chemical_components_dictionary-parsed.tsv',
                       orig_ligand_smiles_index=6, start=1, end=27375):
  """
  :param alt_ligand_file: full path to a tab-delimited file where one of the column values is a
                                molecular structure SMILES string
  :param all_alt_smiles_file: full path to an intermediate file required by "babel" containing just SMILES strings
  :param all_alt_names_file: full path to an intermediate file containing ligand info in the same order as SMILES file
  :param output_file: full path to the output file to write the tab-delimited Tanimoto comparison results
  :param orig_ligand_file: full path to a "ligand" file labeled by mmCIF ID
  :param orig_ligand_smiles_index: the 0-indexed column in orig_ligand_file corresponding to the canonical SMILES string
  :param start: index of original ligand to start with
  :param end: index of original ligand to end on (1 more than actual 0-index, for instance 0:5 includes 0, 1, 2, 3, 4
  :return: None, but print success message
  """

  # Go through the alternate ligand->SMILES file, and keep an ordered list of molecular IDs, names, and SMILES strings
  alt_names_handle = gzip.open(all_alt_names_file) if all_alt_names_file.endswith('gz') else open(all_alt_names_file)
  alt_ligands = [lig_info.strip().split('\t')[:3] for lig_info in alt_names_handle
                 if len(lig_info.strip().split('\t')) > 2]
  alt_names_handle.close()

  # open and begin writing to output handle (header includes all file names and babel call)
  out_handle = gzip.open(output_file, 'w') if output_file.endswith('gz') else open(output_file)
  out_handle.write('\n'.join(['# All-against-all Tanimoto coefficients calculated between SMILES strings found in:',
                              '# (orig) '+orig_ligand_file,
                              '# (alt) '+alt_ligand_file,
                              '# by running from the command-line:',
                              '# babel <file with 1 orig SMILES string> <file with all alt SMILES strings> -ofpt ' +
                              '<output file>',
                              '\t'.join(['#orig_ligand_id', 'orig_ligand_name', 'orig_ligand_smiles',
                                         'alt_ligand_id', 'alt_ligand_name', 'alt_ligand_smiles',
                                         'tanimoto_coefficient'])])+'\n')

  # Now, go through the original list of possible ligands, and calculate the pairwise Tanimoto coefficients
  # between those SMILES strings and the alternate SMILES strings

  orig_ligand_index = 1
  orig_ligand_handle = gzip.open(orig_ligand_file) if orig_ligand_file.endswith('gz') else open(orig_ligand_file)
  for ligand_line in orig_ligand_handle:

    # make sure that this line may contain a SMILES string
    if ligand_line.startswith('#') or len(ligand_line[:-1].split('\t')) < orig_ligand_smiles_index + 1:
      continue

    # make sure that our ligand is in the proper range!
    if orig_ligand_index in range(start, end+1):

      ligand_id, ligand_name = ligand_line[:-1].split('\t')[:2]

      # the entry for canonical SMILES string may contain multiple:
      if len(ligand_line[:-1].split('\t')[orig_ligand_smiles_index]) < 1:
        ligand_smiles = ''
      else:
        ligand_smiles = ligand_line[:-1].split('\t')[orig_ligand_smiles_index].split(',')[0].split(':')[-1]

      if len(ligand_smiles) > 0:

        # print out a temporary file containing just this SMILES string
        current_orig_smiles_file = '/tmp/' + random_filename() + '.smi'
        orig_smiles_handle = open(current_orig_smiles_file, 'w')
        orig_smiles_handle.write(ligand_smiles + '\n')
        orig_smiles_handle.close()

        # and store the output from babel including the Tanimoto coefficient
        tanimoto_file = '/tmp/' + random_filename() + '.out'
        call(' '.join(['babel', current_orig_smiles_file, all_alt_smiles_file,
                       '-ofpt', tanimoto_file, '2>', '/dev/null']),
             shell=True)
        call(['rm', current_orig_smiles_file])

        # parse the babel output
        tanimoto_index = 0
        for taminoto_line in open(tanimoto_file):
          if 'Tanimoto' in taminoto_line:
            if float(taminoto_line.strip().split()[-1]) > 0:

              alt_ligand_id, alt_ligand_name, alt_ligand_smiles = alt_ligands[tanimoto_index][:3]

              out_handle.write('\t'.join([ligand_id, ligand_name, ligand_smiles,
                                          alt_ligand_id, alt_ligand_name, alt_ligand_smiles,
                                          taminoto_line.strip().split()[-1]]) + '\n')
            tanimoto_index += 1

        call(['rm', tanimoto_file])
        sys.stderr.write(str(orig_ligand_index) + '/27375\n')  # print progress
    orig_ligand_index += 1

    if orig_ligand_index > end:
      break

  orig_ligand_handle.close()

  out_handle.close()
  sys.stderr.write('Successfully calculated all-against-all pairwise Tanimoto coefficients between ' +
                   "{:,}".format(orig_ligand_index) + ' original ligand SMILES strings and ' +
                   "{:,}".format(len(alt_ligands))+' alternate ligand SMILES strings! \n' +
                   'Output in: \n'+output_file+'\n')


########################################################################################################
# FIND DRUG-LIKE, METABOLITE-LIKE, and ION-LIKE LIGANDS
########################################################################################################

def check_inputs(tanimoto_file=DATAPATH+'drugbank/drugbank_tanimoto.tsv.gz',
                 check_parsed_file=DATAPATH+'drugbank/drugbank-parsed.tsv',
                 parse_function=download_and_parse_drugbank,
                 smiles_index=6):
  """
  :param tanimoto_file: full path to a tab-delimited file containing the all-against-all pairwise Tanimoto
                        coefficients between SMILES strings
  :param check_parsed_file: full path to a parsed tab-delimited file of all ligand information downloaded from an
                      online source (e.g., DrugBank, HMDB)
  :param parse_function: function to download and parse required raw data about alternate ligands
  :param smiles_index: 0-indexed column of the canonical SMILES string in the parsed_file
  :return: None, but run appropriate functions to guarantee that necessary input files are available
           for the next step (compare_ligands_to_alternate_molecules)
  """

  # parsed information about all ligands in wwPDB's Chemical Components Dictionary
  original_ligand_file = DATAPATH + 'ccd/chemical_components_dictionary-parsed.tsv'
  original_smiles_index = 6  # 0-indexed column of original_ligand_file corresponding to the canonical SMILES string

  # check the Tanimoto coefficient file:
  if not os.path.isfile(tanimoto_file):
    sys.stderr.write('No such file: ' + tanimoto_file + '\n')
    sys.stderr.write('Calculating Tanimoto coefficients...\n')

    # check the parsed raw data file:
    if not os.path.isfile(check_parsed_file):
      sys.stderr.write('Downloading and parsing raw data into ' + check_parsed_file + '\n')
      parse_function()

    # and the original ligand file:
    if not os.path.isfile(original_ligand_file):
      sys.stderr.write('Downloading and parsing raw data into ' + original_ligand_file + '\n')
      download_and_parse_ccd()

    # and calculate the Tanimoto coefficient:
    calculate_tanimoto(parsed_file, smiles_index, parsed_file[:parsed_file.rfind('-parsed.tsv')]+'-smiles.smi',
                       tanimoto_file, original_ligand_file, original_smiles_index)


########################################################################################################

def similar_ligands(tanimoto_files=(DATAPATH+'drugbank/drugbank_tanimoto.tsv.gz',), tanimoto_cutoff=0.9,
                    restriction_group=None, total_processed=0):
  """
  :param tanimoto_files: full paths to tab-delimited files containing ligand IDs, their alternate IDs,
                        and their Tanimoto coefficients
  :param tanimoto_cutoff: float 0->1 specifying the minimum required Tanimoto coefficient to consider
  :param restriction_group: set of alternate ligand IDs to restrict results to
  :param total_processed: int corresponding to the number of lines we have processed so far
  :return: a set of original ligand IDs that passed the Tanimoto cutoff
  """

  ligand_group = set()
  for tanimoto_file in tanimoto_files:
    tanimoto_handle = gzip.open(tanimoto_file) if tanimoto_file.endswith('gz') else open(tanimoto_file)
    try:
      for tanimoto_line in tanimoto_handle:
        if tanimoto_line.startswith('#'):
          continue

        ligand_id, _, _, alt_id, _, _, tanimoto_coefficient = tanimoto_line[:-1].split('\t')[:7]

        if tanimoto_coefficient >= tanimoto_cutoff and (not restriction_group or alt_id in restriction_group):
          ligand_group.add(ligand_id)

        total_processed += 1
        if total_processed % 100000 == 0:
          sys.stderr.write('Processed '+"{:,}".format(total_processed)+' lines\n')
    except IOError:
      sys.stderr.write('Error reading '+tanimoto_file+'\n')
    tanimoto_handle.close()

  return ligand_group, total_processed


########################################################################################################

def edit_ligand_name_string(ligand_name):
  """
  :param ligand_name: string corresponding to a complete ligand name or word from a ligand name descriptor
  :return: a "cleaned" version of the input ligand_name without extraneous characters
  """

  cleaned_ligand_name = ligand_name
  for extraneous_character in ['-', ',', ';', '.', '(', ')']:
    cleaned_ligand_name = cleaned_ligand_name.replace(extraneous_character, '')

  return cleaned_ligand_name


########################################################################################################

def compare_ligands_to_alternate_molecules(metabolite_infiles, drugbank_infiles, tanimoto_cutoff=0.9):
  """
  :return: 3 sets corresponding to mmCIF IDs that can be classified as DRUGLIKE_, METABOLITE_, or ION_
  """

  # DrugBank information up-to-date?
  drug_group, total_processed = similar_ligands(drugbank_infiles, tanimoto_cutoff)

  # HMDB information up-to-date?
  hmdb_raw_data = DATAPATH+'hmdb/human_metabolome_database-parsed.tsv'

  # Extract the subset of ENDOGENOUS human metabolites:
  allowed_ligand_ids = set()
  hmdb_handle = gzip.open(hmdb_raw_data) if hmdb_raw_data.endswith('gz') else open(hmdb_raw_data)
  for hmdb_line in hmdb_handle:
    if hmdb_line.startswith('#'):
      continue

    ligand_id, _, metabolite_category = hmdb_line[:-1].split('\t')[:3]

    if 'endogenous' in metabolite_category:
      allowed_ligand_ids.add(ligand_id)
  hmdb_handle.close()

  metabolite_group = similar_ligands(metabolite_infiles, tanimoto_cutoff, allowed_ligand_ids, total_processed)

  # Finally, extract the ions:
  ion_group = set()

  ligand_raw_data = DATAPATH+'ccd/chemical_components_dictionary-parsed.tsv'
  ligand_handle = gzip.open(ligand_raw_data) if ligand_raw_data.endswith('gz') else open(ligand_raw_data)
  for ligand_line in ligand_handle:
    # lowercase and remove comma separators from the ligand name and the list of synonyms
    if ligand_line.startswith('#') or len(ligand_line[:-1].split('\t')) < 3:
      continue

    ligand_id, ligand_name, ligand_synonyms = [s.lower().replace(',', '') for s in ligand_line[:-1].split('\t')[:3]]
    descriptor = ligand_name.split() + ligand_synonyms.split()  # each word becomes a separate entry

    if 'ion' in [edit_ligand_name_string(ligand_word) for ligand_word in descriptor]:
      ion_group.add(ligand_id)
  ligand_handle.close()

  return drug_group, metabolite_group, ion_group


########################################################################################################

def create_ligand_group_list(metabolite_infiles, drugbank_infiles,
                             lig_grp_outfile=DATAPATH+'ligand_groups.txt', tanimoto_cutoff=0.9):
  """
  :param metabolite_infiles: list of full paths to tab-delimited files containing BioLiP->HMDB Tanimoto coeff.
  :param drugbank_infiles: list of full paths to tab-delimited files containing BioLiP->DrugBank Tanimoto coeff.
  :param lig_grp_outfile: full path to a file where a tab-delimited list of group names will be written
  :param tanimoto_cutoff: cutoff to determine whether mmCIF ligands matched another
  :return: None, but write to the output file a tab-delimited list of group names to mmCIF identifiers
  """

  # get the sets of drugs, metabolites, and ions from the corresponding downloaded files
  drugs, metabolites, ions = compare_ligands_to_alternate_molecules(metabolite_infiles,
                                                                    drugbank_infiles,
                                                                    tanimoto_cutoff)

  out_handle = gzip.open(lig_grp_outfile, 'w') if lig_grp_outfile.endswith('gz') else open(lig_grp_outfile, 'w')
  out_handle.write('\n'.join(['# Groupings of all ligands that may be found in the BioLiP structural database',
                              '# The 4 nucleic acids entries ("NUCACID_") include: ' +
                              'NUCDNA, NUCDNAB, NUCRNA, NUCRNAB',
                              '# The 2 DNA entries ("DNA_") include: NUCDNA, NUCDNAB',
                              '# The 2 RNA entries ("RNA_") include: NUCRNA, NUCRNAB',
                              '# All '+str(len(ions))+' ion entries ("ION_") were extracted from wwPDB\'s Chemical ' +
                              'Component Dictionary (http://www.wwpdb.org/data/ccd),',
                              '# All '+str(len(drugs))+' drug entries ("DRUGLIKE_") were identified by comparing ' +
                              'ligand SMILES strings (Tanimoto coefficient > '+str(tanimoto_cutoff)+') ',
                              '#  to drug SMILES strings obtained from DrugBank ' +
                              '(https://www.drugbank.ca/releases/latest)',
                              '# All '+str(len(metabolites))+' metabolite entries ("METABOLITE_") were identified by ' +
                              'comparing ligand SMILES strings (Tanimoto coefficient > '+str(tanimoto_cutoff)+') ',
                              '#  to *endogenous metabolite* SMILES strings from the Human Metabolome Database',
                              '#  (http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip)'])+'\n')
  out_handle.write('\t'.join(['#group_name', 'original_ligand_mmCIF_identifier'])+'\n')

  for bligand in ['NUCDNA', 'NUCDNAB', 'NUCRNA', 'NUCRNAB']:
    out_handle.write('NUCACID_' + '\t' + bligand + '\n')
  for bligand in ['NUCDNA', 'NUCDNAB']:
    out_handle.write('DNA_' + '\t' + bligand + '\n')
  for bligand in ['NUCRNA', 'NUCRNAB']:
    out_handle.write('RNA_' + '\t' + bligand + '\n')
  for bligand in sorted(list(ions)):
    out_handle.write('ION_' + '\t' + bligand + '\n')
  for bligand in sorted(list(metabolites)):
    out_handle.write('METABOLITE_' + '\t' + bligand + '\n')
  for bligand in sorted(list(drugs)):
    out_handle.write('DRUGLIKE_' + '\t' + bligand + '\n')
  out_handle.close()

  sys.stderr.write('Wrote to '+lig_grp_outfile+'\n')


########################################################################################################

if __name__ == "__main__":

  # Parse the command-line arguments
  parser = argparse.ArgumentParser(description='Compare BioLiP ligands (mmCIFs) to endogenous human metabolites ' +
                                               'and known drugs.')
  parser.add_argument('--database', type=str,
                      help='Either "drugbank" or "hmdb" specifying which ligand set we are going to process',
                      choices=['hmdb', 'drugbank'],
                      default='drugbank')
  parser.add_argument('--force', dest='force', action='store_true', default=False,
                      help='Forcibly overwrite output files that already exist; otherwise skip')
  parser.add_argument('--parse_raw', dest='parse_raw', action='store_true', default=False,
                      help='Download and parse "raw" XML files from the given database')
  parser.add_argument('--tanimoto', dest='tanimoto', action='store_true', default=False,
                      help='Make all-against-all pairwise comparisons using the Tanimoto coefficient')
  parser.add_argument('--create_group_list', dest='create_group_list', action='store_true', default=False,
                      help='Create an output file containing all information from ')
  parser.add_argument('--start', type=int, default=1,
                      help='Starting index (i.e., 1-index) of CCD ligand list to run babel on')
  parser.add_argument('--end', type=int, default=27375,
                      help='Ending index (i.e., 1-index) of CCD ligand list to run babel on')
  parser.add_argument('--tanimoto_cutoff', type=float, default=0.9,
                      help='Tanimoto cutoff between SMILES strings to consider two ligands "the same"')
  args = parser.parse_args()

  # ----------------------------------------------------------------------------------------------------
  # Download and parse raw data if it doesn't already exist
  if args.parse_raw:
    if args.database == 'hmdb' and (args.force or
                                    not os.path.isfile(DATAPATH+'hmdb/human_metabolome_database-parsed.tsv')):
      download_and_parse_hmdb(DATAPATH+'hmdb/', DATAPATH+'hmdb/human_metabolome_database-parsed.tsv')
    elif args.database == 'drugbank' and (args.force or not os.path.isfile(DATAPATH+'drugbank/drugbank-parsed.tsv')):
      download_and_parse_drugbank(DATAPATH+'drugbank/', DATAPATH+'drugbank/drugbank-parsed.tsv')

    # and no matter what, download and parse mmCIF ligands (corresponding to BioLiP):
    if args.force or not os.path.isfile(DATAPATH+'ccd/chemical_components_dictionary-parsed.tsv'):
      download_and_parse_ccd(DATAPATH+'ccd/', DATAPATH+'ccd/chemical_components_dictionary-parsed.tsv')

  # ----------------------------------------------------------------------------------------------------
  # Calculate pairwise Tanimoto coefficients for subsets of CCD ligands
  elif args.tanimoto:

    # (1) Make sure that the parsed files are where they belong
    if not os.path.isfile(DATAPATH+'ccd/chemical_components_dictionary-parsed.tsv'):
      sys.stderr.write('Could not find parsed list of ligands from wwPDB\'s Chemical Components Dictionary in ' +
                       DATAPATH+'ccd/chemical_components_dictionary-parsed.tsv.\n' +
                       'Please run: python group_ligand_types.py --parse_raw\n')
      sys.exit(1)
    if args.database == 'hmdb' and not os.path.isfile(DATAPATH+'hmdb/human_metabolome_database-parsed.tsv'):
      sys.stderr.write('Could not find parsed list of ligands from the Human Metabolome Database in '+DATAPATH +
                       'hmdb/human_metabolome_database-parsed.tsv.\n' +
                       'Please run: python group_ligand_types.py --parse_raw --database hmdb\n')
      sys.exit(1)
    if args.database == 'drugbank' and not os.path.isfile(DATAPATH+'drugbank/drugbank-parsed.tsv'):
      sys.stderr.write('Could not find parsed list of ligands from DrugBank in '+DATAPATH +
                       'drugbank/drugbank-parsed.tsv.\n' +
                       'Please run: python group_ligand_types.py --parse_raw --database drugbank\n')
      sys.exit(1)

    # (2) Extract SMILES strings from parsed file (return in order)
    if args.database == 'hmdb':
      parsed_file = DATAPATH+'hmdb/human_metabolome_database-parsed.tsv'
      smiles_index = 5
    elif args.database == 'drugbank':
      parsed_file = DATAPATH+'drugbank/drugbank-parsed.tsv'
      smiles_index = 3
    else:
      sys.stderr.write('Unknown database option: '+args.database+'...Exiting.\n')
      sys.exit(1)

    smiles_file = parsed_file[:parsed_file.rfind('-parsed.tsv')] + '-smiles.smi'
    ligand_names_file = parsed_file[:parsed_file.rfind('-parsed.tsv')] + '-ligand_names.txt'
    if args.force or not os.path.isfile(smiles_file) or not os.path.isfile(ligand_names_file):
      process_alternate_ligands(parsed_file, smiles_index, smiles_file, ligand_names_file)

    # (2) Run babel to compute Tanimoto coefficients!
    tanimoto_output_file = DATAPATH + args.database + '/' + args.database + '_tanimoto_' + str(args.start) + '-' + \
                           str(args.end) + '.tsv.gz'

    if not os.path.isfile(tanimoto_output_file) or args.force:
      calculate_tanimoto(parsed_file,  # tab-delimited file containing SMILES strings and additional info
                         smiles_file,  # temporary file containing all SMILES strings from the selected database
                         ligand_names_file,  # file containing all ligand info (same order as SMILES file)
                         tanimoto_output_file,  # full path to the output file containing the pairwise Tanimoto coeff.
                         # tab-delimited file containing SMILES strings and additional info for all BioLiP ligands
                         DATAPATH+'ccd/chemical_components_dictionary-parsed.tsv',
                         6,  # 0-index of the column containing the SMILES string in the BioLiP ligands file
                         args.start,  # 1-index of the starting (BioLiP) ligand to compute on
                         args.end)  # 1-index of the ending (BioLiP) ligand to compute on

  # ----------------------------------------------------------------------------------------------------
  # Create a single file containing a list of groups and corresponding BioLiP ligands that fall into those groups
  elif args.create_group_list:
    # (1) Make sure that Tanimoto coefficients have been properly calculated:
    dbinfiles = {}
    for db in ['hmdb', 'drugbank']:
      infiles = [DATAPATH+db+'/'+f for f in os.listdir(DATAPATH+db)
                 if f.startswith(db+'_tanimoto') and f.endswith('.tsv.gz')]
      if len(infiles) < 1:
        sys.stderr.write('Could not find '+DATAPATH+db+'/'+db+'_tanimoto*.tsv.gz\n' +
                         'Please run: python group_ligand_types --tanimoto --database '+db+'\n')
        sys.exit(1)
      dbinfiles[db] = infiles

    # (2) Create the output file!
    create_ligand_group_list(dbinfiles['hmdb'], dbinfiles['drugbank'],
                             DATAPATH + 'ligand_groups-20180925-'+str(args.tanimoto_cutoff)+'.txt',
                             args.tanimoto_cutoff)
