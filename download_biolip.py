#!/usr/bin/python

"""
Download the most recent version of the BioLiP database. If you use this data, please cite:

Jianyi Yang, Ambrish Roy and Yang Zhang (2013). "BioLiP: a semi-manually curated database for biologically relevant 
ligand-protein interactions." Nucleic Acids Research, 41: D1096-D1103.

Contact Shilpa Nadimpalli Kobren (snadimpa@princeton.edu) with questions
"""

import os
import sys
import argparse
from subprocess import call

# path to where this script is currently located (and to where all data should be stored) -- this can
# be updated
DATAPATH = os.path.dirname(os.path.abspath(__file__))+'/'

########################################################################################################

def download_initial_biolip_data():
  """
  :return: None, but download the initial release of BioLiP data (version 2013-03-6) if it has not yet
           been downloaded
  """

  for basedata in [DATAPATH + 'downloaded_data/annotations/BioLiP.tar.bz2',
                   DATAPATH + 'downloaded_data/receptor_2013-03-6.tar.bz2',
                   DATAPATH + 'downloaded_data/ligand_2013-03-6.tar.bz2']:
    if not os.path.isfile(basedata):
      call(['wget', 'http://zhanglab.ccmb.med.umich.edu/BioLiP/download/' + basedata.split('/')[-1],
            '-O', basedata])
    else:
      sys.stderr.write(basedata.split('/')[-1] + ' has already been downloaded. Continuing...\n')

    # Expand the downloaded tar.bz2 files into the appropriate directories
    sys.stderr.write('Expanding ' + basedata + '\n')
    call(['tar', '-jxvf', basedata, '-C', '/'.join(basedata.split('/')[:-1]) + '/'])


########################################################################################################

def update_biolip_data():
  """
  Download all available weekly updates to the original BioLiP data release
  :return: ordered list strings corresponding to all available weekly updates
  """

  weekly_updates_filename = DATAPATH + 'downloaded_data/weekly_updates.html'
  call(['wget', 'http://zhanglab.ccmb.med.umich.edu/BioLiP/weekly.html', '-O', weekly_updates_filename])

  weekly_updates = []
  weekly_updates_handle = open(weekly_updates_filename)
  for html_line in weekly_updates_handle:
    # check if the date is contained on this line
    if html_line.startswith('<tr><td>') and html_line.strip().endswith('</td>') and len(html_line.strip()) in [22, 23]:
      weekly_updates.append(html_line.strip().replace('<tr><td>', '').replace('</td>', ''))
  weekly_updates_handle.close()
  call(['rm', weekly_updates_filename])

  # Now, download and expand each weekly update that we haven't already downloaded
  for date in weekly_updates:
    if not os.path.isfile(DATAPATH + 'downloaded_data/annotations/BioLiP_' + date + '.txt'):
      call(['wget', 'http://zhanglab.ccmb.med.umich.edu/BioLiP/weekly/BioLiP_' + date + '.txt',
            '-O', DATAPATH + 'downloaded_data/annotations/BioLiP_' + date + '.txt'])

      for tarball_type in ['ligand', 'receptor']:
        call(['wget', 'http://zhanglab.ccmb.med.umich.edu/BioLiP/weekly/' + tarball_type + '_' + date + '.tar.bz2',
              '-O', DATAPATH + 'downloaded_data/' + tarball_type + '_' + date + '.tar.bz2'])
        call(['tar', '-jxvf', DATAPATH + 'downloaded_data/' + tarball_type + '_' + date + '.tar.bz2',
              '-C', DATAPATH + 'downloaded_data/'])
        call(['rm', DATAPATH + 'downloaded_data/' + tarball_type + '_' + date + '.tar.bz2'])

  sys.stderr.write('All BioLiP data through ' + weekly_updates[0] + ' has been downloaded!\n')
  return weekly_updates + ['2013-03-6']  # include the base annotation date


########################################################################################################

if __name__ == "__main__":

  # Parse the command-line arguments
  parser = argparse.ArgumentParser(description='Download and update a local copy of the BioLiP protein--ligand ' +
                                               'database.')

  parser.add_argument('--initialize', dest='initialize', action='store_true',
                      help='Download the starting set of structures and annotations (March 6, 2013)?')
  parser.set_defaults(initialize=False)

  args = parser.parse_args()

  # ----------------------------------------------------------------------------------------------------
  # Download the initial set of structures and annotations, *if need be*
  if args.initialize:
    sys.stderr.write('Downloading original BioLiP data...\n')

    # Create the subdirectory structure (as expected by scripts in this repository)
    for subdirectory in [DATAPATH+'downloaded_data',
                         DATAPATH+'downloaded_data/annotations',
                         DATAPATH+'downloaded_data/ligand',
                         DATAPATH+'downloaded_data/receptor']:
      if not os.path.isdir(subdirectory):
        call(['mkdir', subdirectory])

    # Check to see if starting tar.bz2 files already exist; download and expand otherwise
    download_initial_biolip_data()

  # ----------------------------------------------------------------------------------------------------
  # Otherwise, download the list of updates to BioLiP and download those updates
  else:
    sys.stderr.write('Downloading BioLiP updates...\n')

    # Make sure that the initialization step has already been run
    for subdirectory in [DATAPATH+'downloaded_data',
                         DATAPATH+'downloaded_data/annotations',
                         DATAPATH+'downloaded_data/ligand',
                         DATAPATH+'downloaded_data/receptor']:
      if not os.path.isdir(subdirectory):
        sys.stderr.write('No such directory: '+subdirectory+'\n')
        sys.stderr.write('Please run python '+sys.argv[0]+' --initialize\n')
        sys.exit(1)

    # Get the list of all weekly updates (this does not include the base data) and download corresponding files
    release_dates = update_biolip_data()

    # Create a current concatenated list of *all* annotations
    for subdirectory in [DATAPATH+'processed_data',
                         DATAPATH+'processed_data/annotations']:
      if not os.path.isdir(subdirectory):
        call(['mkdir', subdirectory])

    concatenated_annotation_file = DATAPATH+'processed_data/annotations/current_annotations.txt'
    call(['cat'] + [DATAPATH+'downloaded_data/annotations/BioLiP_'+update+'.txt' for update in release_dates] +
         ['>', concatenated_annotation_file])
