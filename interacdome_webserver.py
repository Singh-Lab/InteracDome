#!/usr/bin/python

"""
Generate final InteracDome output files, filtering based on cross-validation step (if required)

Contact Shilpa Nadimpalli Kobren (snadimpa@princeton.edu) with questions
"""

import os
import sys
import gzip
import math
import argparse
from subprocess import call

########################################################################################################
# CONSTANTS
########################################################################################################

# path to where this script is currently located (and to where all data should be stored) -- this can
# be updated
DATAPATH = os.path.dirname(os.path.abspath(__file__)) + '/'

# full path to all per-domain binding weights (created using generate_binding_scores.py)
SCORE_PATH = DATAPATH + 'processed_data/domains/binding_scores/'

# full path to where output files containing cross-validation evaluations will be written
CV_PATH = DATAPATH + 'processed_data/domains/cross_validation/'


########################################################################################################
# DOMAINS THAT BIND LIGANDS CONSISTENTLY
########################################################################################################

def consistent_domain_ligand_pairs(required_overall_precision=0.5,
                                   required_grouped_precision=0.,
                                   required_instance_count=0,
                                   required_unique_instance_count=3,
                                   required_group_count=0,
                                   required_structure_count=3,
                                   distance='mindist'):
  """
  :param required_overall_precision: minimum required (ungrouped) cross-validated precision that a domain-ligand pair
                                     achieved (at some point) in order for the domain to be used to infer ligand-
                                     binding positions
  :param required_grouped_precision: minimum required (grouped by >90% sequence identiy) cross-validated precision that
                                     a domain-ligand pair achieved (at some point) in order for the domain to be used
                                     to infer ligand-binding positions
  :param required_instance_count: minimum required number of domain instances across BioLiP structures
  :param required_group_count: minimum required number of domain instances with <90% sequence identity to each other
                               across BioLiP structures
  :param required_unique_instance_count: minimum required num of NONIDENTICAL domain instances across BioLiP structures
  :param required_structure_count: minimum required number of unique PDB structures containing 1+ domain instances
  :param distance: default 'mindist' (how to score things)
  :return: set of (domain, ligand) type pairs that meet the required restrictions
  """

  # All the cross-validated results files (including real and random curves):
  accuracy_path = CV_PATH + distance + '/precision_threshold/'
  accuracy_files = sorted([accuracy_path + a for a in os.listdir(accuracy_path)
                           if a.endswith('-pt-' + distance + '.txt.gz')])

  passing_domains = set()

  # Keep track of xs (binding propensities) and ys (precisions achieved) for each ligand type
  for cv_file in accuracy_files:
    cv_handle = gzip.open(cv_file) if cv_file.endswith('gz') else open(cv_file)
    for cv_line in cv_handle:
      if cv_line.startswith('#') or len(cv_line[:-1].split('\t')) < 16:
        continue

      (domain_name, ligand_type, num_matchstates, frac_binding, instance_cnt,
       grouped_propensities, grouped_precisions, _, _, _, _,
       ungrouped_propensities, ungrouped_precisions, _, _, _, _) = cv_line[:-1].split('\t')[:17]

      if len(instance_cnt.split('|')) < 4:  # i.e., improperly formatted cross-validation file
        continue
      num_instances, num_unique_instances, num_groups, num_structures = map(int, instance_cnt.split('|')[:4])

      # do we have the minimum ungrouped precision?
      if required_overall_precision > 0:
        try:
          if max(map(float, ungrouped_precisions.split(','))) < required_overall_precision:
            continue
        except ValueError:
          continue

      # do we have the minimum ungrouped precision?
      if required_grouped_precision > 0:
        try:
          if max(map(float, grouped_precisions.split(','))) < required_grouped_precision:
            continue
        except ValueError:
          continue

      # do we have the appropriate number of instances?
      if num_instances < required_instance_count or \
              num_unique_instances < required_unique_instance_count or \
              num_groups < required_group_count or \
              num_structures < required_structure_count:
        continue

      passing_domains.add((domain_name, ligand_type))
    cv_handle.close()

  sys.stderr.write(str(len(passing_domains)) + ' structurally consistent domain-ligand interactions from ' +
                   str(len(set([a[0] for a in passing_domains]))) + ' domains\n')
  return passing_domains


########################################################################################################

def save_consistent_domain_ligand_pairs(required_overall_precision=0.5,
                                        required_grouped_precision=0.,
                                        required_instance_count=0,
                                        required_unique_instance_count=3,
                                        required_group_count=0,
                                        required_structure_count=3,
                                        minimum_passing_threshold=0.5,
                                        distance='mindist'):
  """
  :param required_overall_precision: minimum required (ungrouped) cross-validated precision that a domain-ligand pair
                                     achieved (at some point) in order for the domain to be used to infer ligand-
                                     binding positions
  :param required_grouped_precision: minimum required (grouped by >90% sequence identiy) cross-validated precision that
                                     a domain-ligand pair achieved (at some point) in order for the domain to be used
                                     to infer ligand-binding positions
  :param required_instance_count: minimum required number of domain instances across BioLiP structures
  :param required_group_count: minimum required number of domain instances with <90% sequence identity to each other
                               across BioLiP structures
  :param required_unique_instance_count: minimum required num of NONIDENTICAL domain instances across BioLiP structures
  :param required_structure_count: minimum required number of unique PDB structures containing 1+ domain instances
  :param minimum_passing_threshold: minimum precision received **FOR A PARTICULAR BINDING PROPENSITY** (note this is
                                    different from the required overall precision for a domain-ligand pair)
  :param distance: default 'mindist' (how to score things)
  :return: None, but print the names of the resulting output files as they are written to
  """

  outfile_suffix = []
  if required_overall_precision > 0:
    outfile_suffix.append('overall-precision-' + str(required_overall_precision))
  if required_grouped_precision > 0:
    outfile_suffix.append('grouped-precision-' + str(required_grouped_precision))
  if required_instance_count > 0:
    outfile_suffix.append('instances-' + str(int(required_instance_count)))
  if required_unique_instance_count > 0:
    outfile_suffix.append('unique-instances-' + str(int(required_unique_instance_count)))
  if required_group_count > 0:
    outfile_suffix.append('grouped-instances-' + str(int(required_group_count)))
  if required_structure_count > 0:
    outfile_suffix.append('structures-' + str(int(required_structure_count)))

  # full path to a tab-delimited file containing a list of domain-ligand pairs that are "conservatively assessed"
  dom_outfile = DATAPATH + 'processed_data/domains/binding_scores/' + distance + '_passing-domains' + \
                ('_' if len(outfile_suffix) > 0 else '') + '_'.join(outfile_suffix) + '.txt'

  if minimum_passing_threshold > 0:
    outfile_suffix.append('sitebased-precision-' + str(minimum_passing_threshold))

  wts_outfile = DATAPATH + 'processed_data/domains/binding_scores/' + distance + '_passing-binding-propensities' + \
                ('_' if len(outfile_suffix) > 0 else '') + '_'.join(outfile_suffix) + '.txt'

  # first, get the set of DOMAINS that pass:
  passing_domains = consistent_domain_ligand_pairs(required_overall_precision,
                                                   required_grouped_precision,
                                                   required_instance_count,
                                                   required_unique_instance_count,
                                                   required_group_count,
                                                   required_structure_count,
                                                   distance)

  # write out the list of DOMAINS that pass to the proper file:
  out_handle = gzip.open(dom_outfile, 'w') if dom_outfile.endswith('gz') else open(dom_outfile, 'w')
  out_handle.write('\n'.join(['# All ' + "{:,}".format(len(passing_domains)) + ' domain-ligand type pairs from ' +
                              "{:,}".format(len(set([a[0] for a in passing_domains]))) +
                              ' domains with binding frequencies found in ',
                              '# ' + SCORE_PATH + distance + '/',
                              '# with a (ungrouped) precision >= ' + str(required_overall_precision),
                              '# a (grouped) precision >= ' + str(required_grouped_precision),
                              '# and ' + str(required_instance_count) + '+ instances, ' +
                              str(required_unique_instance_count) + '+ nonidentical sequences, and ' +
                              str(required_group_count) + '+ groups with <90% sequence identity from ' +
                              str(required_structure_count) + '+ distinct structures',
                              '# of the domain in complex with the corresponding ligand in BioLiP']) + '\n')
  out_handle.write('\t'.join(['#domain_name', 'ligand_type']) + '\n')
  for domname, ligtype in sorted(passing_domains):
    out_handle.write(domname + '\t' + ligtype + '\n')
  out_handle.close()

  sys.stderr.write('All passing domains found in ' + dom_outfile + '\n')

  # now, get the corresponding binding propensities thresholds for each passing domain-ligand pair:
  passing_domain_ligand_pairs = {}
  accuracy_files = sorted([CV_PATH+distance+'/precision_threshold/'+a
                           for a in os.listdir(CV_PATH+distance+'/precision_threshold/')
                           if a.endswith('-pt-' + distance + '.txt.gz')])

  for cv_file in accuracy_files:
    cv_handle = gzip.open(cv_file) if cv_file.endswith('gz') else open(cv_file)
    for cv_line in cv_handle:
      if cv_line.startswith('#') or len(cv_line[:-1].split('\t')) < 13:
        continue

      (domain_name, ligand_type, _, _, instance_cnt, grouped_propensities, grouped_precisions, _, _, _, _,
       ungrouped_propensities, ungrouped_precisions) = cv_line[:-1].split('\t')[:13]

      if (domain_name, ligand_type) not in passing_domains:
        continue

      try:
        curve = zip(map(float, ungrouped_propensities.split(',')), map(float, ungrouped_precisions.split(',')))
      except ValueError:
        continue

      # lowest binding propensity for this particular domain-ligand pair that resulted in the minimum required precision
      truncated_curve = [bp for bp, p in curve if p >= minimum_passing_threshold]
      if len(truncated_curve) > 0:
        passing_domain_ligand_pairs[(domain_name, ligand_type)] = truncated_curve[0]
    cv_handle.close()

  # write out the list of BINDING PROPENSITIES that pass:
  out_handle = gzip.open(wts_outfile, 'w') if wts_outfile.endswith('gz') else open(wts_outfile, 'w')
  out_handle.write('\n'.join(['# All binding frequencies that resulted in a (ungrouped) cross-validated precision >= '+
                              str(minimum_passing_threshold) + ' are retained from ',
                              '# ' + "{:,}".format(len(passing_domains)) + ' domain-ligand type pairs from ' +
                              "{:,}".format(len(set([a[0] for a in passing_domains]))) + ' domains',
                              '# that achieved an overall (ungrouped) precision >= ' + str(required_overall_precision),
                              '# a (grouped) precision >= ' + str(required_grouped_precision),
                              '# and there were ' + str(required_instance_count) + '+ instances, ' +
                              str(required_unique_instance_count) + '+ nonidentical sequences, and ' +
                              str(required_group_count) + '+ groups with <90% sequence identity from ' +
                              str(required_structure_count) + '+ distinct structures',
                              '# of the domain in complex with the corresponding ligand in BioLiP',
                              '# Full list of passing domains found in ' + dom_outfile,
                              '# All (unfiltered) binding frequencies found in ' + SCORE_PATH + distance + '/',
                              '\t'.join(['#domain_name', 'ligand_type', 'match_state', 'binding_frequency'])]) + '\n')

  binding_files = sorted([SCORE_PATH + distance + '/' + a for a in os.listdir(SCORE_PATH + distance)
                          if a.endswith('_binding-scores_' + distance + '.txt.gz')])

  for bind_file in binding_files:
    domain_name = bind_file.split('/')[-1].replace('_binding-scores_' + distance + '.txt.gz', '')

    binding_file_handle = gzip.open(bind_file) if bind_file.endswith('gz') else open(bind_file)
    for bind_info in binding_file_handle:
      if bind_info.startswith('#'):
        continue

      ligand_type, match_state, binding_propensity = bind_info[:-1].split('\t')[:3]

      if (domain_name, ligand_type) in passing_domain_ligand_pairs and \
              float(binding_propensity) >= passing_domain_ligand_pairs[(domain_name, ligand_type)]:
        out_handle.write('\t'.join([domain_name, ligand_type, match_state, binding_propensity]) + '\n')
    binding_file_handle.close()

  out_handle.close()

  sys.stderr.write('Wrote all binding information to ' + wts_outfile + '\n')


########################################################################################################
# INTERACDOME WEBSERVER FILES
########################################################################################################

def datasource_website_input(outfile, pfam_path, distance='mindist', for_webserver_display=True):
  """
  :param outfile: full path to a tab-delimited file where all results should be written out to
  :param pfam_path: full path to a directory containing Pfam-formatted HMMs in
  :param distance: how residue-ligand distances were calculated
  :param for_webserver_display: whether to truncate to only DNA, RNA, peptide, ion, metabolite, small molecule
                                results (as displayed on website) or not
  :return: none, but print results to specified outfile along with success message upon completion
  """

  # NOTE: these domains interacted only via their *insertion* states (never via a match state), listed here
  #   just for record-keeping and to remind you of a reason for potential discrepancies in counts for website
  skipped_insertion_states = set()

  outhandle = open(outfile, 'w')
  outhandle.write('\t'.join(['pfam_id', 'domain_length', 'ligand_type', 'num_nonidentical_instances', 'num_structures',
                             'binding_frequencies',
                             'max_achieved_precision',
                             'frequency_at_precision_0.1',
                             'frequency_at_precision_0.25',
                             'frequency_at_precision_0.5',
                             'frequency_at_precision_0.75',
                             'pdb_ids']) + '\n')

  binding_propensity_files = sorted([SCORE_PATH + distance + '/' + a for a in os.listdir(SCORE_PATH + distance)
                                     if a.endswith('_binding-scores_' + distance + '.txt.gz')])
  for bp_file in binding_propensity_files:

    domain_name = bp_file.split('/')[-1].replace('_binding-scores_' + distance + '.txt.gz', '')
    scores = {}  # ligand_type -> {(x,y), (x,y)...}
    structures = {}  # ligand_type -> {pdb_id, pdb_id...}
    unique_instances = {}  # ligand_type -> number of unique instances
    unique_structures = {}  # ligand_type -> number of unique structures
    best_precisions = {}  # ligand_type -> (max, 0.1, 0.25, 0.5, 0.75) -> maximum achieved cross-validated precision

    # --------------------------------------------------------------------------------------------------
    # get the length of the domain (number of MATCH STATES)
    try:
      with open(pfam_path + domain_name + '.hmm') as hmm:
        for hmm_line in hmm:
          if hmm_line.startswith('LENG'):
            domain_length = hmm_line.strip().split()[-1]
            break
    except:
      sys.stderr.write('Could not find length for '+domain_name+'!\n')
      continue

    # --------------------------------------------------------------------------------------------------
    # get the corresponding PDB structure IDs and binding propensities:
    bp_handle = gzip.open(bp_file) if bp_file.endswith('gz') else open(bp_file)
    for bp_line in bp_handle:
      if bp_line.startswith('#'):
        continue

      ligand_type, match_state, binding_propensity, structure_scores = bp_line[:-1].split('\t')[:4]

      if for_webserver_display:
        if ligand_type not in ['DNA_', 'DNABASE_', 'DNABACKBONE_', 'RNA_', 'RNABASE_', 'RNABACKBONE_',
                               'PEPTIDE_', 'ION_', 'METABOLITE_', 'SM_']:
          continue
        ligand_type = ligand_type.replace('_', '').lower()

      # some domains only had modeled interactions in their insertion states....
      if len(skipped_insertion_states) < 1:
        try:
          match_state = int(match_state)
        except ValueError:
          skipped_insertion_states.add(domain_name)
          continue

      if ligand_type not in scores:
        scores[ligand_type] = {}
      scores[ligand_type][match_state] = binding_propensity

      if ligand_type not in structures:
        structures[ligand_type] = set()
      for struct_id in structure_scores.split(','):
        structures[ligand_type].add(struct_id[:4])
    bp_handle.close()

    # --------------------------------------------------------------------------------------------------
    # finally, get the accuracies for this domain's interactions
    cv_file = CV_PATH + distance + '/precision_threshold/' + domain_name + '-pt-' + distance + '.txt.gz'
    cv_handle = gzip.open(cv_file) if cv_file.endswith('gz') else open(cv_file)
    for cv_line in cv_handle:
      if cv_line.startswith('#') or len(cv_line[:-1].split('\t')) < 16:
        continue

      (domain_name, ligand_type, num_matchstates, frac_binding, instance_cnt,
       grouped_propensities, grouped_precisions, _, _, _, _,
       ungrouped_propensities, ungrouped_precisions, _, _, _, _) = cv_line[:-1].split('\t')[:17]

      if for_webserver_display:
        if ligand_type not in ['DNA_', 'DNABASE_', 'DNABACKBONE_', 'RNA_', 'RNABASE_', 'RNABACKBONE_',
                               'PEPTIDE_', 'ION_', 'METABOLITE_', 'SM_']:
          continue
        ligand_type = ligand_type.replace('_', '').lower()

      if len(instance_cnt.split('|')) < 4:  # i.e., improperly formatted cross-validation file
        continue
      _, unique_instances[ligand_type], _, unique_structures[ligand_type] = instance_cnt.split('|')[:4]

      # store the maximum achieved precision and the per domain-ligand pair cutoffs for different precisions
      best_precisions[ligand_type] = {'max': 0., '0.1': '--', '0.25': '--', '0.5': '--', '0.75': '--'}
      try:
        propensity_to_precision = dict(zip(map(float, ungrouped_propensities.split(',')),
                                           map(float, ungrouped_precisions.split(','))))
        best_precisions[ligand_type]['max'] = max(propensity_to_precision.values())
        for domain_threshold in [0.1, 0.25, 0.5, 0.75]:
          passing_propensities = [propensity for propensity, precision in propensity_to_precision.items()
                                  if precision >= domain_threshold and propensity > 0]
          if len(passing_propensities) > 0:
            best_precisions[ligand_type][str(domain_threshold)] = min(passing_propensities)
      except ValueError:
        pass

    # --------------------------------------------------------------------------------------------------
    # write results to specified output file
    for ligand_type, score_set in scores.items():
      binding_frequencies = [str(score_set.get(str(mstate), 0)) for mstate in xrange(1, int(domain_length) + 1)]
      if set(binding_frequencies) == set(['0']):
        continue
      outhandle.write('\t'.join([domain_name,  # pfam_id
                                 domain_length,  # domain_length
                                 ligand_type,  # ligand_type
                                 unique_instances.get(ligand_type, '--'),  # num_nonidentical_instances
                                 unique_structures.get(ligand_type, '--'),  # num_structures
                                 ','.join(binding_frequencies),  # binding_propensities
                                 str(best_precisions.get(ligand_type, {}).get('max',
                                                                              '--')),  # max_achieved_precision
                                 str(best_precisions.get(ligand_type, {}).get('0.1',
                                                                              '--')),  # propensity_at_precision_0.1
                                 str(best_precisions.get(ligand_type, {}).get('0.25',
                                                                              '--')),  # propensity_at_precision_0.25
                                 str(best_precisions.get(ligand_type, {}).get('0.5',
                                                                              '--')),  # propensity_at_precision_0.5
                                 str(best_precisions.get(ligand_type, {}).get('0.75',
                                                                              '--')),  # propensity_at_precision_0.75
                                 ','.join(sorted(structures[ligand_type]))  # pdb_ids
                                 ]) + '\n')
  outhandle.close()
  sys.stderr.write('Wrote to ' + outfile + '\n')

  print skipped_insertion_states


########################################################################################################

def datasource_website_hmms(pfam_hmm_infile, gglogo_hmm_outfile):
  """
  :param pfam_hmm_infile: full path to a Pfam-formatted HMM file
  :param gglogo_hmm_outfile: corresponding .pfm file
  :return: None, but reformat input Pfam-formatted HMM file and output .pfm formatted HMM file
  """

  if not os.path.isfile(pfam_hmm_infile):
    sys.stderr.write('Could not find ' + pfam_hmm_infile + '\n')
    return

  match_state_starts = False
  match_states = {}  # match_state -> amino_acid -> frequency
  domain_length = 0

  hmm_handle = gzip.open(pfam_hmm_infile) if pfam_hmm_infile.endswith('gz') else open(pfam_hmm_infile)
  for hmm_line in hmm_handle:
    if hmm_line.startswith('LENG'):
      domain_length = int(hmm_line.strip().split()[-1])

    elif hmm_line.startswith('HMM '):
      aas = hmm_line.strip().split()[1:]
      match_state_starts = True
      for _ in xrange(4):
        hmm_handle.next()  # skip the rest

    elif match_state_starts:
      mstate = hmm_line.strip().split()[0]
      try:
        if int(mstate) not in xrange(1, domain_length + 1):
          break
      except ValueError:
        break
      match_states[mstate] = {aas[i]: math.exp(-1 * float(hmm_line.strip().split()[i + 1])) for i in xrange(len(aas))}
      for _ in xrange(2):
        hmm_handle.next()
  hmm_handle.close()

  out_handle = open(gglogo_hmm_outfile, 'w')
  for aa in aas:
    out_handle.write(
      '\t'.join([aa] + [str(match_states[m][aa]) for m in map(str, xrange(1, domain_length + 1))]) + '\n')
  out_handle.close()
  # sys.stderr.write('Converted '+pfam_hmm_infile+' to '+gglogo_hmm_outfile+'\n')


########################################################################################################
# MAIN
########################################################################################################

if __name__ == "__main__":

  # --------------------------------------------------------------------------------------------------
  # parse the command-line arguments
  parser = argparse.ArgumentParser(description='Create lists of filtered binding frequencies that can be ' +
                                               'used to infer ligand-binding residues across new sequences.')

  parser.add_argument('--precision', type=float,
                      default=0.5,
                      help='minimum cross-validated (ungrouped) precision required')
  parser.add_argument('--grouped_precision', type=float,
                      default=0.,
                      help='minimum cross-validated (grouped by >90% sequence identity) precision required')
  parser.add_argument('--instances', type=int,
                      default=0,
                      help='minimum domain instances required')
  parser.add_argument('--unique_instances', type=int,
                      default=3,
                      help='minimum nonidentical domain instances required')
  parser.add_argument('--groups', type=int,
                      default=0,
                      help='minimum domain instances with <90% sequence identity required')
  parser.add_argument('--structures', type=int,
                      default=3,
                      help='minimum distinct PDB structures with 1+ domain instances required')
  parser.add_argument('--threshold_precision', type=float,
                      default=0.5,
                      help='minimum achieved (ungrouped) cross-validated precision for a binding frequency to be ' +
                           'used to infer ligand-binding positions')
  parser.add_argument('--distance', type=str,
                      help='How to record the distance between receptor and ligand?',
                      default='mindist',
                      choices={'fracin4', 'mindist', 'meandist', 'maxstd', 'meanstd', 'sumstd',
                               'maxvdw', 'meanvdw', 'sumstd'})
  parser.add_argument('--pfam_path', type=str,
                      help='full path to a directory containing Pfam-formatted HMMs',
                      default=DATAPATH + 'pfam/hmms-v31/')
  parser.add_argument('--webserver', dest='create_webserver_files', action='store_true', default=False,
                      help='Calculate distance-to-ligand consistencies across 50-50 splits of domain-ligand instances')

  args = parser.parse_args()

  # --------------------------------------------------------------------------------------------------
  # check all appropriate directories / subdirectories
  if not os.path.isdir(SCORE_PATH + args.distance + '/'):
    sys.stderr.write('No such directory: ' + SCORE_PATH + args.distance + '/\n' +
                     'Please run: python generate_domain_scores.py\n')
    sys.exit(1)

  if not os.path.isdir(CV_PATH + args.distance + '/precision_threshold/'):
    sys.stderr.write('No such directory: ' + CV_PATH + args.distance + '/precision_threshold/\n' +
                     'Please run: python cross_validate_scores.py\n')
    sys.exit(1)

  # --------------------------------------------------------------------------------------------------
  # write out results files as specified
  if not args.create_webserver_files:

    sys.stderr.write('\n'.join(['Finding ('+args.distance+') domain-ligand pairs with...',
                                '   minimum ungrouped precision >= '+str(args.precision),
                                '   minimum grouped precision >= '+str(args.grouped_precision),
                                '   '+str(args.instances)+'+ instances in BioLiP structures',
                                '   '+str(args.unique_instances)+'+ nonidentical instances in BioLiP structures',
                                '   '+str(args.groups)+'+ instances with <90% sequence identity in BioLiP structures',
                                '   '+str(args.structures)+'+ distinct PDB structures',
                                '   binding frequencies that resulted in an ungrouped precision >= '+
                                str(args.threshold_precision)])+'\n')

    save_consistent_domain_ligand_pairs(args.precision,
                                        args.grouped_precision,
                                        args.instances,
                                        args.unique_instances,
                                        args.groups,
                                        args.structures,
                                        args.threshold_precision,
                                        args.distance)

  # --------------------------------------------------------------------------------------------------
  # create webserver files if directed to do so:
  elif args.create_webserver_files:

    sys.stderr.write('Creating input files required for InteracDome webserver...\n')

    # reformat Pfam's HMM files in order to be used by gglogo:
    if not os.path.isdir(args.pfam_path):
      sys.stderr.write('No such directory: ' + args.pfam_path + '\n' +
                       'Please specify the full path to all Pfam HMMs using the pfam_path flag as:\n' +
                       'python interacdome_webserver.py --pfam_path <full_path_to_directory>\n')
      sys.exit(1)

    for subdir in ['interacdome-webserver', 'interacdome-webserver/pfms']:
      if not os.path.isdir(DATAPATH + subdir):
        call(['mkdir', DATAPATH + subdir])

    # create the two required input files (to display plots and to create downloads:
    datasource_website_input(DATAPATH + 'interacdome-webserver/interacdome_allresults.tsv',
                             args.pfam_path + ('/' if not args.pfam_path.endswith('/') else ''),
                             args.distance,
                             True)  # restrict output to DNA/RNA/peptide/ion/metabolite/small molecule
    datasource_website_input(DATAPATH + 'interacdome-webserver/interacdome_fordownload.tsv',
                             args.pfam_path + ('/' if not args.pfam_path.endswith('/') else ''),
                             args.distance,
                             False)  # don't restrict any output

    for hmm in sorted(os.listdir(SCORE_PATH + args.distance + '/')):
      if hmm.endswith('_binding-scores_' + args.distance + '.txt.gz'):
        pfam_id = hmm.replace('_binding-scores_' + args.distance + '.txt.gz', '')
        if not os.path.isfile(DATAPATH + 'interacdome-webserver/pfms/' + pfam_id + '.pfm'):
          datasource_website_hmms(args.pfam_path + ('/' if not args.pfam_path.endswith('/') else '') + pfam_id + '.hmm',
                                  DATAPATH + 'interacdome-webserver/pfms/' + pfam_id + '.pfm')
