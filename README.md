# Generate Site-Based Domain Interaction Scores
<div align="center"><img src="interaction-domains-pipeline.png" alt="pipeline figure" title="pipeline to generate site-based domain interaction scores" width="50%" /></div>

In this project, we use information from protein co-complex structures to determine how individual positions 
within domains may be involved in binding different ligands. If you use data or scripts from this repository, 
please cite:

> Kobren, S.N. and Singh, M. (2017) "Structure-informed approach to discovering perturbed interaction interfaces in cancer." *Manuscript in preparation.*

### 1: Downloading preliminary data

* To download the *primary set* of BioLiP data, released March 6, 2013, run the following. **Note that this 
step can take a long time and will require substantial hard drive space.**

```bash
python download_biolip.py --initialize
```


* To update the BioLiP data that you may have already downloaded (as BioLiP releases weekly updates), run:

```bash
python download_biolip.py
```

### 2: Computing distances between atoms

* To calculate pairwise Euclidean distances between receptor residue atoms and ligand atoms, run the following. 
**NOTE: This step also takes a long time; we suggest running in parallel. The --prefix option allows you to specify 
a subset of PDB IDs that begin with a specific 2-character prefix (e.g., 2m).**

```bash
python calculate_distances.py --prefix XX
```

* For some distance calculations, you may want to use the computed overlap area between two partially overlapping 
Gaussian distributions (see the table below for examples). The previously-computed distance files can be updated
to include these *optional* values by running:

```bash
python calculate_distances.py --update_overlap --prefix XX
```

### 3: Finding site-based structural scores

There are multiple ways to assign a continuous ligand-specific binding score to each amino acid residue in a 
protein receptor chain. These scores are summarized in the table below, where "receptor-ligand atom pair" always
refers to a heavy (i.e., non-hydrogen) atom from the *side chain* of the amino acid residue at a particular position
in the protein receptor chain, paired with a heavy atom from the ligand.

| Abbreviation | Score Name | Brief Description | 
| ------------ | -------- | :----------------- |
| *mindist* | minimum distance | minimum distance (in &#8491;) across all receptor-ligand atom pairs for each position in the protein receptor chain
| *meandist* | average distance | average distance across all receptor-ligand atom pairs (per protein receptor chain position). **Note that distances >20&#8491; are stored as 20&#8491;, often biasing this score toward 20&#8491;**
| *fracin4* | fraction within 4&#8491; | fraction of heavy, side-chain atoms that are within 4.0&#8491; of the ligand |
| *maxstd* | maximum overlap area between normals with &sigma;=1.5 | maximum total overlap area (integrated from -&#x221e; to &#x221e;) between each heavy side-chain receptor atom and *all* heavy ligand atoms, where the two Gaussian distributions &phi;(&mu;=0; &sigma;=1.5) and &phi;(&mu;=D; &sigma;=1.5) are centered at each atom in the receptor-ligand atom pair, and where D is the Euclidean distance between the two specified atoms | 
| *meanstd* | average " &sigma;=1.5 | as above, but the **average** (across all heavy side-chain receptor residue atoms) total overlap area  | 
| *sumstd* | total " &sigma;=1.5 | as above, but the overall **total** overlap area between all pairs of heavy side-chain receptor residue atoms and all heavy ligand atoms. **Note that this score is biased toward larger amino acids.** | 
| *maxvdw* | maximum " &sigma;=vdw radii | the maximum total overlap area (as in *maxstd*) but where the standard deviations of each Gaussian distribution in a receptor-ligand atom pair are set to the van der Waals interaction radii of the respective atoms (rather than a uniform 1.5).
| *meanvdw* | average " &sigma;=vdw radii | as above, but the **average** total overlap area |
| *sumvdw* | total " &sigma;=vdw radii | as above, but the overall **total** overlap area |

We suggest using the simplest and most intuitive *mindist* score, which performed as well as the other scores for our purposes and ran in a fraction of the time.

* To calculate a binding score for each protein chain receptor amino acid residue for all appropriate ligand types, run:

```bash
python create_fasta.py --distance mindist --prefix XX
```

### 4: Finding domains in protein receptor chains

You can use whatever method you prefer to find domains in your protein receptor chain sequences. 

In our paper, we downloaded all [Pfam-A (version 31)](http://pfam.xfam.org) domains and ran [HMMER](http://hmmer.org) 
locally to find domain hits. General code to run this step of the pipeline can be found at 
<https://github.com/Singh-Lab/run-hmmer>.

The output of these steps, run on BioLiP (version 2017-06-28) can be found in the following file. **Note:** If you 
choose to run this domain-finding step independently, you must format the results to match the tab-delineated
formatting in the file below (to run subsequent steps of the pipeline).

```
processed_data/domains/BioLiP_2017-06-28-domains-pfam_v31.tsv.gz
```

### 5: Assigning site-based domain binding potential scores

* The "uniqueness" of each domain sequence must be assessed to account for structural redundancies across 
PDB entries. To do this, we generate multiple sequence alignments for each domain in contact with each type of ligand, 
and then assign per-sequence scores as in 
[S Henikoff & JG Henikoff (1994)](http://dx.doi.org/10.1016/0022-2836(94)90032-9).

```bash
python evaluate_uniqueness.py --create_alignments --distance mindist
python evaluate_uniqueness.py --distance mindist
```

* Finally, we use the per-domain sequence uniqueness evaluations generated in the previous step to assign positional 
weights, for each domain, for each ligand type:

```bash
python generate_domain_scores.py --distance mindist
```

---

*NOTE: The following scripts included in this repository run optional steps.*

### Grouping small molecule ligand types

The previous steps automatically group ligands into the *ION_*, *DRUGLIKE_* and *METABOLITE_* 
super groups using molecular information from the [Chemical Component Dictionary](http://www.wwpdb.org/data/ccd), 
[DrugBank](http://www.drugbank.ca) and the [Human Metabolome Database](http://www.hmdb.ca) respectively, resulting in
 the following file (provided in this repository):
 
```
downloaded_data/ligand_groups.txt
```
 
To repeat these steps to generate an up-to-date version of the file, the following script can be run 
(this code likely requires additional configuration). You first need to install the command-line program 
`babel` from <http://openbabel.org/wiki/Main_Page>.

```bash
python group_ligand_types.py
```

### Plotting domain-based scores

To create violin plots depicting the per-position *mindist* domain-based scores (useful to gain intuition about 
how domains are generally involved in binding particular ligands), run:

```bash
python plot_domain_scores.py --domain PF00096_zf-C2H2 --ligand ZN
```