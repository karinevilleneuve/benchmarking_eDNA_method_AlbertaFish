# General description 

Script and associated data used to benchmark sets of methods frequently employed to analyse eDNA sequencing outputs.

# Required files, scripts and programs

## Programs

- seqkit v2.1 ([website](https://bioinf.shenwei.me/seqkit/))
- RDP tool classifier v2.11 ([Github](https://github.com/rdpstaff/classifier) / [SourceForge](https://sourceforge.net/projects/rdp-classifier/))
- usearch v11.0.667 ([website](http://www.drive5.com/usearch/))
- vsearch v2.15.2 ([Github](https://github.com/torognes/vsearch))
- Barque V1.8.5 ([Github](https://github.com/enormandeau/barque)) (N.B. Barque will only work on GNU Linux or OSX). Requires the following dependencies : 
  - bash 4+
  - python 3.5+
  - python distutils package
  - R 3+ (ubuntu/mint: sudo apt-get install r-base-core)
  - java (ubuntu/mint: sudo apt-get install default-jre)
  - gnu parallel
  - flash (read merger) v1.2.11+
  - vsearch v2.14.2+ (Barque will not work with older versions of vsearch)

## Databases

The following databases can be downloaded with clickiing the hyperlinks

- [COI formatted for Barque](https://www.ibis.ulaval.ca/services/bioinformatique/barque_databases/)
- [Eukaryote CO1 Classifier for RDP V5.1](https://github.com/terrimporter/CO1Classifier/releases/tag/RDP-COI-v5.1.0) 
- [12S fish Classifier v1.0.1](https://github.com/terrimporter/12SfishClassifier/releases/tag/v1.0.1)
- [12S formatted for barque](https://github.com/enormandeau/barque/blob/master/03_databases/12S.fasta.gz)

## Scripts

In order to train the RDP classifier on any database the following scripts are required and can be found [here](https://github.com/karinevilleneuve/fish_eDNA/tree/master/scripts). 

- lineage2taxTrain.py
- addFullLineage.py

# Contact 

For any questions or comments : karine.villeneuve12@gmail.com
