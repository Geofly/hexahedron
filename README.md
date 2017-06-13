[![Hexahedron][hh_logo]][hh_git-repo-url]
___

# Description
Hexahedron is a novel deterministic algorithm, based on HTS data that can reconstruct local and global sequences and determine their relative frequency at a massive scale. The input of the algorithm is the result of the alignment of the reads against a set of reference sequences which are stably remapped to the multiple sequence alignment common coordinate system. The algorithm proceeds in a step-wise manner following a 5’ to 3’ directionality of the references. In every step the variant calling profile along the reference frame is constructed and mutations trigger branching of the variant calling profile and re-allocation of the alignments.

We also offer a novel visualization technique that comprehensively represents the dynamic nature of the results with a simple interactive interface. The graph contains, per position depth of coverage, information for each reconstructed sequence and the first (bifurcation) and last (merging) position that differentiates a sequence and its closest sequence from which it has been derived. Annotations of the provided references, if included in the input set, will also be plotted in the same coordinate system as the reconstructed sequences. Finally, the interface allows sequences to be combined with each other following the paths that connect different bifurcation and merging positions as well as their consensuses and path compositions to be extracted.

![Pipeline][hh_pipeline]


### Requirement

The algorithm has been implemented as part of [HIVE][hive_git-repo-url] platform, a robust infrastructure for next-generation sequence (NGS) data analysis co-developed by Food and Drug Administration and George Washington University
[![HIVE][hive_logo]][hive_git-repo-url] 

Hexahedron itself is open source with a [public repository][hh_git-repo-url] on GitHub.

### Installation

Hexahedron includes the required code base for HIVE and is installed as an application during HIVE installation. 
HIVE installation instructions can be [here][hive_readme].


   [hh_pipeline]: <https://raw.githubusercontent.com/kkaragiannis/hexahedron/master/doc/images/Resized_overview_final_fig.png>
   [hive_readme]: <https://raw.githubusercontent.com/kkaragiannis/hexahedron/master/HIVE_README.md>
   [hive_logo]: <https://raw.githubusercontent.com/FDA/fda-hive/master/doc/images/hive_logo.png>
   [hive_git-repo-url]: <https://github.com/FDA/fda-hive>
   [hh_logo]: <https://raw.githubusercontent.com/kkaragiannis/hexahedron/master/doc/images/hexahedron_logo.png>
   [hh_git-repo-url]: <https://github.com/kkaragiannis/hexahedron>
