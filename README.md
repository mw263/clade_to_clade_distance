# clade_to_clade_distance
                           Clade to Clade Distance

The aim of these scripts is, given a Newick format input tree, to compute the shortest distance
between two clades. The top script, clade_to_clade_dist.py, expects 3 arguments. The first two
arguments are strings (in fact regular expressions) that denote the respective clades. The
final argument is the input tree.

The ksh script run_clade_to_clade_dist provides an example of how the clade_to_clade_dist.py is being
used to find the shortest distance between H pylori clades hpEuropeSahul and hpSahul versus
hpEuropeSahul and hpEurope.

You will need to edit clade_to_clade_dist.py to make sure that it looks in the right place for
the script that it calls, tree_to_dist.py. Please note that clade_to_clade_dist.py is written in
Python 3, but  tree_to_dist.py and a script it imports, newick, are older script written in Python 2.

The newick module is fragile, and if I had time  I would port  tree_to_dist.py to Python 3 and
replace newick with dentropy. newick is hard to come by, so I have included it here.

The code (other than the newick module) is copyright Michael J. Wise, and provided on a CC BY-NC license.
