#!/usr/bin/env python

"""
Given a simple Newick formated phylo tree and string markers of two clades, compute the mean
clade to clade distance between members of each clade and the first member of the other
clade respectively
"""

import sys, re, os

# Edit these to reflect where tree_to_dist.py really lives
HOME = os.path.expanduser("~")
PROCESS_DATA=HOME+"/etseq/process_data"

PAIR_PATN = re.compile("\((.*)\)")

OVERLAP_TEST = False  # useful to make sure clade labels don't overlap, e.g hpEurope, hpEuropeSahul
                      # but not such use thereafter
PRINT_ENDS = False

def init() :
  if not os.path.isfile(PROCESS_DATA+"/tree_to_dist.py") :
    print("The application tree_to_dist.py cannot be found", file=sys.stderr)
    sys.exit(1)

  args = sys.argv
  if len(args) != 4 :
    print("Usage: {0:s} <Clade string> <Clade string> <Newick format tree>".format(args[0]), file=sys.stderr)
    print("where the Clade strings are string (incl regular expressions) that differentiate member of the two clades", file=sys.stderr)
    print("both from each other, but also from other clades", file=sys.stderr)
    sys.exit(0)

  if not os.path.isfile(args[3]) :
    print("The Newick formated phylo tree file {0:s} does not exist".format(args[3]), file=sys.stderr)
    sys.exit(1)

  cladestrs = [args[1], args[2]]
  clade_REs = [re.compile(cladestrs[0]), re.compile(cladestrs[1])]
  if OVERLAP_TEST :
    cladestr1_count, cladestr2_count = test_tree(clade_REs, args[3])
    if cladestr1_count * cladestr2_count == 0 :
      if cladestr1_count == 0 :
        print("There are no taxa matching {0:s}".format(cladestrs[0]), file=sys.stderr)
      if cladestr2_count == 0 :
        print("There are no taxa matching {0:s}".format(cladestrs[1]), file=sys.stderr)
      sys.exit(1)
      print("Count of taxa for each clade: {0:s}:{1:d}, {2:s}:{3:d}".format(cladestrs[0], cladestr1_count, cladestrs[1], cladestr2_count), file=sys.stderr)
  return(clade_REs, cladestrs, args[3])

def test_tree(clade_REs, treefile) :
  infile = open(treefile, 'r')
  first_tree = infile.readline()
  infile.close()
  taxa_list = clade_REs[0].findall(first_tree)
  count0 = len(taxa_list)
  for t in taxa_list:
    if clade_REs[1].search(t) :  # overlap, bad
      print("The two clade strings overlap at string {0:s}".format(t), file=sys.stderr)
      sys.exit(1)
  taxa_list = clade_REs[1].findall(first_tree)
  count1 = len(taxa_list)
  for t in taxa_list:
    if clade_REs[0].search(t) :  # overlap, bad
      print("The two clade strings overlap at string {0:s}".format(t), file=sys.stderr)
      sys.exit(1)
  return(count0, count1)

# Get average distance from Clade1 to first instance Clade2, and each Clade2
# to first instance of Clade1.
def get_dists(clade_REs, cladestrs, treefile) :
  dists = {0:0.0, 1:0.0}
  ntaxa =  {0:0, 1:0}
  infile = os.popen("{0:s}/tree_to_dist.py matrix -sort_matrix -matrix_sep ',' {1:s}".format(PROCESS_DATA, treefile))
  for line in infile : # other than blank lines, each line is a node and dists to that node
    line = line.strip()
    if line == "" :
      continue
    fields = line.split(',')
    if clade_REs[0].search(fields[0]) :
      head = 0
      other = 1
    elif clade_REs[1].search(fields[0]) :
      head = 1
      other = 0
    else:
      continue
    # print("head", head, fields[0])
    # Look for the first instance of other's RE and thence distance to head
    for pair in fields[1:] :
      mobj = PAIR_PATN.match(pair)
      node, dist = mobj.group(1).split()
      # print(node, dist)
      if clade_REs[other].search(node) != None :
        # print("got it  {0:s} to {1:s} {2:0.3f}\n".format(cladestrs[head], cladestrs[other], float(dist)))
        dists[head] += float(dist)
        ntaxa[head] += 1
        break
  infile.close()
  return(dists[0]/ntaxa[0], dists[1]/ntaxa[1])
      
  

if __name__ == "__main__" :
  clade_REs, cladestrs, treefile = init() 
  c0t0c1, c1toc0 = get_dists(clade_REs, cladestrs, treefile)
  if PRINT_ENDS :
    print("{0:s}-{1:s}\t{2:0.3f}".format(cladestrs[0], cladestrs[1], min(c0t0c1, c1toc0)))
  else:
    print("{0:0.3f}".format( min(c0t0c1, c1toc0)))



