#!/usr/bin/env python2

"""
This app takes the tree and either returns a matrix of evolutionary distances from
each taxon to very other, or it returns for each node the average distance to all
other nodes, or it returns a single value which is the average of all the averages.

[ lh=-7360.433800 ](3_ACH1_LONAC:1.06464,(3_TRYP_STRGR:0.47569,3_PRTA_STRGR:1.34968)
99:0.59535,3_PLMN_HUMAN:0.67404,3_THRB_BOVIN:0.84114,3_EL1_BOVIN:0.74261,
3_CTRA_BOVIN:0.68527,33__TRY1_BOVIN:0.52911);

or 
[ lh=-893.893018 ](H_EU263984_H5N1_1:0.00001,(((((H_CY016237_H1N1_1:0.00386,H_CY016229_H1N1_1:0.00386,
H_CY013814_H1N1_1:0.00001)100:0.01171,H_CY020574_H1N1_1:0.00001)
99:0.00386,H_CY006764_H3N2_1:0.00001)99:0.00386,(((H_DQ508836_H3N2_1:0.00387,
H_CY003753_H3N2_1:0.00001)89:0.00387,(H_CY006852_H3N2_1:0.00386,
H_CY006204_H3N2_1:0.00001)78:0.00001)95:0.00387,H_CY007972_H3N2_1:0.00001)
94:0.00387)52:0.00001,H_CY008157_H3N2_1:0.00001)81:0.05100,H_EU434687_H5N1_1:0.00386);

Input is a text file containing the tree. It is assumed that branch lengths follow
taxon IDs following a colon, or they follow the closing bracket of a subtree. In the
latter case the branch length can be preceded by a bootstrap percentage, separated by
a hash
"""

import sys, re, string, types, math
import newick.tree

ID_PATN = re.compile("([^,:() ]*[a-zA-Z_][^,:() ]*)\s*:")

DEBUG = False

# For the file summaries, average the per sequence median distances rather than means
AVERAGE_MEDIANS = True

LARGE_NUMBER = 100000000

# multiply MAD by this to give unbiased estimator of SD
MAD_MULTIPLIER = 1.48260221850560

params_dict = {"action": "per_seq", "nominalN": None, "MISL":None, "median":False, "sort_matrix":False,
		"sort_matrix_descending":False, "include_branches":False, "matrix_sep":" "}

"""
node_map is a global var
node_map: node -> (parent node, brlen)
where nodes are either identifiers or (identifiers), where the latter denote internal nodes 
"""
node_map = {}

def init():
  global params_dict # so it can be imported by other modules

  args = sys.argv
  if len(args) == 1:
    sys.stderr.write("Usage: %s [action] [OPTIONS] <newick format text file>\n" % args[0])
    sys.stderr.write("The possible actions are:\n")
    sys.stderr.write("\tper_seq - (Default) The average distance between each leaf/taxon and every other is returned\n")
    sys.stderr.write("\tall_ave - A single value for the mean of all the per sequence averages is returned, together with the standard deviation\n")
    sys.stderr.write("\tmatrix - Rather than computing values the matrix of distances is reported\n")
    sys.stderr.write("\nOPTIONS\n\t-nominalN N - the nominal count of values, greater than the actual count\n")
    sys.stderr.write("\t-MISL N - Median Input Sequence Length (to get distance per input site)\n")
    sys.stderr.write("\t-median - Use median and median absolute deviation in place of mean and SD\n")
    sys.stderr.write("\t-sort_matrix - If matrix output requested, sort each list of nodes by ascending distance (rather than input order)\n")
    sys.stderr.write("\t-sort_matrix_descending - Sort matrix in order of descending distance from target node\n")
    sys.stderr.write("\t-include_branches - when creating matrix, include branches (as well as leaves) (default: False)\n")
    sys.stderr.write("\t-matrix_sep S - When outputing the matix, use S as the separator between Node-lenth pairs (default: blank)\n")
    sys.exit(0)

  use_median = False
  i = 1
  while i < len(args) :
    if args[i] in ["per_seq", "all_ave", "matrix"] :
      params_dict["action"] = args[i]
      i +=1
    elif args[i] == "-matrix_sep" :
      params_dict[args[i][1:]] = args[i+1]
      i +=2
    elif args[i] in ["-nominalN", "-MISL"] :
      try:
        params_dict[args[i][1:]] = int(args[i+1])
        if  params_dict[args[i][1:]] <=0 :
          raise ValueError
      except ValueError:
        sys.stderr.write("The value given for %s (%s) is not a a positive integer\n" % (args[i], args[i+1]))
        sys.exit(1)
      i += 2
    elif  args[i] in ["-median", "-sort_matrix", "-sort_matrix_descending", "-include_branches"] :
      params_dict[args[i][1:]] = True
      i += 1
    else:
      break

  infilename = args[i]
  try:
    infile = open(infilename, 'r')
  except:
    sys.stderr.write("Cannot open file %s\n" % infilename)
    sys.exit(1)

  tree = get_tree_from_file(infile)

  return(tree, params_dict)

# This has been split off from init allow it to be called separately
def get_tree_from_file(infile) :
  # The optional statement of overall likelihood breaks the parser to remove it (at least for now)
  tree_string = "".join(map(lambda x: x[:-1], infile))
  i = tree_string.find(']')
  if i > -1 :
    tree_string = tree_string[i+1:]

  # taxa beginning with a number will break the parser, so quote IDs (must have at least one
  # alphabetic character or '-'
  newstring = tree_string
  for ID in ID_PATN.findall(tree_string) :
    newstring = newstring.replace("%s:" % ID, "'%s':"%ID)
  tree_string = newstring
  if DEBUG :
    sys.stderr.write("Tree string: %s\n" % tree_string)

  try:
    tree = newick.tree.parse_tree(tree_string)
  except:
    print sys.exc_info()
    sys.stderr.write("Something went wrong with the parse of the Newick formatted string\n%s\n"\
				% tree_string)
    sys.exit(1)
  return(tree)

def map_the_nodes(tree) :
  global node_map

  if node_map == {}  : # root 
    node_map[tree] = (None, None)
  edges = tree.get_edges()
  for (branch, boot, brlen) in edges :
    if type(branch) == newick.tree.Tree :
      map_the_nodes(branch)
    node_map[branch] = (tree, brlen)
  return

def print_node(node, brlen = None, brackets = False, prepend=False, prependtype = True, return_string=False, width=None) :
  if type(node)  == newick.tree.Leaf :
    if prependtype :
      printstring = "L_%s" %  node.identifier
    else:
      printstring = "%s" %  node.identifier
  else:
    if prependtype :
      printstring = "B_%d" %  hash(node)
    else:
      printstring = "%d" %  hash(node)
  if brlen != None :
    printstring = "%s %0.4f" % (printstring, brlen)
  if brackets :
    printstring = "(%s)" % printstring
  if prepend :
    sep = params_dict["matrix_sep"]
    printstring = "%s%s" % (sep, printstring)
  if return_string:
    return(printstring)
  if width != None :
    printstring = printstring.ljust(width)
  sys.stdout.write(printstring)
  return

# node_map points up from child node (possbily leaf) to parent (always a tree)
def print_node_map(nm) :
  for node in nm.keys() :
    branch, brlen = nm[node]
    print_node(branch, brlen, False, True)
    print_node(node, prepend = True)
    sys.stdout.write("\n")

def print_matrix(matrix, params_dict) :
  for node in matrix.keys() :
    if type(node) == newick.tree.Leaf :
      print_node(node, prependtype=True)
      nodelist = map(lambda x: (x[1], x[0]),  matrix[node].items())
      if params_dict["sort_matrix"] or params_dict["sort_matrix_descending"] :
	nodelist.sort()
      if params_dict["sort_matrix_descending"] :
	nodelist.reverse()
      for dist, node1 in nodelist :
	if type(node1) == newick.tree.Leaf :
	  print_node(node1, dist, True, True)
      sys.stdout.write("\n")
  return
      
def node_map_to_matrix(node_map) :
  matrix = {}
  nodelist = node_map.keys()
  for node in nodelist :
    matrix[node] = {}
  for node in nodelist :
    parent, brlen = node_map[node]
    if parent != None :
      matrix[node][parent] = brlen
      matrix[parent][node] = brlen
  return(matrix)


"""
After node_map_to_matrix, the matrix reflects the nodes that are reachable from any given
node in one step. Now need to iterate that to complete the matrix so it shows the distances
from a given node to each of the other nodes. Will require at worst N-2 iterations for N
nodes (including internal nodes)
"""
def commplete_matrix(firstmatrix, tree) :
  matrix = {}
  for node1 in firstmatrix.keys() :
    matrix[node1] = {}
    for node2, brlen in firstmatrix[node1].items() :
      matrix[node1][node2] = brlen

  nodelist = matrix.keys() # include both leaves and internal nodes
  pending = []
  nodelistlen = len(nodelist)
  for i in range(nodelistlen -1) :
    node1 = nodelist[i]
    for j in range(i+1, nodelistlen) :
      node2 = nodelist[j]
      if matrix[node1].has_key(node2) :
	matrix[node2][node1] = matrix[node1][node2] # just in case
      elif matrix[node2].has_key(node1) :
	 matrix[node1][node2] = matrix[node2][node1]
      else:
	pending.append((node1, node2))
  for node1, node2 in pending :
    if DEBUG :
      print 'searching for node 1', 
      print_node(node1, prepend=True)
      print " node 2",
      print_node(node2, prepend=True)
      print
    dist = recursive_complete(node1, node2, matrix, firstmatrix, set([node1]))
    if DEBUG :
      print "matrix", 
      print_node(node1, prepend=True)
      print_node(node2, prepend=True)
      print " ", dist
    matrix[node1][node2] = matrix[node2][node1] = dist
  return(matrix)

"""
Was not able to get from node1 to node2. Can we find a path from the nodes reachable from
node 1
"""
def recursive_complete(node1, node2, matrix, firstmatrix, seen_set) :
  if matrix[node1].has_key(node2) :  
    matrix[node2][node1] = matrix[node1][node2] # in case not already in place
    return(matrix[node1][node2])
  if matrix[node2].has_key(node1) :  # done on the way to somewhere else
    matrix[node1][node2] = matrix[node2][node1] 
    return(matrix[node2][node1])
    
  distlist = []
  for node11 in firstmatrix[node1].keys() :
    if DEBUG :
      print ' node1', 
      print_node(node1, prepend=True)
      print ' node11',
      print_node(node11, prepend=True)
    if not node11 in seen_set:
      if DEBUG:
        print '  not yet seen'
      seen_set.add(node11)
      dist = recursive_complete(node11, node2, matrix, firstmatrix, seen_set)
      if dist != LARGE_NUMBER :
	if DEBUG :
          print 'dist', dist
        matrix[node11][node2] = matrix[node2][node11] = dist
        distlist.append(matrix[node1][node11] + dist)
    else:
      if DEBUG:
	print '  seen'
 
  if distlist == []:
    return(LARGE_NUMBER) 
  dist = min(distlist)
  matrix[node1][node2] = matrix[node2][node1] = dist
  return(dist)

# remove branches from matrix once matrix computation has completed
def delete_branches_from_matrix(matrix) :
  new_matrix = {}
  for node in matrix.keys() :
    if type(node) == newick.tree.Leaf :
      new_matrix[node] = {}
      for  node1,dist in matrix[node].items() :
	if type(node1) == newick.tree.Leaf :
	  new_matrix[node][node1] = dist
  return(new_matrix)

def tree_to_matrix(tree, params_dict):
  global node_map

  node_map = {}
  map_the_nodes(tree)
  if DEBUG :
    print 'node_map'
    print_node_map(node_map)
  matrix = node_map_to_matrix(node_map)
  if DEBUG :
    print 'matrix after node_map_to_matrix'
    print_matrix(matrix, params_dict)
  matrix = commplete_matrix(matrix, tree)
  # Branches nodes needed for matrix computation,but, generally, not later
  if not params_dict['include_branches'] :
    matrix = delete_branches_from_matrix(matrix)
  if DEBUG :
    print 'matrix after commplete_matrix'
    print_matrix(matrix, params_dict)
  return(matrix)

"""
Compute mean and SD on a list of numerical values. If a  nominal database 
size is given, it is assumed that the list of values to be analysed has been packed
out with zero values, reducing the mean. 

The previous version reduced the SD by assuming all duplicates are the mean. However
SD is not a useful measure in some cases, e.g. a symmetric star will have SD of 0
as will a similar star that is double the size. Will now produce a mean and an 
adjusted mean  if nominalN specified, otherwise mean and SD

Another alternative is to compute the median and median absolute deviation (adjusted)
"""
def stats(list, params_dict) :
  list_len = len(list)
  if params_dict["median"] :
    list.sort()
    halfway = list_len / 2
    if list_len % 2 == 0 :
       median = (list[halfway - 1] + list[halfway]) / 2.0
    else:
       median = list[halfway]

    madlist = []
    for i in list :
      madlist.append(abs(i - median))
    madlist.sort()
    if list_len % 2 == 0 :
       mad = (madlist[halfway - 1] + madlist[halfway]) / 2.0
    else:
       mad = madlist[halfway]
    mad *= MAD_MULTIPLIER
    return(median, mad)
    
  # mean and either adjusted mean or SD returned
  nominalN = params_dict["nominalN"]
  N = list_len = len(list)
  if nominalN != None :
    if nominalN >= N :
      N = nominalN
    else:
      sys.stderr.write("Warning: The nominal size of the dabase %d is less than the actual count of sequences %d\n"\
                                % (nominalN, N))

  sum = 0.0
  for i in list :
    sum += i
  mean = float(sum) / list_len

  if  nominalN == None:
    sumsqr = 0
    for i in list :
      diff = i - mean
      sumsqr += (diff * diff)
    SD_or_adj_mean = math.sqrt(sumsqr/(N-1))
  else:
    SD_or_adj_mean = float(sum) / N 
  return(mean, SD_or_adj_mean)


def leaf_dist_stats(matrix, params_dict, tree_or_leaflist):
  avelist = []
  if type(tree_or_leaflist) == newick.tree.Tree :
    leaflist = tree_or_leaflist.get_leaves()
  else:
    leaflist = tree_or_leaflist
  longest_ID_len = 0
  for leaf in leaflist :
    ID = print_node(leaf, return_string=True, prependtype=False)
    if len(ID) > longest_ID_len:
      longest_ID_len = len(ID)
  for leaf in leaflist :
    distlist = []
    for node, dist in matrix[leaf].items() :
      if type(node) ==  newick.tree.Leaf :
	distlist.append(dist)
    # "mean" may be median, and "SD" may be adjusted mean (if nominalN specified), SD or MAD
    mean, SD = stats(distlist, params_dict)  
    if params_dict["action"] == "per_seq" :
      print_node(leaf, prependtype=False, width=longest_ID_len+4)
      sys.stdout.write("%0.6f  %0.6f\n" % (mean, SD))
    avelist.append(mean)
  # will be median and MAD if -median set
  mean, adj_mean = stats(avelist, params_dict)
  N = len(avelist)
  if  params_dict["median"] :
    adj_mean_per_site100 = mean
    if params_dict["MISL"] != None:
      adj_mean_per_site100 = 100.0 * adj_mean_per_site100 / params_dict["MISL"]
  else:
    if params_dict["nominalN"] == None : # Otherwise reports SD, which is irrelevant here
      adj_mean = mean
    if params_dict["MISL"] != None:
      adj_mean_per_site100 = 100.0 * adj_mean / params_dict["MISL"]
    else:
      adj_mean_per_site100 = adj_mean
  return(mean, adj_mean, adj_mean_per_site100, N)

if __name__ == "__main__" :
  tree, params_dict = init()
  if DEBUG :
    print 'tree', tree
    print "Leaves", tree.get_leaves()
  matrix = tree_to_matrix(tree, params_dict)
  if params_dict["action"] == "matrix" :
    print_matrix(matrix, params_dict)
  else:
    mean, adj_mean, adj_mean_per_site100, N = leaf_dist_stats(matrix, params_dict, tree)
    if params_dict["nominalN"] == None:
      nominalN = N
    else:
      nominalN = params_dict["nominalN"]
    if params_dict["median"] :
      sys.stdout.write("N %d  nominalN %d  median %0.4f  MAD %0.4f  adj_median_per_site100 %0.4f\n"\
	  % (N, nominalN, mean, adj_mean, adj_mean_per_site100))
    else:
      sys.stdout.write("N %d  nominalN %d  mean %0.4f  adj_mean %0.4f  adj_mean_per_site100 %0.4f\n"\
	  % (N, nominalN, mean, adj_mean, adj_mean_per_site100))

