
import sys

# Note: The Matrix and QuadTree implementations are heavily based on original implementations by Christian Storm Pedersen.

class Matrix:
	"""
	A matrix stored as a list of lists
	"""
	def __init__(self,dim,elm=sys.maxint):
		"""
		Builds empty nxn matrix
		"""
		self.undef = sys.maxint
		self.dim = dim
		self.data = []
		for i in xrange(self.dim):
			self.data.append(self.dim * [elm])
			
	def __getitem__(self, i):
		"""
		Get row i from Matrix
		"""
		return self.data[i]
				
	def __str__(self):
		s = ""
		for i in xrange(self.dim):
			for j in xrange(self.dim):
				s += str(self.data[i][j]) + " "
			s += "\n"
		return s.strip("\n")

class QuadTree:
	"""
	A basic QuadTree with names
	"""
	def __init__(self, dim, names=None):
		"""
		Constructs a quad tree of dimension dim and fills it with 0's
		"""
		self.undef = sys.maxint
		self.dim = dim
		n = self.dim + self.dim % 2
		if names is None:
			names = [str(x) for x in xrange(dim)]
		self.names = names
		self.level = []
		while n > 1:
			n += (n % 2)
			self.level.append(Matrix(n))
			n = (n+1) / 2
			
	def get_elm(self,i,j):
		"""
		Returns the element at position (i,j) in the quad tree
		"""
		return self.level[0][i][j]
		
	def insert_row(self, i, row):
		"""
		Inserts row (of dim elements) as row number i
		"""
		curr_row = row
		for l in self.level:
			if len(curr_row) % 2 == 1:
				curr_row.append(self.undef)
			next_row = []
			for j in xrange(len(curr_row)):
				l[i][j] = curr_row[j]
				if j % 2 == 1:
					next_row.append(self.quad_min(i,j,l))
			i /= 2
			curr_row = next_row
			
	def insert_col(self,j,col):
		"""
		Inserts col (of dim elements) as col number j
		"""
		curr_col = col
		for l in self.level:
			if len(curr_col) % 2 == 1:
				curr_col.append(self.undef)
			next_col = []
			for i in xrange(len(curr_col)):
				l[i][j] = curr_col[i]
				if i % 2 == 1:
					next_col.append(self.quad_min(i,j,l))
			j /= 2
			curr_col = next_col
		
	def min(self):
		"""
		Returns minimum element stored in tree
		"""
		return self.quad_min(0,0,self.level[-1])
		
	def argmin(self):
		"""
		Returns coordinates of minimum element in tree
		"""
		i = j = 0
		for l in reversed(self.level[1:]):
			i, j = self.quad_argmin(i,j,l)
			i *= 2
			j *= 2
		return self.quad_argmin(i,j,self.level[0])
		
	def quad_min_all(self, i, j, l):
		"""
		Returns the minimum element stored in the quad (i,j) and its coordinates
		"""
		# Need even numbers
		i = (i/2) * 2
		j = (j/2) * 2
		return min((l[i][j],i,j), (l[i+1][j],i+1,j), (l[i][j+1],i,j+1), (l[i+1][j+1],i+1,j+1))
		
	def quad_min(self,i,j,l):
		"""
		Returns the minimum element stored in the quad containing (i,j)
		"""
		return self.quad_min_all(i,j,l)[0]
		
	def quad_argmin(self,i,j,l):
		"""
		Returns the coordinates of the minimum element in the quad containing (i,j)
		"""
		return self.quad_min_all(i,j,l)[1:]
			
class PhyloTree:
	"""
	A class that represents a binary tree. Phylotrees can be nested. They can also contain tips at left, right or both nodes.
	Stores the max number of paths under each of the following 5 conditions:
	1. Free path to AB. 
	2. Free path to Ab. 
	3. Free path to aB. 
	4. Free path to ab. 
	5. No free path.
	"""
	def __init__(self, leftnode, rightnode, GTC):
		"""
		Constructs a phylotree and links it to its left and right nodes
		"""
		# First check if left and right are tips (single-element list). Then we need to look up the gene-trait combination value at that node.
		# Also check if left and right nodes are PhyloTrees. If not, recursively create PhyloTrees.
		if len(leftnode) == 1:
			self.leftnode = Tip(GTC[leftnode[0]])
		elif isinstance(leftnode, PhyloTree):
			self.leftnode = leftnode
		else:
			self.leftnode = PhyloTree(leftnode=leftnode[0], rightnode=leftnode[1], GTC=GTC)
		if len(rightnode) == 1:
			self.rightnode = Tip(GTC[rightnode[0]])
		elif isinstance(rightnode, PhyloTree):
			self.rightnode = rightnode
		else:
			self.rightnode = PhyloTree(leftnode=rightnode[0], rightnode=rightnode[1], GTC=GTC)
		
		# Initialize the max number of paths. Set to -1 meaning they cannot be reached
		self.maxvalues = {"AB": -1, "Ab": -1, "aB": -1, "ab": -1, "0": -1}
		self.calculate_max()
		self.max_contrasting_pairs = max(self.maxvalues.values())
		
	def calculate_max(self):
		"""
		A method for calculating the max number of pairings under the 5 conditions
		"""
		for condition in ["AB","Ab","aB","ab","nofree"]:
			if condition in ["AB", "Ab", "aB", "ab"]:
				self.maxvalues[condition] = self.calculate_max_condition(condition)
			else: # Condition == nofree
				self.maxvalues["0"] = self.calculate_max_nofree()
				
		
	def calculate_max_condition(self,condition):
		"""
		When passed for example 'AB', finds out the 9 distinct conditions and calculates the max
		"""
		Possible_conditions = set(["AB", "Ab", "aB", "ab"])
		Possible_conditions.remove(condition)
		Otherconditions = list(Possible_conditions) # Now we have a list of the elements that are NOT condition
		max_pairs = -1
		if self.leftnode.maxvalues[condition] > -1 and self.rightnode.maxvalues["0"] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues[condition] + self.rightnode.maxvalues["0"] )
		if self.leftnode.maxvalues[condition] > -1 and self.rightnode.maxvalues[Otherconditions[0]] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues[condition] + self.rightnode.maxvalues[Otherconditions[0]] )
		if self.leftnode.maxvalues[condition] > -1 and self.rightnode.maxvalues[Otherconditions[1]] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues[condition] + self.rightnode.maxvalues[Otherconditions[1]] )
		if self.leftnode.maxvalues[condition] > -1 and self.rightnode.maxvalues[Otherconditions[2]] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues[condition] + self.rightnode.maxvalues[Otherconditions[2]] )
		if self.leftnode.maxvalues[condition] > -1 and self.rightnode.maxvalues[condition] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues[condition] + self.rightnode.maxvalues[condition] )
		if self.leftnode.maxvalues["0"] > -1 and self.rightnode.maxvalues[condition] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues["0"] + self.rightnode.maxvalues[condition] )
		if self.leftnode.maxvalues[Otherconditions[0]] > -1 and self.rightnode.maxvalues[condition] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues[Otherconditions[0]] + self.rightnode.maxvalues[condition] )
		if self.leftnode.maxvalues[Otherconditions[1]] > -1 and self.rightnode.maxvalues[condition] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues[Otherconditions[1]] + self.rightnode.maxvalues[condition] )
		if self.leftnode.maxvalues[Otherconditions[2]] > -1 and self.rightnode.maxvalues[condition] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues[Otherconditions[2]] + self.rightnode.maxvalues[condition] )
			
		return max_pairs
		
	def calculate_max_nofree(self):
		"""
		Under the condition of no free paths, only 5 distinct possibilities exits
		"""
		max_pairs = -1
		# No free pairs in either:
		if self.leftnode.maxvalues["0"] > -1 and self.rightnode.maxvalues["0"] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues["0"] + self.rightnode.maxvalues["0"] )
		if self.leftnode.maxvalues["AB"] > -1 and self.rightnode.maxvalues["ab"] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues["AB"] + self.rightnode.maxvalues["ab"] + 1 )
		if self.leftnode.maxvalues["ab"] > -1 and self.rightnode.maxvalues["AB"] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues["ab"] + self.rightnode.maxvalues["AB"] + 1 )
		if self.leftnode.maxvalues["Ab"] > -1 and self.rightnode.maxvalues["aB"] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues["Ab"] + self.rightnode.maxvalues["aB"] + 1 )
		if self.leftnode.maxvalues["aB"] > -1 and self.rightnode.maxvalues["Ab"] > -1:
			max_pairs = max( max_pairs, self.leftnode.maxvalues["aB"] + self.rightnode.maxvalues["Ab"] + 1 )

		return max_pairs
	
class Tip:
	"""
	A class that references a single tip, which can only be AB, Ab, aB or ab
	"""
	def __init__(self, tipvalue):
		"""
		Sets up the tip
		"""
		self.tipvalue = tipvalue
		self.maxvalues = {}
		for condition in ["AB","Ab","aB","ab","0"]:
			if condition == tipvalue:
				self.maxvalues[condition] = 0
			else:
				self.maxvalues[condition] = -1
