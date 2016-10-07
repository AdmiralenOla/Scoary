import sys

# Note: The Matrix and QuadTree implementations are heavily based on original
# implementations by Christian Storm Pedersen.

# Python 2/3 annoyances
try:
    xrange
except NameError:
    xrange = range


class Matrix:
    """
    A matrix stored as a list of lists
    """
    def __init__(self, dim, elm=sys.maxsize):
        """
        Builds empty nxn matrix
        """
        self.undef = sys.maxsize
        self.dim = int(dim)
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
        self.undef = sys.maxsize
        self.dim = dim
        n = self.dim + self.dim % 2
        if names is None:
            names = [str(x) for x in xrange(dim)]
        self.names = names
        self.level = []
        while n > 1:
            n += (n % 2)
            self.level.append(Matrix(n))
            n = (n+1) // 2

    def get_elm(self, i, j):
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
                    next_row.append(self.quad_min(i, j, l))
            i //= 2
            curr_row = next_row

    def insert_col(self, j, col):
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
                    next_col.append(self.quad_min(i, j, l))
            j //= 2
            curr_col = next_col

    def min(self):
        """
        Returns minimum element stored in tree
        """
        return self.quad_min(0, 0, self.level[-1])

    def argmin(self):
        """
        Returns coordinates of minimum element in tree
        """
        i = j = 0
        for l in reversed(self.level[1:]):
            i, j = self.quad_argmin(i, j, l)
            i *= 2
            j *= 2
        return self.quad_argmin(i, j, self.level[0])

    def quad_min_all(self, i, j, l):
        """
        Returns the minimum element stored in the quad (i,j) and its coordinates
        """
        # Need even numbers
        i = (i//2) * 2
        j = (j//2) * 2
        return min((l[i][j], i, j),
                   (l[i+1][j], i+1, j),
                   (l[i][j+1], i, j+1),
                   (l[i+1][j+1], i+1, j+1))

    def quad_min(self, i, j, l):
        """
        Returns the minimum element stored in the quad containing (i,j)
        """
        return self.quad_min_all(i, j, l)[0]

    def quad_argmin(self, i, j, l):
        """
        Returns the coordinates of the minimum element in the quad containing (i,j)
        """
        return self.quad_min_all(i, j, l)[1:]


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
    def __init__(self, leftnode, rightnode, GTC, OR):
        """
        Constructs a phylotree and links it to its left and right nodes
        """
        # First check if left and right are tips (type = str). Then we need to look up the gene-trait combination value at that node.
        # Also check if left and right nodes are PhyloTrees. If not, recursively create PhyloTrees.

        if type(leftnode) is str:
            self.leftnode = Tip(GTC[leftnode])
        elif isinstance(leftnode, PhyloTree):
            self.leftnode = leftnode
        else:
            self.leftnode = PhyloTree(leftnode=leftnode[0], rightnode=leftnode[1], GTC=GTC, OR=OR)

        if type(rightnode) is str:
            self.rightnode = Tip(GTC[rightnode])
        elif isinstance(rightnode, PhyloTree):
            self.rightnode = rightnode
        else:
            self.rightnode = PhyloTree(leftnode=rightnode[0], rightnode=rightnode[1], GTC=GTC, OR=OR)

        # Initialize the max number of paths. Set to -1 meaning they cannot be reached
        self.OR = OR
        self.maxvalues = {"AB": -1, "Ab": -1, "aB": -1, "ab": -1, "0": -1}
        self.max_propairs = {"AB": -1, "Ab": -1, "aB": -1, "ab": -1, "0": -1}
        self.max_antipairs = {"AB": -1, "Ab": -1, "aB": -1, "ab": -1, "0": -1}
        self.calculate_max()
        self.max_contrasting_pairs = max(self.maxvalues.values())
        self.max_contrasting_propairs = max(self.max_propairs.values())
        self.max_contrasting_antipairs = max(self.max_antipairs.values())

    def calculate_max(self):
        """
        A method for calculating the max number of pairings under the 5 conditions
        """
        for condition in ["AB", "Ab", "aB", "ab", "nofree"]:
            if condition in ["AB", "Ab", "aB", "ab"]:
                pairings = self.calculate_max_condition(condition)
                self.maxvalues[condition] = pairings["Total"]
                self.max_propairs[condition] = pairings["Pro"]
                self.max_antipairs[condition] = pairings["Anti"]
            else:  # Condition == nofree
                pairings = self.calculate_max_nofree()
                self.maxvalues["0"] = pairings["Total"]
                self.max_propairs["0"] = pairings["Pro"]
                self.max_antipairs["0"] = pairings["Anti"]

    def calculate_max_condition(self, condition):
        """
        When passed for example 'AB', finds out the 9 distinct conditions and calculates the max
        """
        Possible_conditions = set(["AB", "Ab", "aB", "ab"])
        Possible_conditions.remove(condition)
        Otherconditions = list(Possible_conditions)  # Now we have a list of the elements that are NOT condition
        max_pairs_1 = -1
        max_pairs_2 = -1
        max_pairs_3 = -1
        max_pairs_4 = -1
        max_pairs_5 = -1
        max_pairs_6 = -1
        max_pairs_7 = -1
        max_pairs_8 = -1
        max_pairs_9 = -1

        max_propairs_1 = -1
        max_propairs_2 = -1
        max_propairs_3 = -1
        max_propairs_4 = -1
        max_propairs_5 = -1
        max_propairs_6 = -1
        max_propairs_7 = -1
        max_propairs_8 = -1
        max_propairs_9 = -1

        max_antipairs_1 = -1
        max_antipairs_2 = -1
        max_antipairs_3 = -1
        max_antipairs_4 = -1
        max_antipairs_5 = -1
        max_antipairs_6 = -1
        max_antipairs_7 = -1
        max_antipairs_8 = -1
        max_antipairs_9 = -1

        if self.leftnode.maxvalues[condition] > -1 and self.rightnode.maxvalues["0"] > -1:
            max_pairs_1 = self.leftnode.maxvalues[condition] + self.rightnode.maxvalues["0"]
            max_propairs_1 = self.leftnode.max_propairs[condition] + self.rightnode.max_propairs["0"]
            max_antipairs_1 = self.leftnode.max_antipairs[condition] + self.rightnode.max_antipairs["0"]

        if self.leftnode.maxvalues[condition] > -1 and self.rightnode.maxvalues[Otherconditions[0]] > -1:
            max_pairs_2 = self.leftnode.maxvalues[condition] + self.rightnode.maxvalues[Otherconditions[0]]
            max_propairs_2 = self.leftnode.max_propairs[condition] + self.rightnode.max_propairs[Otherconditions[0]]
            max_antipairs_2 = self.leftnode.max_antipairs[condition] + self.rightnode.max_antipairs[Otherconditions[0]]

        if self.leftnode.maxvalues[condition] > -1 and self.rightnode.maxvalues[Otherconditions[1]] > -1:
            max_pairs_3 = self.leftnode.maxvalues[condition] + self.rightnode.maxvalues[Otherconditions[1]]
            max_propairs_3 = self.leftnode.max_propairs[condition] + self.rightnode.max_propairs[Otherconditions[1]]
            max_antipairs_3 = self.leftnode.max_antipairs[condition] + self.rightnode.max_antipairs[Otherconditions[1]]

        if self.leftnode.maxvalues[condition] > -1 and self.rightnode.maxvalues[Otherconditions[2]] > -1:
            max_pairs_4 = self.leftnode.maxvalues[condition] + self.rightnode.maxvalues[Otherconditions[2]]
            max_propairs_4 = self.leftnode.max_propairs[condition] + self.rightnode.max_propairs[Otherconditions[2]]
            max_antipairs_4 = self.leftnode.max_antipairs[condition] + self.rightnode.max_antipairs[Otherconditions[2]]

        if self.leftnode.maxvalues[condition] > -1 and self.rightnode.maxvalues[condition] > -1:
            max_pairs_5 = self.leftnode.maxvalues[condition] + self.rightnode.maxvalues[condition]
            max_propairs_5 = self.leftnode.max_propairs[condition] + self.rightnode.max_propairs[condition]
            max_antipairs_5 = self.leftnode.max_antipairs[condition] + self.rightnode.max_antipairs[condition]

        if self.leftnode.maxvalues["0"] > -1 and self.rightnode.maxvalues[condition] > -1:
            max_pairs_6 = self.leftnode.maxvalues["0"] + self.rightnode.maxvalues[condition]
            max_propairs_6 = self.leftnode.max_propairs["0"] + self.rightnode.max_propairs[condition]
            max_antipairs_6 = self.leftnode.max_antipairs["0"] + self.rightnode.max_antipairs[condition]

        if self.leftnode.maxvalues[Otherconditions[0]] > -1 and self.rightnode.maxvalues[condition] > -1:
            max_pairs_7 = self.leftnode.maxvalues[Otherconditions[0]] + self.rightnode.maxvalues[condition]
            max_propairs_7 = self.leftnode.max_propairs[Otherconditions[0]] + self.rightnode.max_propairs[condition]
            max_antipairs_7 = self.leftnode.max_antipairs[Otherconditions[0]] + self.rightnode.max_antipairs[condition]

        if self.leftnode.maxvalues[Otherconditions[1]] > -1 and self.rightnode.maxvalues[condition] > -1:
            max_pairs_8 = self.leftnode.maxvalues[Otherconditions[1]] + self.rightnode.maxvalues[condition]
            max_propairs_8 = self.leftnode.max_propairs[Otherconditions[1]] + self.rightnode.max_propairs[condition]
            max_antipairs_8 = self.leftnode.max_antipairs[Otherconditions[1]] + self.rightnode.max_antipairs[condition]

        if self.leftnode.maxvalues[Otherconditions[2]] > -1 and self.rightnode.maxvalues[condition] > -1:
            max_pairs_9 = self.leftnode.maxvalues[Otherconditions[2]] + self.rightnode.maxvalues[condition]
            max_propairs_9 = self.leftnode.max_propairs[Otherconditions[2]] + self.rightnode.max_propairs[condition]
            max_antipairs_9 = self.leftnode.max_antipairs[Otherconditions[2]] + self.rightnode.max_antipairs[condition]

        max_pairs = max(max_pairs_1, max_pairs_2, max_pairs_3, max_pairs_4, max_pairs_5, max_pairs_6, max_pairs_7, max_pairs_8, max_pairs_9)

        # Calculate maximum number of propairs, given a maxmimum number of pairs
        max_propairs = -1
        if max_pairs == max_pairs_1:
            max_propairs = max(max_propairs, max_propairs_1)
        if max_pairs == max_pairs_2:
            max_propairs = max(max_propairs, max_propairs_2)
        if max_pairs == max_pairs_3:
            max_propairs = max(max_propairs, max_propairs_3)
        if max_pairs == max_pairs_4:
            max_propairs = max(max_propairs, max_propairs_4)
        if max_pairs == max_pairs_5:
            max_propairs = max(max_propairs, max_propairs_5)
        if max_pairs == max_pairs_6:
            max_propairs = max(max_propairs, max_propairs_6)
        if max_pairs == max_pairs_7:
            max_propairs = max(max_propairs, max_propairs_7)
        if max_pairs == max_pairs_8:
            max_propairs = max(max_propairs, max_propairs_8)
        if max_pairs == max_pairs_9:
            max_propairs = max(max_propairs, max_propairs_9)

        # Calculate maximum number of antipairs, given a maxmimum number of pairs
        max_antipairs = -1
        if max_pairs == max_pairs_1:
            max_antipairs = max(max_antipairs, max_antipairs_1)
        if max_pairs == max_pairs_2:
            max_antipairs = max(max_antipairs, max_antipairs_2)
        if max_pairs == max_pairs_3:
            max_antipairs = max(max_antipairs, max_antipairs_3)
        if max_pairs == max_pairs_4:
            max_antipairs = max(max_antipairs, max_antipairs_4)
        if max_pairs == max_pairs_5:
            max_antipairs = max(max_antipairs, max_antipairs_5)
        if max_pairs == max_pairs_6:
            max_antipairs = max(max_antipairs, max_antipairs_6)
        if max_pairs == max_pairs_7:
            max_antipairs = max(max_antipairs, max_antipairs_7)
        if max_pairs == max_pairs_8:
            max_antipairs = max(max_antipairs, max_antipairs_8)
        if max_pairs == max_pairs_9:
            max_antipairs = max(max_antipairs, max_antipairs_9)

        return {"Total": max_pairs, "Pro": max_propairs, "Anti": max_antipairs}

    def calculate_max_nofree(self):
        """
        Under the condition of no free paths, only 5 distinct possibilities exits
        """
        # No free pairs in either:
        max_pairs_nofree = -1
        max_pairs_1100 = -1
        max_pairs_0011 = -1
        max_pairs_1001 = -1
        max_pairs_0110 = -1

        max_propairs_nofree = -1
        max_propairs_1100 = -1
        max_propairs_0011 = -1
        max_propairs_1001 = -1
        max_propairs_0110 = -1

        max_antipairs_nofree = -1
        max_antipairs_1100 = -1
        max_antipairs_0011 = -1
        max_antipairs_1001 = -1
        max_antipairs_0110 = -1

        if self.leftnode.maxvalues["0"] > -1 and self.rightnode.maxvalues["0"] > -1:
            max_pairs_nofree = self.leftnode.maxvalues["0"] + self.rightnode.maxvalues["0"]
            max_propairs_nofree = self.leftnode.max_propairs["0"] + self.rightnode.max_propairs["0"]
            max_antipairs_nofree = self.leftnode.max_antipairs["0"] + self.rightnode.max_antipairs["0"]

        if self.leftnode.maxvalues["AB"] > -1 and self.rightnode.maxvalues["ab"] > -1:
            max_pairs_1100 = self.leftnode.maxvalues["AB"] + self.rightnode.maxvalues["ab"] + 1
            max_propairs_1100 = self.leftnode.max_propairs["AB"] + self.rightnode.max_propairs["ab"]
            max_antipairs_1100 = self.leftnode.max_antipairs["AB"] + self.rightnode.max_antipairs["ab"]
            #if self.OR < 1:
            #    max_antipairs_1100 += 1
            #else:
            #    max_propairs_1100 += 1
            max_propairs_1100 += 1

        if self.leftnode.maxvalues["ab"] > -1 and self.rightnode.maxvalues["AB"] > -1:
            max_pairs_0011 = self.leftnode.maxvalues["ab"] + self.rightnode.maxvalues["AB"] + 1
            max_propairs_0011 = self.leftnode.max_propairs["ab"] + self.rightnode.max_propairs["AB"]
            max_antipairs_0011 = self.leftnode.max_antipairs["ab"] + self.rightnode.max_antipairs["AB"]
            #if self.OR < 1:
            #    max_antipairs_0011 += 1
            #else:
            #    max_propairs_0011 += 1
            max_propairs_0011 += 1

        if self.leftnode.maxvalues["Ab"] > -1 and self.rightnode.maxvalues["aB"] > -1:
            max_pairs_1001 = self.leftnode.maxvalues["Ab"] + self.rightnode.maxvalues["aB"] + 1
            max_propairs_1001 = self.leftnode.max_propairs["Ab"] + self.rightnode.max_propairs["aB"]
            max_antipairs_1001 = self.leftnode.max_antipairs["Ab"] + self.rightnode.max_antipairs["aB"]# + 1
            #if self.OR < 1:
            #    max_propairs_1001 += 1
            #else:
            #    max_antipairs_1001 += 1
            max_antipairs_1001 += 1

        if self.leftnode.maxvalues["aB"] > -1 and self.rightnode.maxvalues["Ab"] > -1:
            max_pairs_0110 = self.leftnode.maxvalues["aB"] + self.rightnode.maxvalues["Ab"] + 1
            max_propairs_0110 = self.leftnode.max_propairs["aB"] + self.rightnode.max_propairs["Ab"]
            max_antipairs_0110 = self.leftnode.max_antipairs["aB"] + self.rightnode.max_antipairs["Ab"]
            #if self.OR < 1:
            #    max_propairs_0110 += 1
            #else:
            #    max_antipairs_0110 += 1
            max_antipairs_0110 += 1

        max_pairs = max(max_pairs_nofree, max_pairs_1100, max_pairs_0011, max_pairs_1001, max_pairs_0110)

        # Calculate max number of propairs
        max_propairs = -1  # Max_propairs can never go below -1
        if max_pairs == max_pairs_nofree:
            max_propairs = max(max_propairs, max_propairs_nofree)
        if max_pairs == max_pairs_1100:
            max_propairs = max(max_propairs, max_propairs_1100)
        if max_pairs == max_pairs_0011:
            max_propairs = max(max_propairs, max_propairs_0011)
        if max_pairs == max_pairs_1001:
            max_propairs = max(max_propairs, max_propairs_1001)
        if max_pairs == max_pairs_0110:
            max_propairs = max(max_propairs, max_propairs_0110)

        # Calculate max number of antipairs
        max_antipairs = -1  # Max_antipairs can never go below -1
        if max_pairs == max_pairs_nofree:
            max_antipairs = max(max_antipairs, max_antipairs_nofree)
        if max_pairs == max_pairs_1100:
            max_antipairs = max(max_antipairs, max_antipairs_1100)
        if max_pairs == max_pairs_0011:
            max_antipairs = max(max_antipairs, max_antipairs_0011)
        if max_pairs == max_pairs_1001:
            max_antipairs = max(max_antipairs, max_antipairs_1001)
        if max_pairs == max_pairs_0110:
            max_antipairs = max(max_antipairs, max_antipairs_0110)

        return {"Total": max_pairs, "Pro": max_propairs, "Anti": max_antipairs}


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
        for condition in ["AB", "Ab", "aB", "ab", "0"]:
            if condition == tipvalue:
                self.maxvalues[condition] = 0
            else:
                self.maxvalues[condition] = -1
        self.max_propairs = {k: v for (k, v) in self.maxvalues.items()}
        self.max_antipairs = {k: v for (k, v) in self.maxvalues.items()}

if __name__ == '__main__':
    pass
