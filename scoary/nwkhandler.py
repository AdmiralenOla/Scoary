import sys

try:
    from ete3 import Tree
    from ete3.parser.newick import NewickError
    import six
except ImportError:
    sys.exit("ERROR: Could not import ete3. You need 'ete3' AND 'six' installed to read custom trees.")

def ReadTreeFromFile(filepath):
    """
    Uses ete3 to read a newick tree file, and converts this to a Scoary-readable nested list
    """
    try:
        myTree = Tree(filepath)
    except NewickError as e:
        sys.exit("Corrupted or non-existing custom tree file? %s" % e)
        
    myTree.resolve_polytomy(recursive=True)
    myTreeList, members = RecTree2List(myTree,Members=None)

    return myTreeList, members

def RecTree2List(Tree, Members=None):
    """
    Recursive function that at each node create a list of the children nodes. Can be nested
    """
    List = []
    if Members is None:
        Members = []

    if len(Tree._children) == 0:
        shavedname = Tree.name.lstrip("'\"").rstrip("'\"")
        Members += [shavedname]
        return str(shavedname), Members
    else:
        for node in Tree._children:
            mynode, Members = RecTree2List(node, Members)
            List.append(mynode)
        return List, Members
