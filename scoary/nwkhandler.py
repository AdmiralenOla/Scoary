try:
    from ete3 import Tree
    from ete3.parser.newick import NewickError
except ImportError:
    sys.exit("ERROR: Could not import ete3. You need ete3 installed to read custom trees.")

def ReadTreeFromFile(filepath):
    """
    Uses ete3 to read a newick tree file, and converts this to a Scoary-readable nested list
    """
    try:
        myTree = Tree(filepath)
    except NewickError as e:
        sys.exit("Corrupted or non-existing custom tree file? %s" % e)
        
    myTree.resolve_polytomy(recursive=True)
    myTreeList = RecTree2List(myTree)

    return myTreeList

def RecTree2List(Tree):
    """
    Recursive function that at each node create a list of the children nodes. Can be nested
    """
    List = []
    if len(Tree._children) == 0:
        shavedname = Tree.name.lstrip("'\"").rstrip("'\"")
        return str(shavedname)
    else:
        for node in Tree._children:
            List.append(RecTree2List(node))
        return List
