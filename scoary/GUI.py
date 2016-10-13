#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Scoary - Microbial Pan-GWAS GUI
By Ola Brynildsrud
Norwegian Institute of Public Health
olbb@fhi.no
"""

import sys, os
import threading

try:
    import Tkinter
except ImportError:
    #Python 3 issues
    try:
        import tkinter as Tkinter
    except:
        sys.exit("Need to have Tkinter / tkinter installed")
    
try:
    import tkFileDialog
except ImportError:
    # Python 3 issues
    try:
        from tkinter import filedialog
        tkFileDialog = filedialog
    except:
        sys.exit("Could not find tkFileDialog / filedialog")

try:
    import ttk
except ImportError:
    # Python 3 issues
    try:
        from tkinter import ttk as ttk
    except:
        sys.exit("Could not find ttk / tkinter.ttk")
    
try:
    ttk
    Tkinter
except NameError:
    sys.exit("Need the following installed: Tkinter, tkFileDialog, ttk")
    
try:
    import os
except ImportError:
    sys.exit("Need to be able to import os")
    
try:
    import argparse
except ImportError:
    sys.exit("Need to be able to import argparse")
    
try:
    import scoary.methods as sm
    import scoary
except ImportError:
    sys.exit("Could not find the main Scoary executable") 

from pkg_resources import resource_string, resource_filename

class ScoaryGUI(Tkinter.Tk):
    """
    Create the main GUI window
    """
    def __init__(self,parent):
        Tkinter.Tk.__init__(self,parent)
        self.parent = parent

        self.Scoary_parameters = {"GPA": None, "Trait":None,"Tree":None,"Restrict":None,"Writetree": False, "Delimiter":",", 
                                "Startcol": "15", "Maxhits":None,"Notime":False, "Outdir":None, "Permutations":0,
                                "Cutoffs": {"I": [1,0.05], "B": [0,1.0], "BH": [0,1.0], "PW": [0,1.0], "EPW": [1,0.05], "P":[0,1.0]},                                
                                }
        self.initialize_menu()
        
        self.toppart = Tkinter.Frame(self,height="350",width="800")
        self.bottompart = Tkinter.Frame(self,height="50",width="800")
        self.toppart.pack(side='top',expand=False)
        self.bottompart.pack(side='bottom',expand=True,fill='both')
        
        self.nwpart = Tkinter.Frame(self.toppart,height="350",width="250")
        self.nepart = Tkinter.LabelFrame(self.toppart,height="350",width="550",text="Control panel")
        self.nwpart.pack(side='left',expand=False)
        self.nepart.pack(side='right',expand=True)
        
        # Add further frames to top or bottom
        
        self.logocanvas = Tkinter.Canvas(self.nwpart,height="250",width="250",relief="ridge",bd=0,highlightthickness=0)
        self.logocanvas.pack(side='top',expand=False)
        
        self.citationframe = Tkinter.Frame(self.nwpart, height="100",width="250")
        self.citationframe.pack(side='bottom',expand=False,fill='none')
        
        self.initialize_controlboard()
        self.initialize_logo()
        self.initialize_citation()
        self.initialize_statusframe()
        
        ######## MENUS ########
        
    def initialize_menu(self):
        """
        Initialize the menu at the top
        """
        self.menubar = Tkinter.Menu(self,relief="flat")
        filemenu = Tkinter.Menu(self.menubar,tearoff=0)
        filemenu.add_command(label="About",command=self.AboutScoary)
        filemenu.add_separator()
        filemenu.add_command(label="Quit",command=self.quit)
        self.filemenu = filemenu
        self.menubar.add_cascade(label="File",menu=self.filemenu)
        
        optsmenu = Tkinter.Menu(self.menubar,tearoff=0)
        optsmenu.add_command(label="Clear", command=self.ClearAll)
        optsmenu.add_command(label="Test example", command=self.TestExample)
        self.optsmenu = optsmenu
        self.menubar.add_cascade(label="Options",menu=self.optsmenu)
        
        self.config(menu=self.menubar)        
        
    def initialize_citation(self):
        """
        Initialize the citation frame - below the logo
        """
        myfontstyle = ("Arial",8)
        self.citation = Tkinter.Label(self.citationframe, text=self.citationtext(),anchor='center',justify='center', font=myfontstyle)
        self.citation.pack(expand=True)
        
    def initialize_statusframe(self):
        """
        Initialize the frame and statusbar occupying the bottom
        """
        frame = self.bottompart
        
        frame.pb = ttk.Progressbar(frame,orient='horizontal',mode='determinate',maximum=100)
        frame.pb.pack(fill='both',expand=True,side='top')

        frame.lab = Tkinter.Label(frame,text=u"Awaiting input options")
        frame.lab.pack(in_=frame.pb,expand=True)
        sys.stdout = StdoutToLabel(frame.lab, progressbar=frame.pb)
    
    def initialize_controlboard(self):
        """
        Initialize the controlboard - where all the settings are set.
        """
        board = self.nepart
        
        self.GPAentryVariable = Tkinter.StringVar()
        board.GPAentry = Tkinter.Entry(board, textvariable = self.GPAentryVariable,width=60)
        board.GPAentry.grid(column=0,row=0,sticky='W',columnspan=2)
        self.GPAentryVariable.set("Path to gene presence absence file")
        
        self.TraitsentryVariable = Tkinter.StringVar()
        board.Traitsentry = Tkinter.Entry(board, textvariable=self.TraitsentryVariable,width=60)
        board.Traitsentry.grid(column=0,row=1,sticky='W',columnspan=2)
        self.TraitsentryVariable.set("Path to traits/phenotype file")
        
        self.TreeentryVariable = Tkinter.StringVar()
        board.Treeentry = Tkinter.Entry(board, textvariable=self.TreeentryVariable,width=60)
        board.Treeentry.grid(column=0,row=2,sticky='W',columnspan=2)
        self.TreeentryVariable.set("(Optional) Path to custom tree file")
        
        self.RestrictVariable = Tkinter.StringVar()
        board.Restrictentry = Tkinter.Entry(board, textvariable=self.RestrictVariable,width=60)
        board.Restrictentry.grid(column=0,row=3,sticky='W',columnspan=2)
        self.RestrictVariable.set("(Optional) Path to file naming isolates to include")
        
        self.Outputdir = Tkinter.StringVar()
        board.Outputentry = Tkinter.Entry(board, textvariable=self.Outputdir,width=60)
        board.Outputentry.grid(column=0,row=4,sticky='W',columnspan=2)
        self.Outputdir.set("(Optional) Output directory")

        browsebuttonGPA = Tkinter.Button(board, text=u"Browse...", command=self.BrowseButtonClickGPA)
        browsebuttonGPA.grid(column=2,row=0,sticky='e')
        
        browsebuttonTraits = Tkinter.Button(board,text=u"Browse...", command=self.BrowseButtonClickTraits)
        browsebuttonTraits.grid(column=2,row=1,sticky='e')
        
        browsebuttonTreeFile = Tkinter.Button(board,text=u"Browse...", command=self.BrowseButtonClickTreeFile)
        browsebuttonTreeFile.grid(column=2,row=2,sticky='e')
        
        browsebuttonRestrict = Tkinter.Button(board,text=u"Browse...",command=self.BrowseButtonClickRestrict)
        browsebuttonRestrict.grid(column=2,row=3,sticky='e')
        
        browsebuttonOutput = Tkinter.Button(board,text=u"Browse...",command=self.BrowseButtonClickOutput)
        browsebuttonOutput.grid(column=2,row=4,sticky='e')
        
        # Initialize frame for cutoffs
        board.pframe = Tkinter.LabelFrame(board,text="Cut-offs",relief='ridge')
        board.pframe.grid(column=0,row=5,sticky='w')
        self.initialize_pvalueframe()
        
        # Initialize frame for misc options
        board.miscframe = Tkinter.LabelFrame(board,text="Misc options",relief='ridge')
        board.miscframe.grid(column=1,row=5,sticky='e',columnspan=2)
        self.initialize_miscopts()
        
        ## Create extra space
        
        board.emptyspace = Tkinter.Frame(board)
        board.emptyspace.grid(column=0,row=6,columnspan=3,pady=(20,20))
        
        masterbuttonfont = ("Courier", 16)
        
        RunButton = Tkinter.Button(board.emptyspace,text=u"Run analysis",font=masterbuttonfont,command=self.RunAnalysis,padx=15,pady=15)
        RunButton.grid(column=2,row=0)
        
        QuitButton = Tkinter.Button(board.emptyspace,text=u"Quit",font=masterbuttonfont, command=self.quit,padx=15,pady=15)
        QuitButton.grid(column=1,row=0)
        
        HelpButton = Tkinter.Button(board.emptyspace,text=u"Help",font=masterbuttonfont,command=self.HelpButton,padx=15,pady=15)
        HelpButton.grid(column=0,row=0)
        
    def initialize_miscopts(self):
        """
        Initialize the miscellaneous options
        """
        board = self.nepart
        miscframe = board.miscframe
        
        # Max hits
        self.mhtext = Tkinter.StringVar()
        miscframe.mhlab = Tkinter.Label(miscframe,textvariable=self.mhtext)
        miscframe.mhlab.grid(column=0,row=0,sticky='w')
        self.mhtext.set("Max hits")
        
        self.maxhitsvar = Tkinter.StringVar()
        miscframe.maxhitsentry = Tkinter.Entry(miscframe,textvariable=self.maxhitsvar,width=8)
        miscframe.maxhitsentry.grid(column=1,row=0)
        self.maxhitsvar.set("")
        
        # Delimiter
        self.delimtext = Tkinter.StringVar()
        miscframe.delim = Tkinter.Label(miscframe,textvariable=self.delimtext)
        miscframe.delim.grid(column=0,row=1,sticky='w')
        self.delimtext.set("Delimiter")
        
        self.delimvar=Tkinter.StringVar()
        miscframe.delimentry=Tkinter.Entry(miscframe,textvariable=self.delimvar,width=8)
        miscframe.delimentry.grid(column=1,row=1)
        self.delimvar.set(",")
        
        # Starting column
        self.sctext = Tkinter.StringVar()
        miscframe.sclab = Tkinter.Label(miscframe,textvariable=self.sctext)
        miscframe.sclab.grid(row=2,column=0,sticky='w')
        self.sctext.set("Startcol GPA file")
        
        self.scvar = Tkinter.StringVar()
        miscframe.sc = Tkinter.Entry(miscframe,textvariable=self.scvar,width=8)
        miscframe.sc.grid(row=2,column=1)
        self.scvar.set("15")
        
        # Permutations
        self.permtext = Tkinter.StringVar()
        miscframe.permlab = Tkinter.Label(miscframe,textvariable=self.permtext)
        miscframe.permlab.grid(row=3,column=0,sticky='w')
        self.permtext.set("Permutations")
        
        self.permvar = Tkinter.StringVar()
        miscframe.perm = Tkinter.Entry(miscframe,textvariable=self.permvar,width=8)
        miscframe.perm.grid(row=3,column=1)
        self.permvar.set("0")
        
        # No timestamp
        self.notimevar = Tkinter.IntVar()
        miscframe.notime = Tkinter.Checkbutton(miscframe,text=u"No timestamp",onvalue=1,offvalue=0,variable=self.notimevar)
        miscframe.notime.grid(row=5,column=0,sticky='w')
        
        # Write tree
        self.writetreevar = Tkinter.IntVar()
        miscframe.writetree = Tkinter.Checkbutton(miscframe,text=u"Write tree",onvalue=1,offvalue=0,variable=self.writetreevar)
        miscframe.writetree.grid(row=6,column=0,sticky='w')
        
    def initialize_pvalueframe(self):
        """
        Initialize the filtration controlboard
        """
        board = self.nepart
        pframe = board.pframe
        
        # Add p-value checkboxes
        
        self.pVar = Tkinter.IntVar()
        self.pVar.set(1)
        self.pBVar = Tkinter.IntVar()
        self.pBHVar = Tkinter.IntVar()
        self.pPWVar = Tkinter.IntVar()
        self.pEPWVar = Tkinter.IntVar()
        self.pEPWVar.set(1)
        self.pPermVar = Tkinter.IntVar()
        
        pframe.pNcheck = Tkinter.Checkbutton(pframe,text=u"Naive (Fisher's)",onvalue=1,offvalue=0, variable=self.pVar)
        pframe.pBcheck = Tkinter.Checkbutton(pframe,text=u"Bonferroni",onvalue=1,offvalue=0, variable=self.pBVar)
        pframe.pBHcheck = Tkinter.Checkbutton(pframe,text=u"Benjamini-Hochberg",onvalue=1,offvalue=0, variable=self.pBHVar)
        pframe.pPWcheck = Tkinter.Checkbutton(pframe,text=u"Pairwise comparison (Best)",onvalue=1,offvalue=0, variable=self.pPWVar)
        pframe.pEPWcheck = Tkinter.Checkbutton(pframe,text=u"Pairwise comparison (Entire)",onvalue=1,offvalue=0, variable=self.pEPWVar)
        pframe.pPermcheck = Tkinter.Checkbutton(pframe,text=u"Empirical p-value (Permutation)", onvalue=1,offvalue=0, variable=self.pPermVar)
        
        pframe.pNcheck.grid(column=0,row=0,sticky='w')
        pframe.pBcheck.grid(column=0,row=1,sticky='w')
        pframe.pBHcheck.grid(column=0,row=2,sticky='w')
        pframe.pPWcheck.grid(column=0,row=3,sticky='w')
        pframe.pEPWcheck.grid(column=0,row=4,sticky='w')
        pframe.pPermcheck.grid(column=0,row=5,sticky='w')
        
        # Add p-value entry cells
        
        self.pNaive = Tkinter.StringVar()
        self.pBonf = Tkinter.StringVar()
        self.pBH = Tkinter.StringVar()
        self.pPW = Tkinter.StringVar()
        self.pEPW = Tkinter.StringVar()
        self.pPerm = Tkinter.StringVar()
        
        pframe.pNaiveEntry = Tkinter.Entry(pframe,textvariable=self.pNaive,width=8)
        pframe.pBonfEntry = Tkinter.Entry(pframe,textvariable=self.pBonf,width=8)
        pframe.pBHEntry = Tkinter.Entry(pframe,textvariable=self.pBH,width=8)
        pframe.pPWEntry = Tkinter.Entry(pframe,textvariable=self.pPW,width=8)
        pframe.pEPWEntry = Tkinter.Entry(pframe,textvariable=self.pEPW,width=8)
        pframe.pPermEntry = Tkinter.Entry(pframe,textvariable=self.pPerm,width=8)
        
        pframe.pNaiveEntry.grid(column=1,row=0,sticky='w')
        pframe.pBonfEntry.grid(column=1,row=1,sticky='w')
        pframe.pBHEntry.grid(column=1,row=2,sticky='w')
        pframe.pPWEntry.grid(column=1,row=3,sticky='w')
        pframe.pEPWEntry.grid(column=1,row=4,sticky='w')
        pframe.pPermEntry.grid(column=1,row=5,sticky='w')
        
        self.pNaive.set("0.05")
        self.pBonf.set("1.0")
        self.pBH.set("1.0")
        self.pPW.set("1.0")
        self.pEPW.set("0.05")
        self.pPerm.set("1.0")
        
        ######## EVENTS ########
        
    def AboutScoary(self):
        """
        Placeholder button. Planned short information about the method
        """
        topwin = Tkinter.Toplevel(self)
        button = Tkinter.Button(topwin,text=str("https://github.com/AdmiralenOla/Scoary"))
        button.pack()
        
    def BrowseButtonClickGPA(self):
        """
        Browse button for gene presence absence field
        """
        myfile = tkFileDialog.askopenfilename(filetypes=[('comma-separated values', '.csv'), ('all files','.*')])
        self.GPAentryVariable.set(myfile)
        
    def BrowseButtonClickTraits(self):
        """
        Browse button for traits field
        """
        myfile = tkFileDialog.askopenfilename(filetypes=[('comma-separated values', '.csv'), ('all files','.*')])
        self.TraitsentryVariable.set(myfile)
        
    def BrowseButtonClickTreeFile(self):
        """
        Browse button for tree field
        """
        myfile = tkFileDialog.askopenfilename(filetypes=[('newick tree files', '.nwk'), ('all files','.*')])
        self.TreeentryVariable.set(myfile)
        
    def BrowseButtonClickRestrict(self):
        """
        Browse button for isolate restriction field
        """
        myfile = tkFileDialog.askopenfilename(filetypes=[('comma-separated values','.csv'),('all files','.*')])
        self.RestrictVariable.set(myfile)            
        
    def BrowseButtonClickOutput(self):
        """
        Browse button for choosing output dir
        """
        mydir = tkFileDialog.askdirectory(mustexist=True)
        self.Outputdir.set(mydir)
    
    def HelpButton(self):
        """
        Placeholder button. Redirects to website
        """
        print("Visit https://github.com/AdmiralenOla/Scoary for help")
        
    def ClearAll(self):
        """
        Sets all variables to defaults
        """
        self.GPAentryVariable.set("Path to gene presence absence file")
        self.TraitsentryVariable.set("Path to traits/phenotype file")
        self.TreeentryVariable.set("(Optional) Path to custom tree file")
        self.RestrictVariable.set("(Optional) Path to file naming isolates to include")
        self.Outputdir.set("")
        self.maxhitsvar.set("")
        self.delimvar.set(",")
        self.scvar.set("15")
        self.notimevar.set(0)
        self.writetreevar.set(0)
        self.permvar.set("0")
        self.pVar.set(1)
        self.pBVar.set(0)
        self.pBHVar.set(0)
        self.pPWVar.set(0)
        self.pEPWVar.set(1)
        self.pNaive.set("0.05")
        self.pBonf.set("1.0")
        self.pBH.set("1.0")
        self.pPW.set("1.0")
        self.pEPW.set("0.05")
        self.pPerm.set("1.0")
        
    def TestExample(self):
        """
        Sets all variables corresponding to --test in the methods script
        """
        self.GPAentryVariable.set(str(os.path.join(resource_filename(__name__, 'exampledata'), 'Gene_presence_absence.csv')))
        self.TraitsentryVariable.set(str(os.path.join(resource_filename(__name__, 'exampledata'), 'Tetracycline_resistance.csv')))
        self.TreeentryVariable.set("")
        self.RestrictVariable.set("")
        self.Outputdir.set("./")
        self.maxhitsvar.set("")
        self.delimvar.set(",")
        self.scvar.set("15")
        self.notimevar.set(0)
        self.writetreevar.set(0)
        self.permvar.set("0")
        self.pVar.set(1)
        self.pBVar.set(0)
        self.pBHVar.set(0)
        self.pPWVar.set(0)
        self.pEPWVar.set(1)
        self.pNaive.set("0.05")
        self.pBonf.set("1.0")
        self.pBH.set("1.0")
        self.pPW.set("1.0")
        self.pEPW.set("0.05")
        self.pPerm.set("1.0")
        
        ######## RUNNING THE ANALYSIS ########
        
    def RunAnalysis(self):
        """
        Upon click "Run analysis" - Reads all the set parameters and calls run method
        """
        self.Scoary_parameters["GPA"] = self.GPAentryVariable.get()
        self.Scoary_parameters["Traits"] = self.TraitsentryVariable.get()
        self.Scoary_parameters["Tree"] = self.TreeentryVariable.get()
        self.Scoary_parameters["Outdir"] = self.Outputdir.get()
        self.Scoary_parameters["Restrict"] = self.RestrictVariable.get()
        self.Scoary_parameters["Writetree"] = self.writetreevar.get()
        self.Scoary_parameters["Delimiter"] = self.delimvar.get()
        self.Scoary_parameters["Startcol"] = self.scvar.get()
        self.Scoary_parameters["Maxhits"] = self.maxhitsvar.get()
        self.Scoary_parameters["Notime"] = self.notimevar.get()
        self.Scoary_parameters["Permutations"] = self.permvar.get()
        self.Scoary_parameters["Cutoffs"]["I"][0] = self.pVar.get()
        self.Scoary_parameters["Cutoffs"]["I"][1] = self.pNaive.get()
        self.Scoary_parameters["Cutoffs"]["B"][0] = self.pBVar.get()
        self.Scoary_parameters["Cutoffs"]["B"][1] = self.pBonf.get()
        self.Scoary_parameters["Cutoffs"]["BH"][0] = self.pBHVar.get()
        self.Scoary_parameters["Cutoffs"]["BH"][1] = self.pBH.get()
        self.Scoary_parameters["Cutoffs"]["PW"][0] = self.pPWVar.get()
        self.Scoary_parameters["Cutoffs"]["PW"][1] = self.pPW.get()
        self.Scoary_parameters["Cutoffs"]["EPW"][0] = self.pEPWVar.get()
        self.Scoary_parameters["Cutoffs"]["EPW"][1] = self.pEPW.get()
        self.Scoary_parameters["Cutoffs"]["P"][0] = self.pPermVar.get()
        self.Scoary_parameters["Cutoffs"]["P"][1] = self.pPerm.get()
        
        self.PrepareScoaryCMDline()
        
    def PrepareScoaryCMDline(self):
        """
        Prepares arguments to correspond with argparse namespace and runs.
        Listens to sys.stdout and updates statusbar while scoary runs
        """
        citation=False
        RunScoary = True
        correction = []
        p_value_cutoff = []
        for m in self.Scoary_parameters["Cutoffs"]:
            if self.Scoary_parameters["Cutoffs"][m][0] == 1:
                correction.append(m)
                try:
                    p_value_cutoff.append(float(self.Scoary_parameters["Cutoffs"][m][1]))
                except ValueError:
                    print("Please enter real numbers in the p value fields")
                    RunScoary = False
                    break
                    
        delimiter = self.Scoary_parameters["Delimiter"]
        genes = self.Scoary_parameters["GPA"] if self.Scoary_parameters["GPA"] not in ["","Path to gene presence absence file"] else None
        try:
            max_hits = int(self.Scoary_parameters["Maxhits"]) if self.Scoary_parameters["Maxhits"] not in [""] else None
        except ValueError:
            print("Please enter a real number (or nothing) in the max hits field")
            max_hits = None
            RunScoary = False
        try:    
            permutations = abs(int(self.Scoary_parameters["Permutations"])) if self.Scoary_parameters["Permutations"] not in [""] else 0
        except ValueError:
            print("Please enter a real number (or nothing) in the max hits field")
            permutations = 0
            RunScoary = False
        newicktree = self.Scoary_parameters["Tree"] if self.Scoary_parameters["Tree"] not in ["", "(Optional) Path to custom tree file"] else None
        no_time = True if self.Scoary_parameters["Notime"] == 1 else False
        outdir = self.Scoary_parameters["Outdir"] if self.Scoary_parameters["Outdir"] not in ["","(Optional) Output directory"] else "./"
        restrict_to = self.Scoary_parameters["Restrict"] if self.Scoary_parameters["Restrict"] not in ["", "(Optional) Path to file naming isolates to include"] else None
        start_col = self.Scoary_parameters["Startcol"]
        test = False
        threads = 1
        traits = self.Scoary_parameters["Traits"] if self.Scoary_parameters["Traits"] not in ["","Path to traits/phenotype file"] else None
        upgma_tree = True if self.Scoary_parameters["Writetree"] == 1 else False
        write_reduced=False
        collapse=False
               
        myargs = argparse.Namespace(citation=citation, correction=correction, p_value_cutoff=p_value_cutoff, delimiter=delimiter, genes=genes, max_hits=max_hits, newicktree=newicktree,no_time=no_time,restrict_to=restrict_to,
        outdir=outdir,permute=permutations,start_col=start_col,test=test,threads=threads,traits=traits,upgma_tree=upgma_tree,write_reduced=write_reduced,collapse=collapse)
        
        if RunScoary:
            try:
                
                sm.main(args=myargs, cutoffs=dict(list(zip(correction, p_value_cutoff))),statusbar=sys.stdout)

            except SystemExit as SE:
                # Set status bar color to red?
                if str(SE) == "0":
                    print("Analysis complete!")
                else:
                    print("Fatal error: %s" % str(SE))
        
        # Listen to stdout and update statusbar
        
    ######## MISC METHODS ########
        
    def initialize_logo(self):
        """
        Initialize logo
        """
        photo=Tkinter.PhotoImage(data=self.Photobase64())
        self.logocanvas.img = photo
        self.logocanvas.create_image(0,0,anchor='nw',image=photo)
        
    def citationtext(self):
        """
        Returns citation info
        """
        text = "SCOARY version %s \n\n" \
        "Please cite as: \n" \
        "Brynildsrud O. Scoary: \n" \
        "Microbial Pan-GWAS. \n" \
        "https://github.com/AdmiralenOla/Scoary" % scoary.__version__
        return text
        
    def Photobase64(self):
        """
        base64 encoding of the logo
        """
        photo= """
        R0lGODlh+gD6AOf/AAABAAEEAAMGAgcJBQoMCAwPCw8RDhIUERYYFhkbGRwdGx8gHiIkISYoJSstKjsrLkArL0gqJ00pI0QsJzEyMEQtLT8wMkYvL0IxLkUwM0svMDQ2M0IyNE0wMUQzMEMzNUoyMjY4NUQ1N0Y1\
        MkI3ODk7OEc3OUk5PEU7PD4/PU07OUU9QkE/Q04+QD9DRUFDQERETVBBSEtDSE1DQ0dFSEZHRVRDRlZDQUFJUUlLSD1PWVRKSjhRWkxOSzRVY1xMUEdSXj5VZFBST0BWYFNVUjRcbl9VWFZYVWJVUT9ebEtdaDZj\
        dVtdWkBjd2lcVyNsh2hcXSxrgV9hXjlsg2VnZDV0inFmZzF3kmpsaXZpZyl8lh5/nit9nm5wbSeEo3hzcjiCnXN1ch6JroBzcnZ4dYF2di+KqRaQuoB4cnh6dyeOs4d5dBiUuASaw3x+e4CBf4mBeyyXwACj0wCh\
        8YyFgIWHhJCEhBugygCm8ACo8QCp8h6lzhGo2ACr7h+i7QCr9IuNigCu4wCw3gqt8I6QjRyr2wCy9BCu8QC16ZqSjJKUkQC28Q2z6QC2+Bew8wC59AC67wC6/Be166GWip6XkQC895eZlgC9/iC18qSclSa49QjB\
        9qKdnJ2fnBDA+yi78aWgny2+9KKkoTLA9qWnpKynpamrqEfD9Kyuq7OtqjvJ+EvG97iysbK0sE7M9lvK9be5tl/N+bu9umDS9768wMS8tWzQ9sDCv8bBwHvT+3rW94XW+MbIxM/HwHza+4fY+snLyM7KzIjd95Hb\
        +M3PzJXf/J3e/NfSytjS0dPV0ZHl+aXg+drWx57j+tXX1Kjj/N3a0dnb2Kjn+bLo+7nn++Te2t7h3ePi2bvu+8Ps+uPl4ubm3Mbv/c3t/u7o5+jq5svy+uzr4tTx++3s4+zu6/Pt7Nj0/+7w7d7z//Hw5t32++T1\
        /PL08Pjz8ff27ej5//T38/D4/+z6+/P4++77/PX6/fj69/r8+ff9//n+///9+/z++//9////9v7//P///yH+C1Njb2FyeSBsb2dvACH5BAEKAP8ALAAAAAD6APoAAAj+AP0JHEiwoMGDCBMqXMiwocOHECNKnEix\
        osWLGDNq3Mixo8ePIEOKHEmypMmTKFOqXMmypcuXMGPKnEmzps2bOHPq3Mmzp8+fQIMKHUq0qNGjSJMqXcq0qdOnUKNKnUq1qtWrWLNq3cq1q9evYMOKHUu2rNmzaNOqXcu2rdu3cOPKnUu3rt27ePPq3cu3r9+/\
        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
        """
        return photo
        
class StdoutToLabel(Tkinter.Label):
    """
    The special widget that listens to sys.stdout and updates its
    label accordingly. This widget also owns the progressbar it is
    placed on and updates it every time sys.stdout is called.
    Use stop() to stop retrieving from sys.stdout.
    Credit: http://stackoverflow.com/questions/2914603/
    """

    def __init__(self, widget, progressbar=None, width='default'):
        Tkinter.Label.__init__(self)
        self.defstdout = sys.stdout
        self.widget = widget
        self.progressbar = progressbar
        
        if width == 'default':
            try:
                self.width = widget.cget('width')
            except:
                self.width = None
        else:
            self.width = width   

    def flush(self):
        """
        Frame sys.stdout's flush method.
        """
        self.defstdout.flush()

    def write(self, string, end=None):
        """
        Frame sys.stdout's write method. Puts input
        strings in the widget.
        """
        self.defstdout.write(string)
        if string == "\n":
            pass
        elif len(string) > self.width and self.width is not None and self.width > 0:
            self.widget.config(text="See log file for details")
        elif string is not None:

            last_char = string[-1]

            if last_char != '\n' and string.startswith('\r'):
                number = string.lstrip("\r").rstrip("%")
                number = number.split(".")[0]
                number = int(number)
                self.progressbar["value"] = number                
                
            elif last_char != '\n':
                self.widget.config(text=string)

            elif last_char == '\n' and string.startswith('\r'):
                number = string.lstrip("\r").rstrip("%")
                number = number.split(".")[0]
                number = int(number)
                self.progressbar["value"] = number
                
            else:
                self.widget.config(text=string[:-1])
                
        self.update()
        self.progressbar.update()

    def start(self):
        """
        Starts retrieving.
        """
        sys.stdout = self

    def stop(self):
        """
        Stops retrieving.
        """
        sys.stdout = self.defstdout

    @property
    def errors(self):
        return self.defstdout.errors

    @property
    def encoding(self):
        return self.defstdout.encoding
        
def main():
    root = ScoaryGUI(None)
    root.title("Scoary")
    root.geometry("800x400")
    root.resizable(0,0)
    root.mainloop()
    
if __name__ == "__main__":
    pass
