#!/usr/bin/env python
# coding: utf8

from Tkinter import *
from ttk import *
#import Tix
import calculate_energy as calc
import add_structure as adds
import pwm as pwmpy
import occupancy as occupancypy
from tkFileDialog import askopenfilename
import os
#from Tkconstants import *
import tkFont
import tkMessageBox


class Dnapp(Tk):
    def __init__(self,parent):
        Tk.__init__(self)
        self.parent = parent
        self.initialize()

    def HelloWorld(self):
        print ""
        
    def initialize(self):
        self.grid()
        self.rowconfigure(0, pad=30)
        self.fastas=StringVar()
        self.cnfs=StringVar()
        self.fastas.trace("w",self.enable_launch)
        self.cnfs.trace("w",self.enable_launch)

        self.fasta=""
        self.cnf=""
        
        # fonts
        tf=tkFont.Font(family="Sans",size=12,weight="bold")
        bf=tkFont.Font(family="Sans",size=10,weight="bold")
        mf=tkFont.Font(family="Sans",size=10)
        
        s=Style()
        s.theme_use('classic')

        # ---------------------------------------

        exf=Frame(self,borderwidth="3m",relief="ridge")
        exf.grid(column=0,row=0,columnspan=3,sticky="N")
        
        lcalc=Label(exf,text=u"Main program: \ncompute energy along sequence        ",font=tf)
        lcalc.grid(column=0,row=0,columnspan=2)
        raj=Label(exf,text=u"",font=("Sans","6"))
        raj.grid(column=0,row=1,columnspan=3,sticky='W')
        w=Separator(exf,orient=HORIZONTAL)
        w.grid(column=0,row=1,columnspan=3,sticky='EW')
        #lsep=Label(exf)
        #lsep.grid(column=0,row=1,columnspan=2)
                
        cnfl=Label(exf,text=u"Input configuration file",font=bf)
        cnfl.grid(column=0,row=2,columnspan=2,sticky='W')
        opcnf = Button(exf,text=u"   Open config file   ",command=self.open_cnf)
        opcnf.grid(column=0,row=3,sticky='W')
        cnfh=Label(exf,text="The configuration file describes\nthe protein model with various\noptions: choose an existing file,\nor use a simplified interface to\ncreate a new configuration file")
        cnfh.grid(column=2,row=2,rowspan=3,sticky='W')
        lcnf = Label(exf,textvariable=self.cnfs)
        lcnf.grid(column=1,row=3,sticky='EW')
        self.cfg=Button(exf,text=u"or create config file",command=self.edit_config)
        self.cfg.grid(column=0,row=4,sticky="W")
        #raj=Label(exf,text=u"",font=("Sans","12"))
        #raj.grid(column=0,row=8,columnspan=2,sticky='W')

        fasl=Label(exf,text=u"Optional: input sequence file (fasta)",font=bf)
        fasl.grid(column=0,row=5,columnspan=2,sticky='W')
        opfasta = Button(exf,text=u"Open FASTA",command=self.open_fasta)
        opfasta.grid(column=0,row=6,sticky='W')
        lfasta = Label(exf,textvariable=self.fastas)
        lfasta.grid(column=1,row=6,sticky='EW')
        fash=Label(exf,text="Compute an energy profile\nalong the given sequence.\nOtherwise a PWM file is generated.")
        fash.grid(column=2,row=5,rowspan=2,sticky='W')
        w=Separator(exf,orient=HORIZONTAL)
        w.grid(column=0,row=7,columnspan=3,sticky='EW')
        
        # w=Separator(exf,orient=HORIZONTAL)
        # w.grid(column=0,row=8,columnspan=3,sticky='EW')
        self.run=Button(exf,text=u"Calculate Energy",command=self.launch_calc,state="disabled")
        self.run.grid(column=0,row=9,columnspan=1)
        runh=Label(exf,text="Starts computation of energy profile and/or PWM file.\nYou will get a message once the computation is complete.")
        runh.grid(column=1,row=9,columnspan=2,sticky='W')
        #w=Separator(exf,orient=HORIZONTAL)
        #w.grid(column=0,row=10,columnspan=2,sticky='EW')

        exf2=Frame(self,borderwidth="3m",relief="ridge",width=500,height=180)
        exf2.grid(column=0,row=10,columnspan=2,sticky="N")
        lcalc=Label(exf2,text=u"Subprogram: add protein to database",font=tf)
        lcalc.grid(column=0,row=11,columnspan=2)
        self.add=Button(exf2,text=u"   Add protein   ",command=self.new_prot)
        self.add.grid(column=0,row=12,columnspan=1)
        #w.grid(column=0,row=10,columnspan=13,sticky='EW')
        addh=Label(exf2,text="  For a new protein, you must launch this subprogram     \n  to add it to the database, before any energy profile  \n can be computed by the main program")
        addh.grid(column=1,row=12,columnspan=1,rowspan=1,sticky='W')
        
        exf3=Frame(self,borderwidth="3m",relief="ridge",width=500,height=150)
        exf3.grid(column=0,row=13,columnspan=2,sticky="N")
        lcalc=Label(exf3,text=u"Subprogram: PWM computations   ",font=tf)
        lcalc.grid(column=0,row=14,columnspan=2)
        self.hel=Button(exf3,text=u"PWM computations",command=self.pwm)
        self.hel.grid(column=0,row=15,columnspan=1)
        helph=Label(exf3,text="Manipulate mono/dinucleotide position-weight-matrix \nin combination with ThreaDNA")
        helph.grid(column=1,row=15,columnspan=1,rowspan=1,sticky='W')
        
        exf4=Frame(self,borderwidth="3m",relief="ridge",width=500,height=150)
        exf4.grid(column=0,row=16,columnspan=2,sticky="N")
        lcalc=Label(exf4,text=u"Helper program: occupancy profile      ",font=tf)
        lcalc.grid(column=0,row=17,columnspan=2)
        self.occ=Button(exf4,text=u" Occupancy ",command=self.occupancy)
        self.occ.grid(column=0,row=18,columnspan=1)
        occh=Label(exf4,text="Translate free energy profile, such as provided by ThreaDNA,\ninto occupancy/coverage profile of a protein")
        occh.grid(column=1,row=18,columnspan=1,rowspan=1,sticky='W')



    def new_prot(self):

        self.cfgpan=Toplevel()
        self.cfgpan.grab_set()
        self.cfgpan.grid()

        #self.grid()
        self.rowconfigure(0, pad=30)
        self.prota=StringVar()
        self.strucs=StringVar()
        self.ida=StringVar()
        self.pdbs=StringVar()
        self.ref=IntVar()
        self.ida.trace("w",self.enable_add)
        self.prota.trace("w",self.enable_add)
        self.strucs.trace("w",self.enable_add)
        self.ref.trace("w",self.enable_add)

        self.struc=""
        self.pdb=""
        
        # fonts
        tf=tkFont.Font(family="Sans",size=12,weight="bold")
        bf=tkFont.Font(family="Sans",size=10,weight="bold")
        mf=tkFont.Font(family="Sans",size=10)
        
        s=Style()
        s.theme_use('classic')
        
        #addf=Frame(self,borderwidth="3m",relief="ridge")
        #addf.grid(column=0,row=0,columnspan=2,sticky='N')

        addl=Label(self.cfgpan,text=u"Add a new protein structure to the local database",font=tf)
        addl.grid(column=0,row=0,columnspan=3)
        w=Separator(self.cfgpan,orient=HORIZONTAL)
        w.grid(column=0,row=1,columnspan=3,sticky='EW')

        lprot=Label(self.cfgpan,text=u"Protein name",font=bf)
        lprot.grid(column=0,row=2,sticky='EW')
        prot = Entry(self.cfgpan,textvariable=self.prota,font=mf)
        prot.grid(column=1,row=2,sticky='EW')
        lproth=Label(self.cfgpan,text="Name of the protein in the database (mind the case\nfor different structures of the same protein)")
        lproth.grid(column=2,row=2,sticky='EW')
        
        w=Separator(self.cfgpan,orient=HORIZONTAL)
        w.grid(column=0,row=3,columnspan=3,sticky='EW')

        lprot=Label(self.cfgpan,text=u"DNA structure in complex",font=bf)
        lprot.grid(column=0,row=4,columnspan=2,sticky='EW')
        lid=Label(self.cfgpan,text="Enter NDB ID",font=mf,justify=RIGHT)
        lid.grid(column=0,row=5,sticky='NW')
        eid=Entry(self.cfgpan,textvariable=self.ida)
        eid.grid(column=1,row=5,sticky='EW')
        NDBh=Label(self.cfgpan,text="Provide ID of the structure in the NDB database. \nOtherwise, a base-pair (step) coordinate file,\nfrom software 3DNA (.out) or Curves+ (.lis).\nCAUTION: check that no basepair is missing!")
        NDBh.grid(column=2,row=4,rowspan=3,sticky='EW')
        opstruc = Button(self.cfgpan,text=u"or open coord file",command=self.open_struc)
        opstruc.grid(column=0,row=6,sticky='NW')
        
        lstruc = Label(self.cfgpan,textvariable=self.strucs)
        lstruc.grid(column=1,row=6,sticky='EW')
        w=Separator(self.cfgpan,orient=HORIZONTAL)
        w.grid(column=0,row=7,columnspan=3,sticky='EW')

        lprot=Label(self.cfgpan,text="Recommended optional parameters: ",font=bf)
        lprot.grid(column=0,row=8,columnspan=2,sticky='EW')
        oppdb = Button(self.cfgpan,text=u"Additional PDB file",command=self.open_pdb)
        oppdb.grid(column=0,row=9,sticky='NW')
        lpdb = Label(self.cfgpan,textvariable=self.pdbs)
        lpdb.grid(column=1,row=9,sticky='EW')
        PDBh=Label(self.cfgpan,text="Gives ThreaDNA important information on structure")
        PDBh.grid(column=2,row=9,rowspan=1,sticky='EW')
        lref=Label(self.cfgpan,text=u"Reference position")
        lref.grid(column=0,row=10,sticky='NW')
        eref=Entry(self.cfgpan,textvariable=self.ref)
        eref.grid(column=1,row=10,sticky='EW')
        refh=Label(self.cfgpan,text="Index of the (reference) basepair facing a reference\nposition on the protein. Default (0): central basepair.")
        refh.grid(column=2,row=10,rowspan=1,sticky='EW')
        w=Separator(self.cfgpan,orient=HORIZONTAL)
        w.grid(column=0,row=11,columnspan=3,sticky='EW')

        self.addb=Button(self.cfgpan,text=u"Add Structure",command=self.add_struct,state="disabled")
        self.addb.grid(column=0,row=12,columnspan=2)
        addbh=Label(self.cfgpan,text="Add structure to ThreaDNA database.")
        addbh.grid(column=2,row=12,rowspan=1,sticky='EW')


    """
    def helper(self):
        self.cfgpan=Toplevel()
        self.cfgpan.grab_set()
        self.cfgpan.grid()

        #self.grid()
        #self.ida=StringVar()
        #self.ida.trace("w",self.enable_add)
        
        self.rowconfigure(0, pad=30)
        self.fastas=StringVar()
        self.fastas.trace("w",self.enable_launch_helper)
        self.nuc=StringVar()
        self.nuc.trace("w",self.enable_launch_helper)
        self.ref=IntVar()
        self.ref.trace("w",self.enable_launch_helper)
        self.pos=IntVar()
        self.pos.trace("w",self.enable_launch_helper)
        self.leng=IntVar()
        self.leng.trace("w",self.enable_launch_helper)

        self.fasta=""
        #self.pos=0
        #self.leng=0
        #self.ref=0
        

         # fonts
        tf=tkFont.Font(family="Sans",size=12,weight="bold")
        bf=tkFont.Font(family="Sans",size=10,weight="bold")
        mf=tkFont.Font(family="Sans",size=10)
        
        s=Style()
        s.theme_use('classic')

        #penf=Frame(self,borderwidth="3m",relief="ridge")
        #penf.grid(column=5,row=0,columnspan=2,sticky="N")

        penlcalc=Label(self.cfgpan,text=u"Helper program:\nsequence patterns\nin a protein-DNA complex",font=tf)
        penlcalc.grid(column=5,row=0,columnspan=2)
        penraj=Label(self.cfgpan,text=u"",font=("Sans","3"))
        penraj.grid(column=5,row=1,columnspan=2,sticky='W')
        w=Separator(self.cfgpan,orient=HORIZONTAL)
        w.grid(column=5,row=1,columnspan=2,sticky='EW')

        penopfasta = Button(self.cfgpan,text=u"Open FASTA: ",command=self.open_fasta)
        penopfasta.grid(column=5,row=2,sticky='N')
        penlfasta = Label(self.cfgpan,textvariable=self.fastas)
        penlfasta.grid(column=6,row=2,sticky='N')

        penlnuc=Label(self.cfgpan,text=u"Nucleotide: ",font=mf)
        penlnuc.grid(column=5,row=3,sticky='N')
        penenuc=Entry(self.cfgpan,textvariable=self.nuc)
        penenuc.grid(column=6,row=3,sticky='N')
        
        penlpos=Label(self.cfgpan,text=u"Position in the complex: ",font=mf)
        penlpos.grid(column=5,row=4,sticky='N')
        penepos=Entry(self.cfgpan,textvariable=self.pos)
        penepos.grid(column=6,row=4,sticky='N')
        
        penllen=Label(self.cfgpan,text=u"Length of the complex: ",font=mf)
        penllen.grid(column=5,row=5,sticky='N')
        penelen=Entry(self.cfgpan,textvariable=self.leng)
        penelen.grid(column=6,row=5,sticky='N')

        penlref=Label(self.cfgpan,text=u"Optional:\nReference nucleotide\nof the complex ",font=mf)
        penlref.grid(column=5,row=6,rowspan=3,sticky='N')
        peneref=Entry(self.cfgpan,textvariable=self.ref)
        peneref.grid(column=6,row=8,sticky='N')
        self.helperrun=Button(self.cfgpan,text=u"Compute\npattern",command=self.run_helper,state="disabled")
        self.helperrun.grid(column=5,row=9,columnspan=2,rowspan=2)
    """

    def pwm(self):
        self.cfgpan=Toplevel()
        self.cfgpan.grab_set()
        self.cfgpan.grid()

        #self.grid()
        #self.ida=StringVar()
        #self.ida.trace("w",self.enable_add)
        
        self.rowconfigure(0, pad=30)
        self.fastas=StringVar()
        self.epwms=StringVar()
        self.ppwms=StringVar()
        #self.pwms.trace("w",self.enable_launch_pwm)
        self.energs=StringVar()
        #self.energ.trace("w",self.enable_launch_helper)

        self.energs.set("1.")

        
        self.fasta=""
        self.epwm=""
        self.ppwm=""
        #self.pos=0
        self.energ="1."
        #self.ref=0
        

         # fonts
        tf=tkFont.Font(family="Sans",size=12,weight="bold")
        bf=tkFont.Font(family="Sans",size=10,weight="bold")
        mf=tkFont.Font(family="Sans",size=10)
        
        s=Style()
        s.theme_use('classic')
        
        #penf=Frame(self,borderwidth="3m",relief="ridge")
        #penf.grid(column=5,row=0,columnspan=2,sticky="N")

        penlcalc=Label(self.cfgpan,text=u"Subprogram:\ncomputations on position-weight-matrices\nin combination with ThreaDNA",font=tf)
        penlcalc.grid(column=5,row=0,columnspan=2)
        penraj=Label(self.cfgpan,text=u"",font=("Sans","3"))
        penraj.grid(column=5,row=1,columnspan=2,sticky='W')
        w=Separator(self.cfgpan,orient=HORIZONTAL)
        w.grid(column=5,row=1,columnspan=2,sticky='EW')

        expl=Label(self.cfgpan,text=u"You can provide a sequence file and/or a PWM file (JASPAR format).\nThe computation depends on the input type (see readme file)",font=mf)
        expl.grid(column=5,row=2,columnspan=2)
        penopfasta = Button(self.cfgpan,text=u"Sequence file: ",command=self.open_fasta)
        penopfasta.grid(column=5,row=3,sticky='N')
        penlfasta = Label(self.cfgpan,textvariable=self.fastas)
        penlfasta.grid(column=6,row=3,sticky='N')

        penopfasta = Button(self.cfgpan,text=u"Energy PWM file: ",command=self.open_epwm)
        penopfasta.grid(column=5,row=4,sticky='N')
        penlfasta = Label(self.cfgpan,textvariable=self.epwms)
        penlfasta.grid(column=6,row=4,sticky='N')

        penopfasta = Button(self.cfgpan,text=u"Proba PWM file: ",command=self.open_ppwm)
        penopfasta.grid(column=5,row=5,sticky='N')
        penlfasta = Label(self.cfgpan,textvariable=self.ppwms)
        penlfasta.grid(column=6,row=5,sticky='N')

        
        # penlnuc=Label(self.cfgpan,text=u"Nucleotide: ",font=mf)
        # penlnuc.grid(column=5,row=3,sticky='N')
        # penenuc=Entry(self.cfgpan,textvariable=self.nuc)
        # penenuc.grid(column=6,row=3,sticky='N')
        
        # penlpos=Label(self.cfgpan,text=u"Position in the complex: ",font=mf)
        # penlpos.grid(column=5,row=4,sticky='N')
        # penepos=Entry(self.cfgpan,textvariable=self.pos)
        # penepos.grid(column=6,row=4,sticky='N')
        
        # penllen=Label(self.cfgpan,text=u"Length of the complex: ",font=mf)
        # penllen.grid(column=5,row=5,sticky='N')
        # penelen=Entry(self.cfgpan,textvariable=self.leng)
        # penelen.grid(column=6,row=5,sticky='N')

        penlref=Label(self.cfgpan,text=u"Optional:\nEnergy rescaling factor ",font=mf)
        penlref.grid(column=5,row=6,rowspan=2,sticky='N')
        peneref=Entry(self.cfgpan,textvariable=self.energs)
        peneref.grid(column=6,row=6,sticky='N')

        w=Separator(self.cfgpan,orient=HORIZONTAL)
        w.grid(column=5,row=8,columnspan=2,sticky='EW')

        
        self.pwmrun=Button(self.cfgpan,text=u"Run computation",command=self.run_pwm)
        self.pwmrun.grid(column=5,row=9,columnspan=1,rowspan=2)

        self.pwmplot=Button(self.cfgpan,text=u"Plot PWM",command=self.plot_pwm)
        self.pwmplot.grid(column=6,row=9,columnspan=1,rowspan=2)



        

    def occupancy(self):
        self.cfgpan=Toplevel()
        self.cfgpan.grab_set()
        self.cfgpan.grid()

        #self.grid()
        #self.ida=StringVar()
        #self.ida.trace("w",self.enable_add)
        
        self.rowconfigure(0, pad=30)
        self.infs=StringVar()
        self.infs.trace("w",self.enable_launch_occ)
        self.rowconfigure(0, pad=30)
        self.efftemps=StringVar()
        self.coverages=StringVar()
        
        self.efftemps.set("1.0")


        self.inf=""
        self.efftemp="1.0"
        self.coverage=""
        #self.leng=0
        #self.ref=0

         # fonts
        tf=tkFont.Font(family="Sans",size=12,weight="bold")
        bf=tkFont.Font(family="Sans",size=10,weight="bold")
        mf=tkFont.Font(family="Sans",size=10)
        
        s=Style()
        s.theme_use('classic')

        #penf=Frame(self,borderwidth="3m",relief="ridge")
        #penf.grid(column=5,row=0,columnspan=2,sticky="N")

        penlcalc=Label(self.cfgpan,text=u"Helper program: from energy to occupancy profile",font=tf)
        penlcalc.grid(column=0,row=0,columnspan=3)
        penraj=Label(self.cfgpan,text=u"",font=("Sans","3"))
        penraj.grid(column=0,row=1,columnspan=3,sticky='W')
        w=Separator(self.cfgpan,orient=HORIZONTAL)
        w.grid(column=0,row=1,columnspan=3,sticky='EW')

        occinput = Button(self.cfgpan,text=u"Open input file\n(typically ThreaDNA output file): ",command=self.open_profile)
        occinput.grid(column=0,row=2,sticky='N')
        occinp = Label(self.cfgpan,textvariable=self.infs)
        occinp.grid(column=1,row=2,sticky='N')
        inph=Label(self.cfgpan,text="Computation of occupancy profile, i.e.\npositioning distribution of the protein")
        inph.grid(column=2,row=2,rowspan=1,sticky='EW')

        occtemp=Label(self.cfgpan,text=u"OPTIONAL: effective temperature: ",font=mf)
        occtemp.grid(column=0,row=4,sticky='N')
        occtempq=Entry(self.cfgpan,textvariable=self.efftemps)
        occtempq.grid(column=1,row=4,sticky='N')
        effth=Label(self.cfgpan,text="Effective temperature for Boltzmann\ninversion, in 300K reduced unit")
        effth.grid(column=2,row=4,rowspan=1,sticky='EW')

        occtemp=Label(self.cfgpan,text=u"OPTIONAL: coverage size: ",font=mf)
        occtemp.grid(column=0,row=5,sticky='N')
        occtempq=Entry(self.cfgpan,textvariable=self.coverages)
        occtempq.grid(column=1,row=5,sticky='N')
        covh=Label(self.cfgpan,text="If provided (e.g. 147 for nucleosome),\nan additional coverage profile is\ncomputed (probability that a given\nbp is covered by protein)")
        covh.grid(column=2,row=5,rowspan=1,sticky='EW')
        
        self.occrun=Button(self.cfgpan,text=u"Compute occupancy",command=self.run_occupancy,state="disabled")
        self.occrun.grid(column=0,row=6,columnspan=3,rowspan=1)
      
        


    def edit_config(self):
        self.cfgpan=Toplevel()
        self.cfgpan.grab_set()
        self.cfgpan.grid()
        self.l=list_struct(calc.loc+"structures")

        self.namev=StringVar()
        self.parv=StringVar()
        self.protv=StringVar()
        self.structv=[]
        self.ifsym=IntVar()
        self.iftmp=IntVar()
        self.tfv=DoubleVar()
        self.tempv=IntVar()
        self.patv=StringVar()

        self.tfv.set(1)
        self.parv.set("NP")
        self.tempv.set(300)
        self.namev.set("protein.in")
        
        self.namev.trace("w",self.enable_create)
        self.parv.trace("w",self.enable_create)
        self.tfv.trace("w",self.enable_create)
        self.tempv.trace("w",self.enable_create)


        # fonts
        tf=tkFont.Font(family="Sans",size=12,weight="bold")
        bf=tkFont.Font(family="Sans",size=10,weight="bold")
        mf=tkFont.Font(family="Sans",size=10)
        
        s=Style()
        s.theme_use('classic')


        self.protv.trace("w",self.update_prot)
        self.create=Button(self.cfgpan,text="Create configuration file",command=self.create_cnf,state="disabled")
        self.create.grid(column=0,row=11,columnspan=3)
        namel=Label(self.cfgpan,text=u"Config file name: ",font=bf)
        namel.grid(column=0,row=0,sticky='EW')
        name=Entry(self.cfgpan,textvariable=self.namev)
        name.grid(column=1,row=0,sticky='EW')
        cfnh=Label(self.cfgpan,text="Name of configuration file\n(created in working directory)")
        cfnh.grid(column=2,row=0,rowspan=1,sticky='EW')

        protl=Label(self.cfgpan,text=u"Protein name: ",font=bf)
        protl.grid(column=0,row=1,sticky='EW')
        #prot=Tix.ComboBox(self.cfgpan,variable=self.protv)
        prot=Frame(self.cfgpan)
        prot.grid(column=1,row=1)
        for x in self.l.keys():
            p=Radiobutton(prot,text=x,value=x,variable=self.protv)
            p.pack()
        #for p in self.l.keys():
        #    prot.insert(END,p)
        protnh=Label(self.cfgpan,text="Select the protein from\nall available in the database")
        protnh.grid(column=2,row=1,rowspan=1,sticky='EW')

        strucl=Label(self.cfgpan,text=u"Selected\nstructures:",font=bf)
        strucl.grid(column=0,row=2,sticky='EW')
        self.struc=Frame(self.cfgpan)
        self.struc.grid(column=1,row=2)
        protnh=Label(self.cfgpan,text="Select the structure(s) used in\nthe computation. With several\nstructures, the profiles will be\n computed separately, and then\ncombined into a single profile. ")
        protnh.grid(column=2,row=2,rowspan=1,sticky='EW')
        patl=Label(self.cfgpan,text="Pattern:",font=bf)
        patl.grid(column=0,row=3,sticky='EW')
        pat=Entry(self.cfgpan,textvariable=self.patv)
        pat.grid(column=1,row=3,sticky='EW')
        path=Label(self.cfgpan,text="Size of protein (in DNA bases)\nto include in the computation.\nDefault (empty): full structure.")
        path.grid(column=2,row=3,rowspan=1,sticky='EW')
        

        w=Separator(self.cfgpan,orient=HORIZONTAL)
        w.grid(column=0,row=4,columnspan=3,sticky='EW')
        opt=Label(self.cfgpan,text=u"Advanced parameters",font=bf)
        opt.grid(column=0,row=5,columnspan=3)
        sym=Checkbutton(self.cfgpan,text=u"Include symmetrized\nstructures",variable=self.ifsym)
        sym.grid(column=0,row=6,sticky='EW')
        tmp=Checkbutton(self.cfgpan,text=u"Keep tmp files",variable=self.iftmp)
        tmp.grid(column=1,row=6,sticky='EW')
        path.grid(column=2,row=3,rowspan=1,sticky='EW')
        tfl=Label(self.cfgpan,text=u"Temperature factor",font=mf)
        tfl.grid(column=0,row=7,sticky='EW')
        tf=Entry(self.cfgpan,textvariable=self.tfv)
        tf.grid(column=1,row=7,sticky='EW')
        TFh=Label(self.cfgpan,text="Sets the Boltzmann factor for\nstructure combination (in 300K unit)")
        TFh.grid(column=2,row=7,rowspan=1,sticky='EW')
        templ=Label(self.cfgpan,text=u"Temperature",font=mf)
        templ.grid(column=0,row=8,sticky='EW')
        temp=Entry(self.cfgpan,textvariable=self.tempv)
        temp.grid(column=1,row=8,sticky='EW')
        temph=Label(self.cfgpan,text="Physical temp for DNA flexibility")
        temph.grid(column=2,row=8,rowspan=1,sticky='EW')
        parl=Label(self.cfgpan,text=u"DNA parameter",font=mf)
        parl.grid(column=0,row=9,sticky='EW')
        par=Frame(self.cfgpan)
        par.grid(column=1,row=9)
        for x in ["NP","ABC_s","ABC_i"]:
            p=Radiobutton(par,text=x,value=x,variable=self.parv)
            p.pack()
        parh=Label(self.cfgpan,text="Parameter set for DNA flexibility:\nNP: exp, step, dinucleotide\nABC_s: simul, step, tetranucleotide\nABC_i: simul, intra-bp, trinucleotide")
        parh.grid(column=2,row=9,rowspan=1,sticky='EW')
        w=Separator(self.cfgpan,orient=HORIZONTAL)
        w.grid(column=0,row=10,columnspan=3,sticky='EW')

        #Sequence Energy:No

    def enable_add(self,*args):
        if self.ida.get() and self.strucs.get():
            self.strucs.set("")
            self.struc=""
        if self.ida.get() or self.struc and self.prota.get():
            try:
                self.ref.get()
                self.addb.config(state='normal')
            except:
                self.addb.config(state='disabled')
        else:
            self.addb.config(state='disabled')

    def enable_launch(self,*args):
        if self.cnf:
            self.run.config(state='normal')

    def enable_launch_helper(self,*args):
        if self.fasta and self.nuc and self.pos and self.leng:
            self.helperrun.config(state='normal')

    def enable_launch_occ(self,*args):
        if self.inf:
            self.occrun.config(state='normal')

    def enable_create(self,*args):
        try:
            self.tfv.get()
            self.tempv.get()
        except:
            self.create.config(state='disabled')
            return
        if self.namev.get() and self.protv.get() and self.parv.get() and "".join([x.get() for x in self.structv]):
            self.create.config(state='normal')
        else:
            self.create.config(state='disabled')

    def launch_calc(self):
        seq=self.fasta
        if seq=="":
            seq="None"
        name, namemat=calc.main(params=unicode(self.cnf),sequence=seq,output="None")
        if seq=="None":
            tkMessageBox.showinfo("Execution completed","Computation completed.\nPWM file in %s"%namemat)
        else:
            tkMessageBox.showinfo("Execution completed","Computation completed.\nEnergy profile in %s\nPWM file in %s"%(name,namemat))

    def add_struct(self):
        ref=self.ref.get()
        if ref==0:
            ref=None
        if self.struc:
            s=self.struc
        else:
            s=self.ida.get()
        try:
            mesg=adds.main(unicode(self.prota.get()),unicode(s),unicode(self.pdb),ref)
            tkMessageBox.showinfo("Execution completed",mesg)
        except:
            tkMessageBox.showinfo("Problem", "Error in the execution. A possible cause, if you provided the NDB ID, is that the coordinate file is absent from the database. We recommend you to provide the PDB file to the Web3DNA webserver, and provide the output .out file to ThreaDNA.")

    # def run_helper(self):
    #     if self.ref.get()==0:
    #         ref=self.leng.get()/2
    #     else:
    #         ref=self.ref.get()
    #     os.system("python seqmotifs.py -r %d %s %s %d %d"%(ref,unicode(self.fasta),unicode(self.nuc.get()),self.pos.get(),self.leng.get()))
    #     tkMessageBox.showinfo("Execution successful","Helper program execution completed")

    def run_pwm(self):
        fasta=self.fasta
        if self.fasta=="" or self.fasta==():
            fasta="None"
        epwm=self.epwm
        if self.epwm=="" or self.epwm==():
            epwm="None"
        ppwm=self.ppwm
        if self.ppwm=="" or self.ppwm==():
            ppwm="None"
        outf=pwmpy.main(fasta, epwm, ppwm, unicode(self.energs.get()), output="None")
        tkMessageBox.showinfo("Execution successful","Computation complete. Output in file\n%s"%outf)

    def plot_pwm(self):
        try:
            import plot_dinuc_pwm as plotpwmpy
            if self.ppwm=="":
                if self.epwm=="":
                    tkMessageBox.showinfo("Problem","Provide a probability PWM file to plot")
                    return 1
                else:
                    pwm=unicode(self.epwm)
            else:
                pwm=unicode(self.ppwm)
            outf=plotpwmpy.plot_pwm(pwm, outfile=None)
            tkMessageBox.showinfo("Execution successful","PWM %s was plotted to file\n%s"%(pwm,outf))
        except:
            tkMessageBox.showinfo("Problem","The plot could not be completed. One possible reason is that you don't have the MatPlotLib graphical library, that is generally not installed by default. If you don't want to install the library, you can plot traditional (mononucleotides) PWMs in JASPAR format on the STAMP website (http://www.benoslab.pitt.edu/stamp)")

        
    def run_occupancy(self):
        cov=self.coverages.get()
        if cov=="":
            cov=None
        outf=occupancypy.main(unicode(self.inf), output=None, efftemp=self.efftemps.get(), coverage=cov)
        tkMessageBox.showinfo("Execution successful","Execution completed. The result file is\n%s"%outf)


    def open_struc(self):
        self.struc=askopenfilename()
        self.strucs.set(self.struc.split("/")[-1])

    def open_pdb(self):
        self.pdb=askopenfilename()
        self.pdbs.set(self.pdb.split("/")[-1])

    def open_fasta(self):
        self.fasta=askopenfilename()
        self.fastas.set(self.fasta.split("/")[-1])

    def open_epwm(self):
        self.epwm=askopenfilename()
        self.epwms.set(self.epwm.split("/")[-1])

    def open_ppwm(self):
        self.ppwm=askopenfilename()
        self.ppwms.set(self.ppwm.split("/")[-1])
        
    def open_profile(self):
        self.inf=askopenfilename()
        self.infs.set(self.inf.split("/")[-1])

    def open_cnf(self):
        self.cnf=askopenfilename()
        self.cnfs.set(self.cnf.split("/")[-1])

    def create_cnf(self):
        f=open(os.getcwd()+"/"+self.namev.get(),"w")
        f.write("Protein:"+self.protv.get()+"\n")
        f.write("Structures:"+",".join(filter(None,[x.get() for x in self.structv]))+"\n")
        f.write("DNA parameter:"+self.parv.get()+"\n")
        if self.ifsym.get():
            f.write("Symmetrization:Yes\n")
        else:
            f.write("Symmetrization:No\n")
        if self.iftmp.get():
            f.write("Keep tmp:Yes\n")
        else:
            f.write("Keep tmp:No\n")
        f.write("Temperature factor:"+str(self.tfv.get())+"\n")
        if self.tempv.get()!=300:
            f.write("Temperature:"+str(self.tempv.get())+"\n")
        if self.patv.get():
            f.write("Pattern:"+self.patv.get()+"\n")
        f.close()
        self.cnf=self.namev.get()
        self.cnfs.set(self.namev.get())
        tkMessageBox.showinfo("Configuration file created","The file is located at %s"%(os.getcwd()+"/"+self.namev.get()))
        self.cfgpan.destroy()

    def update_prot(self,*args):
        for old in self.struc.winfo_children(): #clearing old structures
            old.destroy()
        self.structv=[]
        for p in self.l[self.protv.get()]:
            var=StringVar()
            s=Checkbutton(self.struc,text=p, variable=var,onvalue=p,offvalue="",command=self.enable_create)
            s.pack()
            self.structv.append(var)

def list_struct(path):
    struc={}
    for d,ds,fs in os.walk(path):
        if d==path:
            prot=ds
        elif d[len(path)+1:] in prot:
            struc[d[len(path)+1:]]=ds
    return struc

if __name__ == "__main__":
    app = Dnapp(None)
    app.title('ThreaDNA')
    app.mainloop()
