import argparse
import json
from math import ceil

class pDefaults:
    def __init__(self, id, N):
        self.rank        = id
        self.nProcessors = N

        self.tstart     = 0.0
        self.tend       = 0.0

        self.winsize    = None
        self.t1         = 0.0
        self.t2         = 0.0

        self.tol        = float(0.01)
        self.tolcheck   = 'L2'

        self.dtCoarse   = None
        self.dtFine     = None

        self.nCoarse    = 1
        self.nFine      = 10
        self.dbg        = 0
        self.outfile    = None
        self.nSave      = 10

        self.profile    = None

        self.faultfile  = None
        self.fault      = None

        self.distfile  = None
        self.dist      = None

    def parse_arguments(self):
        parser = argparse.ArgumentParser(description='Parareal driver, distributed version')
        parser.add_argument("tstart", type=float, help="start time")
        parser.add_argument("tend",   type=float, help="end time")
        parser.add_argument('--tol',      action="store", dest="tol",      type=float, help="tolerance for convergence, default 0.01")
        parser.add_argument('--tolcheck', action="store", dest="tolcheck",             help="method for convergence check, default L2", choices=('L2', 'maxabs'))
        parser.add_argument('--dtCoarse', action="store", dest="dtCoarse", type=float, help="coarse time increment size")
        parser.add_argument('--dtFine',   action="store", dest="dtFine",   type=float, help="fine time increment size")
        parser.add_argument('--nCoarse',  action="store", dest="nCoarse",  type=int,   help="number of coarse time increments in a coarse interval, default 1")
        parser.add_argument('--nFine',    action="store", dest="nFine",    type=int,   help="number of fine time increments in a coarse interval, default 10")
        parser.add_argument('--debug',    action="store", dest="dbg",      type=int,   help="debug printout", choices=(0, 1, 2))
        parser.add_argument('-o','--out', action="store", dest="outfile",              help="write results to FILE", metavar="FILE")
        parser.add_argument('--nsave',    action="store", dest="nSave",    type=int,   help="number of increments to save in fine sequential simulation, default 10")
        parser.add_argument('--fault',    action="store", dest="faultfile",            help="read fault definitions", metavar="FILE")
        parser.add_argument('--dist',    action="store", dest="distfile",            help="read distribution definitions", metavar="FILE")
#        parser.add_argument('--profile',  action='store_true', dest="prof",            help="profile run using cProfile")

        # Parse the command line
        args = None
        args = parser.parse_args()

        if args.dbg :
            self.dbg = args.dbg

        if args.tend <= args.tstart or args.tend<=0.0:
            if self.rank == 0:
                print("Incorrect time values", args.tstart, args.tend)
            exit()
    
        if args.nCoarse and args.dtCoarse:
            if self.rank == 0:
                print("Invalid options: Both nCoarse and dtCoarse specified.")
            exit()

        if args.nFine and args.dtFine:
            if self.rank == 0:
                print("Invalid options: Both nFine and dtFine specified.")
            exit()

        if args.dtCoarse:
            if args.dtCoarse <= 0.0:
                if self.rank == 0:
                    print("Incorrect dtCoarse increment size", args.dtCoarse)
                exit()
            self.dtCoarse=args.dtCoarse

        if args.dtFine:
            if args.dtFine <= 0.0:
                if self.rank == 0:
                    print("Incorrect dtFine increment size", args.dtFine)
                exit()
            self.dtFine=args.dtFine

        if args.nCoarse:
            if args.nCoarse <= 0:
                if self.rank == 0:
                    print("Incorrect number of coarse time steps", args.nCoarse)
                exit()
            self.nCoarse = args.nCoarse

        if args.nFine:
            if args.nFine <= 0:
                if self.rank == 0:
                    print("Incorrect number of fine time steps", args.nFine)
                exit()
            self.nFine = args.nFine

        if args.tol:
            if args.tol <=0.0:
                if self.rank == 0:
                    print("Incorrect tolerance value", args.tol)
                exit()
            self.tol=args.tol
    
        if args.tolcheck:
            if args.tolcheck != 'L2' and args.tolcheck != 'maxabs':
                if self.rank == 0:
                    print("Incorrect tolerance check value", args.tolcheck)
                exit()
            self.tolcheck=args.tolcheck

#        self.profile = args.prof

        # Time interval and increments
        self.tstart=args.tstart  # Simulation start
        self.tend  =args.tend    # Simulation end

        # Window information
        self.t1    =args.tstart      # Window start
        self.t2    =args.tend        # Window end
        timeCoarse=(self.t2-self.t1)/float(self.nProcessors)  # time interval for each processor

        # Number of increments for each processor interval
        if args.dtCoarse:
            self.dtCoarse=args.dtCoarse
            self.nCoarse=int(ceil(timeCoarse/args.dtCoarse))

        if args.dtFine:
            self.dtFine=args.dtFine
            self.nFine  =int(ceil(timeCoarse/args.dtFine))

        # Final check
        if self.nCoarse <= 0 or self.nFine <= 0:
            if self.rank == 0:
                print("Incorrect time steps", self.nCoarse, self.nFine)
            exit()

        if args.outfile:
            self.outfile=args.outfile

        if args.faultfile:
            self.faultfile=args.faultfile
            with open(self.faultfile, "r") as read_file:
                self.fault = json.load(read_file)

        if args.distfile:
            self.distfile=args.distfile
            with open(self.distfile, "r") as ntwk_file:
                self.dist = json.load(ntwk_file)

        if args.nSave:
            if args.nSave <= 0:
                print("Incorrect number of increments to save", args.nSave)
                exit()
            self.nSave = args.nSave

        if self.dbg == 1:
            if self.rank == 0:
                print("Command line arguments", args)
                if self.fault:
                    print("Number of faults", len(self.fault['faults']))
                    for ifault in self.fault['faults']:  # ifault refers to things inside fault['faults']
                        print("Fault id:", ifault['id'])
                        print("Fault start:", ifault['start'])
                        print("Fault end:", ifault['end'])
                        print("Fault type:", ifault['type'])
                        print("Fault index:", ifault['index'])
                if self.dist:
                    print("Number of Distribution systems", len(self.dist['distribution_network']))
                    for idist in self.dist['distribution_network']:  
                        print("CoSim True/False:", idist['active'])
                        print("Dist index:", idist['index'])

if __name__ == "__main__":
    import pprint
    id=0
    N=10
    md=pDefaults(id=id, N=N)
    md.parse_arguments()
    pp = pprint.PrettyPrinter(indent=4,depth=6)
    pp.pprint(vars(md))
