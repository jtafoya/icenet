# DQCD 2024 steering code

import sys
sys.path.append(".")

# Configure plotting backend
import matplotlib
matplotlib.use('Agg')

from icenet.tools import process
from icedqcd import common_2024 as common, optimize

def main():
    args, runmode = process.generic_flow(rootname='dqcd', func_loader=common.load_root_file, func_factor=common.splitfactor)

    if runmode == 'optimize':
        optimize.optimize_selection(args=args)


if __name__ == '__main__' :
    main()
