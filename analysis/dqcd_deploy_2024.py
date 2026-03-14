# DQCD 2024 deployment steering code

import os
import sys
sys.path.append(".")

# Configure plotting backend
import matplotlib
matplotlib.use('Agg')

from icenet.tools import process
from icedqcd import common_2024 as common
from icedqcd import deploy


def main():
    args, runmode = process.generic_flow(rootname='dqcd', func_loader=common.load_root_file, func_factor=common.splitfactor)

    deploy.process_data(args=args)

if __name__ == '__main__' :
    main()
