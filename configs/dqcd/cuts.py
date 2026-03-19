# Basic kinematic fiducial cuts, use only variables available in real data.
#
# m.mieskolainen@imperial.ac.uk, 2022

import awkward as ak
import numpy as np
import numba
import matplotlib.pyplot as plt

from icenet.tools import stx


def cut_nocut(X, xcorr_flow=False):
    """ No cuts
    """
    return ak.Array(np.ones(len(X), dtype=np.bool_)) # Note datatype np.bool_


def cut_fiducial(X, xcorr_flow=False):
    """ Basic fiducial (kinematic) selections.
    
    Args:
        X:          Awkward jagged array
        xcorr_flow: cut N-point cross-correlations
    
    Returns:
        Passing indices mask (N)
    """
    global O; O = X  # __technical__ recast due to eval() scope
    
    # Create cut strings
    
    names = ['ak.sum(np.logical_or(O.muonSV.mu1pt > 5.0, O.muonSV.mu2pt > 5.0), -1) > 0']
    

    #names = ['ak.sum(np.logical_and(O.Muon.pt >  5.0, np.abs(O.Muon.eta) < 2.4), -1) > 0']
             #'ak.sum(np.logical_and(O.Jet.pt  > 15.0, np.abs(O.Jet.eta)  < 2.4), -1) > 0']
             
    # 'ak.sum(O.muonSV.charge == 0, -1) > 0'
    
    # Evaluate columnar cuts; Compute cutflow
    cuts  = [eval(names[i], globals()) for i in range(len(names))]
    mask  = stx.apply_cutflow(cut=cuts, names=names, xcorr_flow=xcorr_flow)

    return mask


def cut_fiducial_2024(X, xcorr_flow=False):
    """ Basic fiducial (kinematic) selections for 2024 data.
        OR of HLT_Mu10_Barrel_L1HP11_IP6 and HLT_DoubleMu4_3_LowMass legs.

    Args:
        X:          Awkward jagged array
        xcorr_flow: cut N-point cross-correlations

    Returns:
        Passing indices mask (N)
    """
    global O; O = X  # __technical__ recast due to eval() scope

    # Create cut strings
    # OR of the two 2024 HLT paths:
    #   HLT_Mu10_Barrel_L1HP11_IP6:   at least one muon with pT > 10 GeV in the barrel (|eta| < 1.2)
    #   HLT_DoubleMu4_3_LowMass:      leading muon pT > 4 GeV and subleading pT > 3 GeV

    names = [
        # HLT_Mu10_Barrel_L1HP11_IP6 leg
        '(ak.sum('
        '  np.logical_or('
        '    np.logical_and(O.muonSV.mu1pt > 10.0, np.abs(O.muonSV.mu1eta) < 1.2),'
        '    np.logical_and(O.muonSV.mu2pt > 10.0, np.abs(O.muonSV.mu2eta) < 1.2)'
        '  ), -1) > 0)'
        ' | '
        # HLT_DoubleMu4_3_LowMass leg
        '(ak.sum('
        '  np.logical_and('
        '    np.maximum(O.muonSV.mu1pt, O.muonSV.mu2pt) > 4.0,'
        '    np.minimum(O.muonSV.mu1pt, O.muonSV.mu2pt) > 3.0'
        '  ), -1) > 0)'
    ]

    # Evaluate columnar cuts; Compute cutflow
    cuts  = [eval(names[i], globals()) for i in range(len(names))]
    mask  = stx.apply_cutflow(cut=cuts, names=names, xcorr_flow=xcorr_flow)

    return mask


def cut_fiducial_Mu10(X, xcorr_flow=False):
    """ Fiducial cuts for HLT_Mu10_Barrel_L1HP11_IP6:
        at least one muon with pT > 10 GeV in the barrel (|eta| < 1.2).

    Args:
        X:          Awkward jagged array
        xcorr_flow: cut N-point cross-correlations

    Returns:
        Passing indices mask (N)
    """
    global O; O = X  # __technical__ recast due to eval() scope

    names = [
        'ak.sum('
        '  np.logical_or('
        '    np.logical_and(O.muonSV.mu1pt > 10.0, np.abs(O.muonSV.mu1eta) < 1.2),'
        '    np.logical_and(O.muonSV.mu2pt > 10.0, np.abs(O.muonSV.mu2eta) < 1.2)'
        '  ), -1) > 0'
    ]

    # Evaluate columnar cuts; Compute cutflow
    cuts  = [eval(names[i], globals()) for i in range(len(names))]
    mask  = stx.apply_cutflow(cut=cuts, names=names, xcorr_flow=xcorr_flow)

    return mask


def cut_fiducial_DoubleMu(X, xcorr_flow=False):
    """ Fiducial cuts for HLT_DoubleMu4_3_LowMass:
        leading muon pT > 4 GeV and subleading pT > 3 GeV.

    Args:
        X:          Awkward jagged array
        xcorr_flow: cut N-point cross-correlations

    Returns:
        Passing indices mask (N)
    """
    global O; O = X  # __technical__ recast due to eval() scope

    names = [
        'ak.sum('
        '  np.logical_and('
        '    np.maximum(O.muonSV.mu1pt, O.muonSV.mu2pt) > 4.0,'
        '    np.minimum(O.muonSV.mu1pt, O.muonSV.mu2pt) > 3.0'
        '  ), -1) > 0'
    ]

    # Evaluate columnar cuts; Compute cutflow
    cuts  = [eval(names[i], globals()) for i in range(len(names))]
    mask  = stx.apply_cutflow(cut=cuts, names=names, xcorr_flow=xcorr_flow)

    return mask
