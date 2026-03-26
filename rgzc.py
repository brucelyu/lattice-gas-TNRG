#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : findzc.py
# Author            : Xinliang (Bruce) Lyu <xlyu@ihes.fr>
# Date              : 14.05.2025
# Last Modified Date: 14.05.2025
# Last Modified By  : Xinliang (Bruce) Lyu <xlyu@ihes.fr>
import argparse
from datetime import datetime
from dateutil.relativedelta import relativedelta
from wilsonRG import findTc, rgzcflow

# argument parser
parser = argparse.ArgumentParser(
    "Determine critical parameter of a given model"
)

parser.add_argument("--model", dest="model", type=str,
                    help="statsphys model",
                    choices=["ising", "hardsquare1NN", "hardsquareNeg"],
                    default="hardsquare1NN")
parser.add_argument("--scheme", dest="scheme", type=str,
                    help="TNRG scheme",
                    choices=["trg", "efrg", "loop-TNR"],
                    default="trg")
parser.add_argument("--ver", dest="ver", type=str,
                    help="TNRG scheme version",
                    default="general")
parser.add_argument("--chi", dest="chi", type=int,
                    help="bound dimension (default: 10)",
                    default=10)
parser.add_argument("--cgeps", dest="cgeps", type=float,
                    help="a number smaller than which we think the" +
                    " singluar values for the environment" +
                    " in RG spectrum is zero (default: 1e-12)",
                    default=1e-12)
parser.add_argument("--rgn", dest="rgn", type=int,
                    help="maximal RG iteration (default: 35)",
                    default=35)
parser.add_argument("--outDir", type=str,
                    help="output directory to save rg flows and Tc",
                    default="./")
parser.add_argument("--chis", dest="chis", type=int,
                    help="squeezed bond dimension (default: 8)",
                    default=8)
parser.add_argument("--epspinv", dest="epspinv", type=float,
                    help="cutoff value for Moore-Penrose inverse" +
                    " in Loop Optimization (default: 1e-5)",
                    default=1e-5)

# read from argument parser
args = parser.parse_args()

model = args.model
scheme = args.scheme
ver = args.ver

chi = args.chi
cgeps = args.cgeps
chis = args.chis

n_RG = args.rgn
outDir = args.outDir


if scheme == "trg":
    tnrg_pars = {"chi": chi, "dtol": cgeps, "display": True}
elif scheme == "efrg":
    tnrg_pars = {"chi": chi, "dtol": cgeps, "display": True,
                 "chis": chis, "chienv": chis**2, "epsilon": 1e-8}
elif scheme == "loop-TNR":
    LF_max = 200
    LF_min = 10
    softinv = False
    eps_pinv = cgeps * 1.0
    if model == "hardsquareNeg":
        softinv = True
        eps_pinv = args.epspinv
    tnrg_pars = {"chi": chi, "dtol": cgeps, "display": True,
                 "eps_pinv": eps_pinv, "eps_errEF": cgeps * 1.0,
                 "LF_max": LF_max, "LF_min": LF_min, "softinv": softinv}

print("The model:         --{:s}--".format(model))
print("TNRG scheme:       --{:s}--".format(scheme))
print("    with version   --{:s}--".format(ver))
print("Bond dimension:    --{:d}--".format(chi))
if scheme == "efrg":
    print("  The squeezed χs: --{:d}--".format(chis))
elif scheme == "loop-TNR":
    print("      For the loop optimization:")
    print("       ε in pinv of Y environment (softinv={}):   --{:.2e}--.".format(softinv, eps_pinv))

start_now = datetime.now()
current_time = start_now.strftime("%Y-%m-%d. %H:%M:%S")
print("Running Time is", current_time)

rgzcflow(
    model, scheme, ver, tnrg_pars,
    n_RG, outDir="./"
)

now = datetime.now()
current_time = now.strftime("%Y-%m-%d. %H:%M:%S")
print("----------------------------------")
print("Finished Time is", current_time)
diffT = relativedelta(now, start_now)
print("Elapsed wall time is {} days {} hours {} minutes {} seconds.".format(
    diffT.days, diffT.hours, diffT.minutes, diffT.seconds
)
      )


# end of file
