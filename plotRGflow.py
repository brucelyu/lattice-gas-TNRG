#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plotRGflow.py
# Author            : Xinliang (Bruce) Lyu <xlyu@ihes.fr>
# Date              : 26.06.2025
# Last Modified Date: 26.06.2025
# Last Modified By  : Xinliang (Bruce) Lyu <xlyu@ihes.fr>
import argparse
from wilsonRG import saveDirName
import pickle as pkl
import matplotlib.pyplot as plt
import numpy as np


def plotscaleD(rgsteps, scDList, saveDir, scheme, ver, tnrg_pars):
    xShift = [0, 0, 0, -0.05, 0.05, -0.15, -0.05, 0.05, 0.15, -0.1, 0, 0.1]
    D_Color = ["k", "k", "k", "b", "b", "k", "k", "k", "k", "b", "b", "b"]
    D_Marker = [".", "x", "*", "+", "+", "p", "p", "p", "p", "+", "+", "+"]
    fontsz = 12
    plt.figure(figsize=(8, 4))
    for rgn, scD in zip(rgsteps, scDList):
        for k in range(12):
            plt.plot(rgn + xShift[k], scD[k],
                     D_Color[k] + D_Marker[k])
        # relative errors
        sExt = 0.125
        sErr = np.abs(scD[1] - sExt) / sExt
        plt.text(rgn - 0.20, sExt + 0.10, "{:.2%}".format(sErr), size=fontsz)
        eExt = 1.0
        eErr = np.abs(scD[2] - eExt) / eExt
        plt.text(rgn - 0.20, eExt - 0.20, "{:.2%}".format(eErr), size=fontsz)
        dsExt = 1.125
        dsErr1 = np.abs(scD[3] - dsExt) / dsExt
        dsErr2 = np.abs(scD[4] - dsExt) / dsExt
        dsErrAve = (dsErr1 * dsErr2)**(1/2)
        # plt.text(rgn - 0.20, eExt + 0.30, "{:.2%}".format(dsErr1), size=14)
        plt.text(rgn - 0.20, eExt + 0.25, "{:.3%}".format(dsErrAve), size=fontsz)

        TExt = 2.0
        TErr1 = np.abs(scD[5] - TExt) / TExt
        TErr2 = np.abs(scD[6] - TExt) / TExt
        TErr3 = np.abs(scD[7] - TExt) / TExt
        TErr4 = np.abs(scD[8] - TExt) / TExt
        TErrAve = (TErr1 * TErr2 * TErr3 * TErr4)**(1/4)
        plt.text(rgn - 0.20, TExt - 0.20, "{:.2%}".format(TErrAve), size=fontsz)

        ddsExt = 2.125
        ddsErr1 = np.abs(scD[9] - ddsExt) / ddsExt
        ddsErr2 = np.abs(scD[10] - ddsExt) / ddsExt
        ddsErr3 = np.abs(scD[11] - ddsExt) / ddsExt
        ddsErrAve = (ddsErr1 * ddsErr2 * ddsErr3)**(1/3)
        plt.text(rgn - 0.20, ddsExt + 0.15, "{:.2%}".format(ddsErrAve), size=fontsz)

    # best-known values
    plt.hlines(0.125, rgsteps[0]-0.2, rgsteps[-1]+0.2,
               colors="black", linestyles="solid", alpha=0.2)
    plt.hlines(1.0, rgsteps[0]-0.2, rgsteps[-1]+0.2,
               colors="black", linestyles="solid", alpha=0.2)
    plt.hlines(1.125, rgsteps[0]-0.2, rgsteps[-1]+0.2,
               colors="blue", linestyles="solid", alpha=0.2)
    plt.hlines(2.0, rgsteps[0]-0.2, rgsteps[-1]+0.2,
               colors="black", linestyles="solid", alpha=0.2)
    plt.hlines(2.125, rgsteps[0]-0.2, rgsteps[-1]+0.2,
               colors="blue", linestyles="solid", alpha=0.2)
    # set axis ranges
    plt.xticks(rgsteps)
    if len(rgsteps) == 1:
        plt.xlim([rgsteps[0] - 1, rgsteps[0] + 1])
    plt.ylabel("Scaling dimensions", size=fontsz)
    plt.ylim([-0.1, 2.5])
    plt.xlabel(r"RG step $n$", size=fontsz)
    plt.xticks(fontsize=fontsz)
    plt.yticks(fontsize=fontsz)
    if scheme == "efrg":
        plt.title(r"{:s} ({:s}) with $\chi=${:d} ($\chi_s=${:d})".format(
            scheme, ver, chi, chis)
                  )
    elif scheme in ["trg", "loop-TNR"]:
        plt.title(r"{:s} ({:s}) with $\chi=${:d}".format(scheme, ver, chi),
                  size=fontsz)

    plt.savefig(saveDir + "/scDim.png",
                bbox_inches='tight', dpi=300)


def plotscaleDNeg(rgsteps, scDList, saveDir, scheme, ver, tnrg_pars):
    xShift = [0, -0.05, 0.05, -0.10, 0, 0.10, -0.10,  0.10]
    D_Color = ["k", "b", "b", "b", "b", "b", "k", "k"]
    D_Marker = [".", "+", "+", "+", "+", "+", "p", "p"]
    plt.figure(figsize=(8, 4))
    fontsz = 12
    for rgn, scD in zip(rgsteps, scDList):
        for k in range(8):
            plt.plot(rgn + xShift[k], scD[k+1],
                     D_Color[k] + D_Marker[k], markersize=8)
        # relative errors
        sExt = 2/5
        sErr = np.abs(scD[1] - sExt) / np.abs(sExt)
        plt.text(rgn - 0.20, sExt - 0.20, "{:.2%}".format(sErr), size=fontsz)
        dsExt = 1
        dsErr1 = np.abs(scD[2] - dsExt) / np.abs(dsExt)
        dsErr2 = np.abs(scD[3] - dsExt) / np.abs(dsExt)
        dsErrAve = (dsErr1 * dsErr2)**(1/2)
        plt.text(rgn - 0.20, dsExt - 0.20, "{:.2%}".format(dsErrAve), size=fontsz)
        # plt.text(rgn - 0.20, dsExt - 0.30, "{:.1e}".format(dsErr2), size=fontsz)
        ddsExt = 2
        ddsErr1 = np.abs(scD[4] - ddsExt) / ddsExt
        ddsErr2 = np.abs(scD[5] - ddsExt) / ddsExt
        ddsErr3 = np.abs(scD[6] - ddsExt) / ddsExt
        ddsErrAve = (ddsErr1 * ddsErr2 * ddsErr3)**(1/3)
        # plt.text(rgn - 0.20, ddsExt - 0.20, "{:.1e}".format(ddsErr1), size=fontsz)
        plt.text(rgn - 0.20, ddsExt - 0.30, "{:.2%}".format(ddsErrAve), size=fontsz)
        # plt.text(rgn - 0.20, ddsExt - 0.40, "{:.1e}".format(ddsErr3), size=fontsz)
        TExt = 2.4
        TErr1 = np.abs(scD[7] - TExt) / TExt
        TErr2 = np.abs(scD[8] - TExt) / TExt
        TErrAve = (TErr1 * TErr2)**(1/2)
        plt.text(rgn - 0.20, TExt + 0.25, "{:.2%}".format(TErrAve), size=fontsz)
        # plt.text(rgn - 0.20, TExt + 0.10, "{:.1e}".format(TErr2), size=fontsz)
    # best-known values
    plt.hlines(0.4, rgsteps[0]-0.2, rgsteps[-1]+0.2,
               colors="black", linestyles="solid", alpha=0.2)
    plt.hlines(1.0, rgsteps[0]-0.2, rgsteps[-1]+0.2,
               colors="blue", linestyles="solid", alpha=0.2)
    plt.hlines(2.0, rgsteps[0]-0.2, rgsteps[-1]+0.2,
               colors="blue", linestyles="solid", alpha=0.2)
    plt.hlines(2.4, rgsteps[0]-0.2, rgsteps[-1]+0.2,
               colors="black", linestyles="solid", alpha=0.2)
    # set axis ranges
    plt.xticks(rgsteps)
    if len(rgsteps) == 1:
        plt.xlim([rgsteps[0] - 1, rgsteps[0] + 1])
    plt.ylabel("Scaling dimensions")
    plt.xlabel(r"RG step $n$")
    plt.ylabel("Scaling dimensions", size=fontsz)
    plt.xlabel(r"RG step $n$", size=fontsz)
    plt.xticks(fontsize=fontsz)
    plt.yticks(fontsize=fontsz)
    if scheme == "efrg":
        plt.title(r"{:s} ({:s}) with $\chi=${:d} ($\chi_s=${:d})".format(
            scheme, ver, chi, chis)
                  )
    elif scheme in ["trg", "loop-TNR"]:
        plt.title(r"{:s} ({:s}) with $\chi=${:d}".format(scheme, ver, chi))

    plt.ylim(0.0, 3.0)
    plt.savefig(saveDir + "/scDim.png",
                bbox_inches='tight', dpi=200)


# argument parser
parser = argparse.ArgumentParser(
    "Plot various RG flows"
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
parser.add_argument("--outDir", type=str,
                    help="output directory to save rg flows and Tc",
                    default="./")
parser.add_argument("--chis", dest="chis", type=int,
                    help="squeezed bond dimension (default: 8)",
                    default=8)
parser.add_argument("--startn", type=int,
                    help="starting RG step (default: 1)",
                    default=1)
parser.add_argument("--endn", type=int,
                    help="ending RG step (default: 9)",
                    default=9)

# read from argument parser
args = parser.parse_args()

model = args.model
scheme = args.scheme
ver = args.ver

chi = args.chi
cgeps = 1e-10
chis = args.chis

outDir = args.outDir

startn = args.startn
endn = args.endn

if scheme == "trg":
    tnrg_pars = {"chi": chi, "dtol": cgeps, "display": True}
elif scheme == "efrg":
    tnrg_pars = {"chi": chi, "dtol": cgeps, "display": True,
                 "chis": chis, "chienv": chis**2, "epsilon": 1e-8}
elif scheme == "loop-TNR":
    LF_max = 50
    LF_min = 5
    eps_pinv = 1e-5
    tnrg_pars = {"chi": chi, "dtol": cgeps, "display": True,
                 "eps_pinv": eps_pinv, "eps_errEF": cgeps * 1.0,
                 "LF_max": LF_max, "LF_min": LF_min}

# read the RG flow data
dataDir = saveDirName(model, scheme, ver, tnrg_pars, outDir)
tenDir = dataDir + "/tensors"

# I. Plot flow of scaling dimensions
scDfile = tenDir + "/ScD.pkl"
with open(scDfile, "rb") as f:
    rgnList, scDflow = pkl.load(f)

if model == "hardsquare1NN":
    plotscaleD(rgnList[startn-1:endn], scDflow[startn-1:endn],
               dataDir, scheme, ver, tnrg_pars)
elif model == "hardsquareNeg":
    plotscaleDNeg(
        rgnList[startn-1:endn], scDflow[startn-1:endn],
        dataDir, scheme, ver, tnrg_pars
    )

# II. Plot RG flow of errors
rgflowsfile = tenDir + "/tenflows.pkl"
with open(rgflowsfile, "rb") as f:
    rgnList, SPerrsFlow, lrerrsFlow = pkl.load(f)

SPerrsFlow = np.array(SPerrsFlow)


if scheme == "trg":
    fig, axs = plt.subplots(1, 1, figsize=(5, 2))
elif scheme == "loop-TNR":
    fig, axs = plt.subplots(2, 1, figsize=(5, 4))

if scheme == "trg":
    axs.plot(rgnList[startn-1:endn], SPerrsFlow[startn-1:endn, 1],
             "kx-", alpha=0.4
             )
    axs.set_ylabel("2nd TRG error")
    axs.set_title(r"$\chi = ${:d}".format(chi))
    axs.set_yscale("log")
    axs.set_xlabel("RG step")


if scheme == "loop-TNR":
    ax1 = axs[0]
    ax1.plot(rgnList[startn-1:endn], SPerrsFlow[startn-1:endn, 1],
             "kx-", alpha=0.4
             )
    ax1.set_ylabel("2nd TRG error")
    ax1.set_title(r"$\chi = ${:d}".format(chi))
    ax1.set_yscale("log")
    ax2 = axs[1]
    ax2.plot(rgnList[startn-1:endn], lrerrsFlow[startn-1:endn], "g.--")
    ax2.set_ylabel("Loop optimization error")
    ax2.set_yscale("log")
    ax2.set_xlabel("RG step")

plt.savefig(dataDir + "/RGerrs.png", bbox_inches='tight', dpi=200)


# end of file
