#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : wilsonRG.py
# Author            : Xinliang (Bruce) Lyu <xlyu@ihes.fr>
# Date              : 14.05.2025
# Last Modified Date: 14.05.2025
# Last Modified By  : Xinliang (Bruce) Lyu <xlyu@ihes.fr>
"""
Functions for Wilsonian RG prescription
- determine Tc by generating tensor RG flows
"""
import os
import itertools
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
from tensornetworkrg import tnrg


# I. (general purpose) For generating tensor RG flows
# I.1 create model instant and print out all informations
def genRGflow(
    model, scheme, ver, tnrg_pars, n_RG, model_pars,
    isOnSiteSym=True, dataDir=None, determPhase=True,
    Xhi=1, Xlow=2
):
    # create model instance
    modIns = tnrg.TensorNetworkRG2D(model)
    if model == "hardsquareNeg":
        modIns = tnrg.TensorNetworkRG2D("hardsquare1NN")

    modIns.model_parameters = model_pars.copy()
    if model == "hardsquare1NN":
        modIns.generate_initial_tensor(scheme="trg")
    elif model in ["ising2d", "ising3d"]:
        modIns.generate_initial_tensor(onsite_symmetry=isOnSiteSym)
    elif model == "hardsquareNeg":
        if ver == "general":
            # no lattice symmetry; so no bond matrix (like the general TRG)
            modIns.generate_initial_tensor(scheme="trgR")
        elif ver == "rotsym":
            # rotation symmetry representation with bond matrix
            modIns.generate_initial_tensor(scheme="trgRsym")


    # display info of the algorithm and model
    if tnrg_pars["display"]:
        print("The model:           --{:s}--".format(model))
        print("The basic TNRG parameters are")
        print("TNRG bond dimension: --{:d}--".format(tnrg_pars["chi"]))
        print("TNRG epsilon:        --{:.2e}--".format(tnrg_pars["dtol"]))
        print("TNRG iteration step: --{:d}--".format(n_RG))
        print("Onsite symmetry:     --{}--".format(isOnSiteSym))
        print("------")
        # print additional parameters
        if scheme == "efrg":
            print("    ",
                  "EFRG ({:s}) is applied...".format(ver))
            print("    ",
                  "Parameters for EF process:")
            print("    ", "- chis: --{:d}--".format(tnrg_pars["chis"]))
            print("    ", "- chienv: --{:d}--".format(tnrg_pars["chienv"]))
            print("    ", "- epsilon: --{:.2e}--".format(tnrg_pars["epsilon"]))
        print("_/_/_/_/_/_/_/_/_/")
    # generate degeneracy index X flow
    (
        XFlow, SPerrsFlow, lrerrsFlow
    ) = rgIterate(modIns, n_RG, scheme, ver,
                  tnrg_pars=tnrg_pars, dataDir=dataDir,
                  determPhase=determPhase, Xhi=Xhi, Xlow=Xlow)
    return XFlow, SPerrsFlow, lrerrsFlow


# I.2 Apply the RG map repeatedly
def rgIterate(
    modIns, n_RG, scheme, ver, tnrg_pars,
    dataDir=None, determPhase=False, Xhi=1, Xlow=2
):
    XFlow = []
    lrerrsFlow = []
    SPerrsFlow = []
    # save the data related to the initial tensor
    if dataDir is not None:
        # list to save the scaling dimensions a la Gu & Wen
        scDflow = []
        rgnList = []
        # save the initial tensor
        Aorg = modIns.get_tensor()
        tenDir = dataDir + "/tensors"
        fname = "A00.pkl"
        saveData(tenDir, fname,
                 data=Aorg)

    # apply the RG repeatedly
    for k in range(n_RG):
        if tnrg_pars["display"]:
            print()
            print("Applying the RG map...")

        # apply the RG map
        (lrerrs, SPerrs) = modIns.rgmap(tnrg_pars, scheme=scheme, ver=ver)
        # various properties of the current tensor
        curX = modIns.degIndX2()
        if modIns.model == "hardsquare1NN":
            if modIns.get_model_parameters()['activity'] < 0:
                curX = modIns.degIndX2()
        XFlow.append(curX)
        lrerrsFlow.append(lrerrs)
        SPerrsFlow.append(SPerrs)
        if tnrg_pars["display"]:
            print("===")
            print("The RG step {:d} finished!".format(modIns.get_iteration()))
            print("Shape of A is")
            print(" {}".format(modIns.get_tensor().shape))
            print("Degnerate index X = {:.2f}".format(curX))
            print("----------")
            print("----------")
        # exit the RG iteration if already flowed to the trivial phase
        if determPhase and k > 1:
            stop_eps = 0.01
            near1 = (
                (abs(XFlow[-1] - Xhi) < stop_eps) and
                (abs(XFlow[-2] - Xhi) < stop_eps)
            )
            near2 = (
                (abs(XFlow[-1] - Xlow) < stop_eps) and
                (abs(XFlow[-2] - Xlow) < stop_eps)
            )
            if near1 or near2:
                break

        # for saving the RG flows
        if dataDir is not None:
            fname = "A{:02d}.pkl".format(k + 1)
            Anew = modIns.get_tensor()
            saveData(tenDir, fname,
                     data=Anew)
            # extract scaling dimensions a la Gu, Wen and Cardy
            if modIns.model == "hardsquare1NN":
                if modIns.get_model_parameters()['activity'] < 0:
                    ccharge, scDim = modIns.gu_wen_cardy(
                        aspect_ratio=2, num_scale=17, indId=0
                    )
                else:
                    ccharge, scDim = modIns.gu_wen_cardy(
                        aspect_ratio=2, num_scale=17
                    )
            if tnrg_pars["display"]:
                print("Scaling dimensions are")
                with np.printoptions(precision=6, suppress=True):
                    print(scDim)
                print("----------")
                print("----------")
            scDflow.append(scDim)
            rgnList.append(modIns.get_iteration())

    # save scaling dimensions (outside the RG iteration) and other flows
    if dataDir is not None:
        fnameScD = "ScD.pkl"
        fname = "tenflows.pkl"
        saveData(
            tenDir, fnameScD,
            data=[rgnList, scDflow]
        )
        saveData(
            tenDir, fname,
            data=[rgnList, SPerrsFlow, lrerrsFlow]
        )

    return XFlow, SPerrsFlow, lrerrsFlow


# II. For determining the critical parameter
# II.1 Tell whether the degeneracy index flows to ssb or sym phase
def ishiT(X, Xhi, Xlow):
    dist2hi = np.abs(X - Xhi)
    dist2low = np.abs(X - Xlow)
    if dist2hi < dist2low:
        res = True
    else:
        res = False
    return res


# II.2 The main function for determing the critical parameter
def findTc(
    model, scheme, ver, tnrg_pars,
    Tlow, Thi, n_RG, n_bisect,
    outDir="./", Xhi=1, Xlow=2
):
    # directory name for saving the output
    saveDir = saveDirName(model, scheme, ver, tnrg_pars, outDir)
    if (not os.path.exists(saveDir)):
        os.makedirs(saveDir)
    # read Tc if it exists
    TcFile = saveDir + "/Tc.pkl"
    if os.path.exists(TcFile):
        with open(TcFile, "rb") as f:
            Tlow, Thi = pkl.load(f)

    #   - generate two flows of degeneracy index X
    model_pars_hi = getModelPara(model, Thi)
    model_pars_lo = getModelPara(model, Tlow)
    XhiFlow = genRGflow(
        model, scheme, ver, tnrg_pars, n_RG, model_pars_hi,
        determPhase=True, Xhi=Xhi, Xlow=Xlow
    )[0]
    XloFlow = genRGflow(
        model, scheme, ver, tnrg_pars, n_RG, model_pars_lo,
        determPhase=True, Xhi=Xhi, Xlow=Xlow
    )[0]

    # The model at Tlow and Thi should flow to the corresponding fixed point
    errMesgHi = "The model with Thi should flow to the high-T fixed point"
    assert ishiT(XhiFlow[-1], Xhi, Xlow), errMesgHi
    errMesgLo = "The model with Tlow should flow to the low-T fixed point"
    assert not ishiT(XloFlow[-1], Xhi, Xlow), errMesgLo

    # for plotting the flows of degeneracy index X
    plt.figure(figsize=[8, 6])
    lines = itertools.cycle(("--", "-.", "-"))
    # enter the bisection \\\\\\
    for k in range(n_bisect + 1):
        # generate trial degeneracy index X flow
        Ttry = 0.5 * (Tlow + Thi)
        model_pars_try = getModelPara(model, Ttry)
        XtryFlow = genRGflow(
            model, scheme, ver, tnrg_pars, n_RG, model_pars_try,
            determPhase=True, Xhi=Xhi, Xlow=Xlow
        )[0]
        # plot every 3 iteration
        if (k % 3 == 0):
            curline = next(lines)
            saveFigName = saveDir + "/X_iterk{:02d}.png".format(k+1)
            plotXFlows(
                XloFlow, XhiFlow, XtryFlow, Tlow, Thi, Ttry,
                saveFigName, curline, alpha=(k+1)/(n_bisect+1),
                Xhi=Xhi, Xlow=Xlow
            )
            print("This is the {:d}-th bisection step...".format(k + 1))
            print("Critical parameter is bounded by",
                  "{:.8f} and {:.8f}".format(Tlow, Thi)
                  )
            print("---")
        # update low and high bound
        if ishiT(XtryFlow[-1], Xhi, Xlow):
            Thi = Ttry * 1.0
            XhiFlow = XtryFlow.copy()
        else:
            Tlow = Ttry * 1.0
            XloFlow = XtryFlow.copy()
    # outside the bisection iteration //////
    # save the lower and upper bound of Tc
    with open(TcFile, "wb") as f:
        pkl.dump([Tlow, Thi], f)
    # Append all figures
    orgfile = saveDir + '/X_iterk*.png'
    tarfile = saveDir + '/Xflow_all.png'
    syscommand = 'magick ' + orgfile + " -append " + tarfile
    if os.system(syscommand) != 0:
        print("Command, convert, not found in current os")
    else:
        os.system('rm ' + orgfile)


# II.3 Helper functions for the main function `findTc` above
def saveDirName(model, scheme, ver, pars, outDir):
    """generate directory name for saving
    """
    modelDir = (
        outDir + "{:s}_out/{:s}_{:s}/".format(model, scheme, ver)
    )

    if scheme == "efrg":
        # For entanglement-filtering-enhanced methods
        paraDir = paraDir = "chi{:02d}s{:d}".format(
                pars["chi"], pars["chis"]
        )
    else:
        # For schemes without filtering
        paraDir = "chi{:02d}".format(pars["chi"])

    saveDir = modelDir + paraDir
    return saveDir


def getModelPara(model, T):
    if model in ["ising2d", "ising3d"]:
        model_pars = {"temperature": T, "magnetic_field": 0}
    elif model in ["hardsquare1NN", "hardsquareNeg"]:
        model_pars = {"activity": T}
    return model_pars


def plotXFlows(
    XlowFlow, XhiFlow, XtryFlow, Tlow, Thi, Ttry,
    savefile, linesty, alpha=1, Xhi=1, Xlow=2
):
    # precision of the estimate
    precisionT = abs(Thi - Ttry) / Ttry
    rg_n = max(len(XlowFlow), len(XhiFlow), len(XtryFlow))
    plt.title(("Precision of the estimate is {:.2e}, ".format(precisionT)))
    # plot three flows of the degeneracy index X
    plt.plot(
        XlowFlow[0:], "bo" + linesty, alpha=alpha,
        label=r"$z_{{ssb}}$ = {:.5f}".format(Tlow)
    )
    plt.plot(
        XhiFlow[0:], "k." + linesty, alpha=alpha,
        label=r"$z_{{sym}}$ = {:.5f}".format(Thi)
    )
    plt.plot(
        XtryFlow[0:], "gx" + linesty, alpha=alpha,
        label=r"$z_{{try}}$ = {:.5f}".format(Ttry)
    )
    # set plot info
    plt.ylim([Xhi - 0.1, Xlow + 0.1])
    plt.hlines(Xhi, 0, rg_n-1, color="k", linestyles='dashed', alpha=0.2)
    plt.hlines(Xlow, 0, rg_n-1, color="k", linestyles='dashed', alpha=0.2)
    plt.legend()
    plt.xlabel(r"RG step $n$")
    plt.ylabel("Degeneracy index $X$")
    plt.savefig(savefile, dpi=300)


# III. For generating the RG flow at estimated zc
def rgzcflow(
    model, scheme, ver, tnrg_pars, n_RG,
    outDir="./"
):
    """
    Generate the RG flow at the estimated zc
    """
    # read zc file
    saveDir = saveDirName(model, scheme, ver, tnrg_pars, outDir)
    TcFile = saveDir + "/Tc.pkl"
    with open(TcFile, "rb") as f:
        zlow, zhi = pkl.load(f)
        zc = 0.5 * (zlow + zhi)

    if model == "hardsquare1NN":
        zcBest = 3.79625517391234
    elif model == "hardsquareNeg":
        zcBest = -0.119338886
    zcerr = np.abs(zc - zcBest) / np.abs(zcBest)
    print("=======")
    print("Estimate zc = {:.8f} (err = {:.2e})".format(zc, zcerr))
    print("=======")
    # generate the RG flow at zc
    model_pars_c = getModelPara(model, zc)
    (
        XFlow, SPerrsFlow, lrerrsFlow
    ) = genRGflow(
        model, scheme, ver, tnrg_pars, n_RG, model_pars_c,
        isOnSiteSym=True, dataDir=saveDir, determPhase=False
    )


# III.1 For save RG flows
def saveData(outDir, fname, data=None):
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    saveFile = outDir + "/" + fname
    with open(saveFile, "wb") as f:
        pkl.dump(data, f)
