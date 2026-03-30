# A TNRG study of the 1NN hard-square lattice gas model
In this repository, we used a TNRG scheme enhanced by loop optimization to study the hard-square lattice gas model with nearest-neighbor exclusion (1NN).
The scheme we use here incorporates the lattice-reflection, lattice-rotation and $\mathcal{PT}$ symmetry of the model.
These symmetries are spontaneously broken in the two transitions of the model.
In our preprint, [Lattice and PT symmetries in tensor-network renormalization group: a case study of a hard-square lattice gas model](https://arxiv.org/abs/2603.25492), we report the design of the scheme and the details of the numerical calculations conducted below. 
The scheme is a generalization of the symmetric version of an advanced TNRG scheme called [loop-TNR](https://arxiv.org/abs/1512.04938).


## 0. Requirements
Apart from the standard tool kit for scientific computating for python, [Ananconda](https://www.anaconda.com/download), we need the following packages for tensor network manipulations:
- [ncon](https://github.com/mhauru/ncon), for tensor contractions,
- [abeliantensors](https://github.com/mhauru/abeliantensors), for various tensor manipulations, and
- [tn-tools](https://github.com/mhauru/tntools).

The implementation of the tensor RG map is in another repository [tensornetworkrg](https://github.com/brucelyu/tensornetworkrg), which is included here as a submodule.
After cloning this repository to your computer, remember using the following command to pull from the submodules:
```
git submodule update --init --recursive
```


## I. Description of the files
### I.1 Scripts
We will use three scripts to do the numerical calculations:
1. `findzc.py` determines the critical activity of the model using a bisection method.
1. `rgzc.py` generates a tensor RG flow at the estimated critical activity; scaling dimensions are also estimated in this script.
1. `plotRGflow.py` plots RG flows of estimated scaling dimensions, as well as the RG approximation errors.

The file `wilsonRG.py` contains functions that are used in the above three scripts.


### I.2 Jupyter Notebooks
The three notebooks in `jupyterNB` explore the stability of the SSB fixed points, as well as the PT symmetry breaking of the transfer matrix:
- `Sym-SSB-zplus.ipynb` studies the SSB fixed point of the positive-activity transition.
- `Sym-SSB-zminus.ipynb` studies the SSB fixed point of the negative-activity transition.
- `TMspectz.ipynb` studies the PT symmetry breaking of the transfer matrix for the negative-activity transition.


### I.3 Implementation of the RG map
We use the short hand notation `tnrg` to denote the library `tensornetworkrg` that is included here as a submodule.
***The proposed RG map is implemented as the method `tnrg.TensorNetworkRG2D.trg_loopOpt` of the class `tnrg.TensorNetworkRG2D`.***
This method assembles functions in the following two files:
- TODO
- TODO


## II. Apply the EF-enhanced TNRG to the model
### II.1 For positive-activity transition
*Step 1:*
At bond dimension $\chi=10$, use bisection method to find the critical activity $z_c^{+}$:
```python
python findzc.py --model hardsquare1NN --scheme loop-TNR --ver rotsym --chi 10 --rgn 26 --itern 9 --zhi 3.6 --zlo 4.0
```
This script conducts 9 iterations of the bisection with the lower and upper bounds of $z$ being 3.6 and 4.0.
The estimated value of $z_c^{+}$ will be save in the file `./hardsquare1NN_out/loop-TNR_rotsym/chi10/Tc.pkl`.
This script can be run repeatedly to make the window of the estimate smaller.

*Step 2:*
Generate the tensor RG flow at the estimated $z_c^{+}$:
```python
python rgzc.py --model hardsquare1NN --scheme loop-TNR --ver rotsym --chi 10 --rgn 12
```
The tensor RG flow will be save in the folder `./hardsquare1NN_out/loop-TNR_rotsym/chi10/tensors`.
The estimated scaling dimensions will be save in the file `./hardsquare1NN_out/loop-TNR_rotsym/chi10/tensors/ScD.pkl`.

*Step 3:*
Plot RG flows of scaling dimensions and the RG approximation errors from RG step $n=3-9$:
```python
 python plotRGflow.py --model hardsquare1NN --scheme loop-TNR --ver rotsym --chi 10 --startn 3 --endn 9
```
The figures are save in the folder `./hardsquare1NN_out/loop-TNR_rotsym/chi10`.

### II.2 For negative-activity transition
The procedure is the same as the positive-activity transition above.
Choose bond dimension $\chi=10$.

*Step 1:* Find $z_c^{-}$
```python
python findzc.py --model hardsquareNeg --scheme loop-TNR --ver rotsym --chi 10 --epspinv 1e-8 --rgn 26 --itern 9 --zhi -0.119 --zlo -0.121
```

*Step 2:* Generate the tensor RG flow
```python
python rgzc.py --model hardsquareNeg --scheme loop-TNR --ver rotsym --chi 10 --epspinv 1e-8 --rgn 12
```

*Step 3:* Plot RG flows
```python
python plotRGflow.py --model hardsquareNeg --scheme loop-TNR --ver rotsym --chi 10 --startn 3 --endn 7
```

## III. Apply the usual TRG to the model 
Simply change the scheme and its version in the commands above to
`--scheme trg --version general`.


