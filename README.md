# The adaptive dictionary learning (ADL) toolbox

ADL is a modular dictionary learning toolbox for plug-and-play sparse approximation and dictionary update steps.

## Overview

<b>ADL implements:</b>
* several of the popular sparse approximation techniques
* many of the dictionary update steps
* overall full dictionary learning algorithms (like MOD, K-SVD, ITKrM etc.)

<b>The main features of the toolbox are:</b>
* Adaptive dictionary size
* Training of replacement atoms
* Several image/synthetic demos
* Generation of synthetic datasets
* Coherence constraints
* Structured dictionaries
* Adaptive sparsity
* Avoid fitting noisy atoms
* Several initialization techniques

<b>Who could use the toolbox?</b>

I. dictionary learning researchers:
* easily develop new ideas
* easily compare with previous work
* reproducibility of results

II. other researchers:
* apply dictionary learning to your application
* no/little parameter tuning
* easily compare several dictionary learning techniques and choose the best one

## Examples

### Demo 0

In the first demo [blind run](./demo_0_blind_run.m) we just load some image data and then run the adaptive dictionary learning algorithm without setting any of its parameters. The toolbox falls to its default setting and for each choice it outputs a warning to let the user know about the choice that was made.

### Demos 1 and 2

Then we implement two of the most popular dictionary learning algorithms: [K-SVD](./demo_1_ksvd_algorithm.m) and [ITKrM](./demo_2_itkrm_algorithm.m).

In both cases, we fix the sparsity level and the size of the dictionary. 


### Demo 3

The first two demos are just implementations of well-known dictionary learning algorithm. This demo shows how the toolbox can be used to experiment with new ideas and combine old ones to produce new results. The [hybrid demo](./demo_3_hybrid.m) combines the dictionary update step of K-SVD with the thresholding sparse approximation step used in ITKrM.

The purpose of this demo is to highlight how easy it is to use the toolbox to experiment with different combinations of sub-routines.

### Demo 4

That that we have covered the basic examples, let's move to a more interesting case.

Consider that we want to construct a dictionary with a fixed number of atoms (say K = 80 atoms) such that the representation accuracy is improved over a given number of iterations. The [atom replacement demo](./demo_4_fix_sparsity_replacement_dictionary.m) implements ideas of pruning old atoms and adding new atoms such that the representation error is improved.

A notion of importance (or ranking) is attached to each atom and the worse atoms are pruned while a new training procedure takes care of adding fresh blood to the dictionary.

### Demo 5

The mutual coherence of a dictionary plays an important role in the theoretical understanding of sparse coding in that dictionary. It is therefore useful to have a dictionary learning method that caps the mutual coherence to a maximum imposed value. The [mutual coherence demo](./demo_5_fix_sparsity_replacement_dictionary_max_coherence.m) deals with this problem by pruning highly coherent atoms.

### Demos 6 and 7

In most dictionary learning algorithms there are two parameters that need to be chosen, or tunned: the size of the dictionary and the sparsity level of the representations.

The [adaptive dictionary size demo](./demo_6_fix_sparsity_adapt_dictionary.m) keeps the sparsity level fixed and tries to find the dictionary of best size. Having few atoms we might expect bad results while increasing the number of atoms might leads us to fitting noise components. There is a sweet spot in the middle that we are searching for.

Finally, the [full adaptive demo](./demo_7_adaptive_sparsity_and_dictionary.m) deals with both the sparsity and the dictionary size adaptively.


### References

[1]  Michal Aharon; Michael Elad; Alfred Bruckstein (2006), K-SVD: An Algorithm for Designing Overcomplete Dictionaries for Sparse Representation, IEEE Transactions on Signal Processing, 54 (11): 4311–4322, doi:10.1109/TSP.2006.881199
http://www.cs.technion.ac.il/~elad/publications/journals/2004/32_KSVD_IEEE_TSP.pdf

[2] Karin Schnass (2018), Convergence radius and sample complexity of ITKM algorithms for dictionary learning, Applied and Computational Harmonic Analysis, 45 (1): 22-58<br>
https://arxiv.org/pdf/1503.07027.pdf
