# cda-samplers
Examples of fitting a simple Gaussian mixture model with three different Bayesian sampling methods.

Developed for a lecture for Carnegie's Computational Data Analysis program

## Setup
Install these packages.
* Install `dynesty`: `pip install dynesty` https://dynesty.readthedocs.io/en/latest/
  * Install `tqdm` as part of this: `pip install tqdm`
* Install `emcee`: `pip install emcee` https://emcee.readthedocs.io/en/stable/user/install/
* Install `pystan`: `pip install pystan` https://pystan.readthedocs.io/en/latest/getting_started.html
* Install `corner`: `pip install corner` https://corner.readthedocs.io/en/latest/
* Install `schwimmbad`: `pip install schwimmbad` https://schwimmbad.readthedocs.io/en/latest/

2024 update: 
* `conda create -n cda-samplers-2024-05 python jupyter numpy scipy emcee dynesty corner schwimmbad astropy -c conda-forge`


When installing `pystan`, it involves compiling the code Stan (C++ code).
Whatever compiler is used here needs to be in your path when running Stan.
This is the alias Alex needs to call before running the ipython notebook for his version to work:
```
alias stancc="export CC=<your-version-of-clang>; export CXX=<your-version-of-clang++>"
```

If you have not installed standard astronomy libraries (numpy, scipy, matplotlib, astropy) you should do that too.

## Data
The dataset we will use is not yet published but will be soon.
For now, Alex will send you the file.

# Reference

## Alex's recommendation
This is motivated by empirical practice with relatively small problems.
* Use `dynesty` by default. It's slowest but most robust. It's easier to spend computer time than human time.
* Use `emcee` if `dynesty` is too slow and/or you have a complicated model.
* Use `pystan` if `emcee` is too slow and your model can be written in terms of analytic functions

## `dynesty`
Dynamic nested sampling.

Benefits:
* Scientists are better at thinking about what prior they want than the intricacies of sampling
* Good at multimodal stuff
* Computes the evidence (I think that means it's good for model comparison)

Costs:
* Slow and costly, about 10x more than other methods.
* Need to be very careful about choosing your priors

## `emcee`
Affine-invariant multiple walkers

Benefits:
* Probably easiest to use
* Lots of people in astronomy use it and can help you

Costs:
* Need some practice with setting up walkers, how long to burn in, etc.
* It's a bit fiddly how long you need to run
* Not the most efficient sampler

## `pystan`
Hamiltonian Monte Carlo in Stan

Benefits:
* By far the fastest and most efficient sampler discussed here
* Actively developed by statisticians to be good, "bleeding edge"
* A good forum and documentation https://mc-stan.org/users/documentation/

Costs:
* The hardest to learn: you kind of need to know C++, it's a whole new language on top of that, and the language is only "intuitive" for statisticians
* Cannot use very complicated models (e.g. interpolating grids)
* The sampler is so efficient that you have to be a bit careful about priors too
