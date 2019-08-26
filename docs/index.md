_gophy_ is not finished and may never be. It serves as a phylogenetic development toolkit from which to make tests and tools for a variety of purposes. 

## installation

Installation of gophy is much like most other github developed _go_ based packages. If you have not setup or installed _go_ before, I recommend you check installation instructions out here. Once you have setup your environment, then you can obtain _gophy_ like any other _go_ package: `go get github.com/FePhyFoFum/gophy`. This will get all of the dependencies with it. To update it if you already have it installed, you can type `go get -u github.com/FePhyFoFum/gophy` or you can update all of your go packages (this might be good at least for _gonum_ that will come with _gophy_) you can type `go get -u all`.

## some useful packages
While gophy is not finished, there are several tools that you may use. Here is a description of how to compile and use those packages.

- [bp](#bp) : bipartition analyzer
- [lentil](#lentil) : 
- [parsbl](#parsbl) : parsimony branch length estimator
- [sites](#sites) : sites toy

### bp
_bp_ is a bipartition analyzer. It does a lot of things and pretty quickly. First, to build it, you run `go build github.com/FePhyFoFum/gophy/bp/bp.go`. That will make an executable called _bp_. You can move that to your PATH or just type the full path to use it. Some features include 
 - trees don't need overlapping taxa sets
 - trees should be unrooted but don't have to be
 - you can use branch length or support cutoffs
 - you can compare sets of trees to a tree or to themselves
 - you can ignore taxa
 - you can run these multicore
 - you can compare only a set of trees in a larger tree file
 - you can plot the concordance, conflict, and uninformativeness (in numbers and percentages) from a set of trees onto another tree

Here are some examples of some runs that you can do with bp. 

_Pairwise comparisons_: You can compare all the trees in a file to each other `bp -t t -pc` or to a pool of biparts made from the trees in the file `bp -t t -pca`. 

_Map concordance and conflict_: _bp_ can take a list of trees, say in a file `t` and map whether bipartitions are concordant, conflicting, or uninformative (based on support or branch length) to another file with a single tree, say in a file `t1`. You run this with `bp -t t -c t1`. This will give a list of the biparts from `t1` and for each that have conflict list the conflicting ones from `t` with the number that conflict (the number in parentheses is just an id for the bipart). If you want to spit out a set of trees with the internal nodes representing conflict, concordance, and uninformative, do `bp -t t -c t1 -tv`. You will get the same output and then a line that says "TREES WITH CONFLICT (FIRST), CONCORDANCE (SECOND), UNSUPPORTED (THIRD), PROPS AFTER" and 6 trees. The first is the number of conflicting trees, concordant trees, and unsupported trees (because of support or branch length cutoffs). The next three are the same order but instead of number, they are proportions.  If you have a big file you can do `bp -t t -c t1 -w 100 -tv` to create a bunch of workers and use multicore. The number of workers is not the number of cores. I generally do ~100. 

### lentil


### parsbl


### sites