# gophy

gophy is not finished and may never be. It serves as a phylogenetic development toolkit from which to make tests and tools for a variety of purposes. 

## installation

Installation of gophy is much like most other github developed _go_ based packages. If you have not setup or installed _go_ before, I recommend you check installation instructions out here. Once you have setup your environment, then you can obtain _gophy_ like any other _go_ package: `go get github.com/FePhyFoFum/gophy`. This will get all of the dependencies with it. To update it if you already have it installed, you can type `go get -u github.com/FePhyFoFum/gophy` or you can update all of your go packages (this might be good at least for _gonum_ that will come with _gophy_) you can type `go get -u all`.

## some useful packages
While gophy is not finished, there are several tools that you may use. Here is a description of how to compile and use those packages.

### bp
_bp_ is a bipartition analyzer. It does a lot of things and pretty quickly. First, to build it, you run `go build github.com/FePhyFoFum/gophy/bp/bp.go`. That will make an executable called _bp_. You can move that to your PATH or just type the full path to use it. Here are some examples of some runs that you can do with bp.