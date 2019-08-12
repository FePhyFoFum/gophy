package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"time"

	"github.com/carolinetomo/gophy"
)

func main() {
	treeArg := flag.String("t", "", "input tree")
	traitArg := flag.String("m", "", "continuous traits")
	//mclArg := flag.String("start", "", "user-specified starting clusters")
	genArg := flag.Int("gen", 500000, "number of MCMC generations to run")
	kArg := flag.Int("K", 2, "maximum number of clusters")
	minKArg := flag.Int("minK", 1, "minimum number of clusters")
	printFreqArg := flag.Int("pr", 100, "Frequency with which to print to the screen")
	searchArg := flag.Int("f", 3, "0\tOptimize branch lengths for a user-specified clustering\n1\tOutput distance matrices calculated for each cluster provided by the -start argument\n2\tCalculate the log-likelihood of the dataset on a particular topology\n3\tPerform cluster analysis")
	//sampFreqArg := flag.Int("samp", 1, "Frequency with which to sample from the chain")
	runNameArg := flag.String("o", "gophy", "specify the prefix for outfile names")
	critArg := flag.Int("c", 0, "Criterion to use for hill climbing:\n0\tAIC\n1\tBIC\n2\tAICc\n")
	splitGenArg := flag.Int("split", 10, "Number of iterations to run at each splitting step")
	profileArg := flag.Bool("prof", false, "indicate whether to run the go profiler (for development)")
	//threadArg := flag.Int("T", 1, "maximum number of cores to use during run")
	//workersArg := flag.Int("W", 4, "Number of Go workers to use for LL calculation concurrency")
	clustArg := flag.Float64("a", 1.0, "concentration parameter for new cluster penalty")
	flag.Parse()
	if *profileArg == true {
		f, err := os.Create("profile.prof")
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}
	//var ntax,ntraits int
	runtime.GOMAXPROCS(2)
	nwk := gophy.ReadLine(*treeArg)[0]
	rt := gophy.ReadNewickString(nwk)
	t := gophy.NewTree()
	t.Instantiate(rt)

	gophy.MapContinuous(t, *traitArg)
	rand.Seed(time.Now().UTC().UnixNano())
	for _, n := range t.Pre[1:] {
		r := rand.Float64()
		n.BMLen = r
	}
	gophy.InitMissingValues(t.Pre)
	//gophy.BMOptimBLEM(t, 2) //going to optimize branch lengths to set mean parameter for tree length in dirichlet prior
	treeOutFile := *runNameArg
	if *searchArg == 3 {
		search := gophy.InitGreedyHC(t, *genArg, *printFreqArg, *critArg, true, *kArg, treeOutFile, *splitGenArg, *clustArg, *minKArg)
		//fmt.Println(search.ClusterString())
		start := time.Now()
		search.PerturbedRun()
		elapsed := time.Since(start)
		fmt.Println("COMPLETED ", *genArg, "ITERATIONS IN ", elapsed)
	}
}
