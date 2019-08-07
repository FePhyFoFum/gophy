package main

import (
	"bufio"
	"flag"
	"fmt"
	"gophy"
	"log"
	"math"
	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"strconv"
	"strings"
	"time"
)

func main() {
	treeArg := flag.String("t", "", "input tree")
	traitArg := flag.String("m", "", "continuous traits")
	mclArg := flag.String("start", "", "user-specified starting clusters")
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
	tree := gophy.ReadTree(nwk)
	traits, _, ntraits := gophy.ReadContinuous(*traitArg)
	gophy.MapContinuous(tree, traits, ntraits)
	rand.Seed(time.Now().UTC().UnixNano())
	if *searchArg == 2 {
		gophy.InitMissingValues(tree.PreorderArray())
		gophy.MissingTraitsEM(tree, 100)
		LL := 0.0
		for site := range tree.ContData {
			LL += gophy.SingleSiteLL(tree, site)
		}
		fmt.Println(LL)
		fmt.Println(tree.Newick(true) + ";")
		os.Exit(0)
	}
	for _, n := range tree.PreorderArray()[1:] {
		r := rand.Float64()
		n.Len = r
	}
	gophy.InitMissingValues(tree.PreorderArray())
	gophy.MissingTraitsEM(tree, 100) //going to optimize branch lengths to set mean parameter for tree length in dirichlet prior
	//fmt.Println(tree.Newick(true))
	gophy.InitParallelPRNLen(tree.PreorderArray())
	//fmt.Println("Starting tree AIC/BIC:", gophy.CalcTreeAIC(tree, *critArg))
	//fmt.Println(tree.Newick(true))
	treeOutFile := *runNameArg
	if *searchArg == 3 {
		search := gophy.InitGreedyHC(tree, *genArg, *printFreqArg, *critArg, true, *kArg, treeOutFile, *splitGenArg, *clustArg, *minKArg)
		//fmt.Println(search.ClusterString())
		start := time.Now()
		search.PerturbedRun()
		elapsed := time.Since(start)
		fmt.Println("COMPLETED ", *genArg, "ITERATIONS IN ", elapsed)
	} else if *searchArg == 0 {
		if *mclArg == "" {
			fmt.Println("You need to specify a cluster input file to run this option")
			os.Exit(1)
		}
		clusters := gophy.ReadMCLoutput(*mclArg)
		nodes := tree.PreorderArray()
		clmap, err := os.Create("cluster_key_dist")
		if err != nil {
			log.Fatal(err)
		}
		w1 := bufio.NewWriter(clmap)
		for lab, c := range clusters {
			f, err := os.Create(strconv.Itoa(lab) + ".bl.tre")
			if err != nil {
				log.Fatal(err)
			}
			w := bufio.NewWriter(f)

			for _, n := range nodes[1:] {
				r := rand.Float64()
				n.Len = r
			}
			gophy.ClusterMissingTraitsEM(tree, c, 100)
			//sites := ""
			//for _, site := range c.Sites {
			//	sites += strconv.Itoa(site) + "\t"
			//}
			fmt.Fprint(w, tree.Newick(true)+";")
			err = w.Flush()
			if err != nil {
				log.Fatal(err)
			}
			f.Close()
			var strsites []string
			for _, site := range c.Sites {
				strsites = append(strsites, strconv.Itoa(site))
			}
			fmt.Fprint(w1, "CLUSTER"+strconv.Itoa(lab)+"\t"+strings.Join(strsites, "\t")+"\n")
		}
		err = w1.Flush()
		if err != nil {
			log.Fatal(err)
		}
		clmap.Close()
	} else if *searchArg == 1 {
		if *mclArg == "" {
			fmt.Println("You need to specify a cluster input file to run this option")
			os.Exit(1)
		}

		clusters := gophy.ReadMCLoutput(*mclArg)
		nodes := tree.PreorderArray()
		clmap, err := os.Create("cluster_key_dist")
		if err != nil {
			log.Fatal(err)
		}
		w1 := bufio.NewWriter(clmap)
		for lab, c := range clusters {
			f, err := os.Create(strconv.Itoa(lab) + ".dist.phy")
			if err != nil {
				log.Fatal(err)
			}
			f1, err := os.Create(strconv.Itoa(lab) + ".phy")
			if err != nil {
				log.Fatal(err)
			}
			w1 := bufio.NewWriter(f1)
			w := bufio.NewWriter(f)
			dm := gophy.SubDM(nodes, c)
			out := gophy.DMtoPhylip(dm, nodes)
			fmt.Fprint(w, strings.Join(out, "\n"))
			fmt.Fprint(w1, c.WriteClusterPhylip(nodes))
			err = w1.Flush()
			if err != nil {
				log.Fatal(err)
			}
			f1.Close()
			err = w.Flush()
			if err != nil {
				log.Fatal(err)
			}
			f.Close()
			var strsites []string
			for _, site := range c.Sites {
				strsites = append(strsites, strconv.Itoa(site))
			}
			fmt.Fprint(w1, "CLUSTER"+strconv.Itoa(lab)+"\t"+strings.Join(strsites, "\t")+"\n")
		}
		err = w1.Flush()
		if err != nil {
			log.Fatal(err)
		}
		clmap.Close()
	} else if *searchArg == 4 {
		if *mclArg == "" {
			fmt.Println("You need to specify a cluster input file to run this option")
			os.Exit(1)
		}
		clusters := gophy.ReadMCLoutput(*mclArg)
		nodes := tree.PreorderArray()
		clustSiteLikes := make(map[int][]float64)
		f, err := os.Create(*runNameArg + ".tab")
		if err != nil {
			log.Fatal(err)
		}
		w := bufio.NewWriter(f)
		for lab, c := range clusters {
			for _, n := range nodes[1:] {
				r := rand.Float64()
				n.Len = r
			}
			gophy.ClusterMissingTraitsEM(tree, c, 10)
			sitelikes := gophy.SitewiseLogLike(tree)
			clustSiteLikes[lab] = sitelikes
		}
		for lab, c := range clusters {
			for _, site := range c.Sites {
				assignLL := clustSiteLikes[lab][site]
				totalDens := math.Exp(assignLL)
				for kk, likes := range clustSiteLikes {
					if kk != lab {
						curlike := likes[site]
						totalDens += math.Exp(curlike)
					}
				}
				fmt.Println(strconv.FormatFloat(math.Exp(assignLL)/totalDens, 'f', -1, 64))
				fmt.Fprintln(w, strconv.FormatFloat(math.Exp(assignLL)/totalDens, 'f', -1, 64))
			}
		}
		err = w.Flush()
		if err != nil {
			log.Fatal(err)
		}
		f.Close()
	}
}
