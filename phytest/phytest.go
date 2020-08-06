// phytest is just for testing different functions in the phylogenetics stuff
//
package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"runtime/pprof"
	"strconv"
	"strings"
	"time"

	"github.com/FePhyFoFum/gophy"

	"golang.org/x/exp/rand"
)

func slidingWindow(x, w float64) float64 {
	return math.Abs((rand.Float64() * ((x + w/2) - (x - w/2))) + (x - w/2.))
}

// MCMC simple MCMC for branch lengths
func MCMC(t *gophy.Tree, x *gophy.DNAModel, patternval []float64, wks int, outfilename string) {
	f, err := os.Create(outfilename)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	wr := bufio.NewWriter(f)
	w := 0.1
	printf := 1000
	iter := 100000
	curlike := gophy.PCalcLogLikePatterns(t, x, patternval, wks)
	newlike := 0.0
	start := time.Now()
	//fullcalc := false
	for i := 0; i < iter; i++ {
		nd := rand.Intn(len(t.Post))
		for t.Post[nd] == t.Rt {
			nd = rand.Intn(len(t.Post))
		}
		odv := t.Post[nd].Len
		newv := slidingWindow(odv, w)
		t.Post[nd].Len = newv
		t.Post[nd].Marked = true
		if len(t.Post[nd].Chs) == 0 {
			t.Post[nd].Par.Marked = true
		}
		//if fullcalc {
		//newlike = PCalcLogLike(t, x, nsites, wks)
		//} else {
		//	newlike = PCalcLogLikeBack(t, t.Post[nd], x, nsites, wks)
		//}
		//newlike = gophy.PCalcLogLikeMarked(t, x, nsites, wks)
		//make pattern marked
		newlike = gophy.PCalcLikePatternsMarked(t, x, patternval, wks)
		for _, n := range t.Post {
			n.Marked = false
		}
		r := math.Log(rand.Float64())
		if r < newlike-curlike {
			//if newlike-curlike > 0 {
			curlike = newlike
			//fullcalc = false
		} else {
			t.Post[nd].Len = odv
			t.Post[nd].Marked = true
			if len(t.Post[nd].Chs) == 0 {
				t.Post[nd].Par.Marked = true
			}
			//fullcalc = true
		}
		if i%printf == 0 {
			curt := time.Now()
			fmt.Fprintln(os.Stderr, i, curlike, "(last", printf, ":", curt.Sub(start), ")")
			fmt.Fprint(wr, t.Rt.Newick(true)+";\n")
			start = time.Now()
		}
	}
	end := time.Now()
	fmt.Println(end.Sub(start))
	err = wr.Flush()
	if err != nil {
		log.Fatal(err)
	}
}

func printModel(modelparams []float64, basefreqs []float64) {
	fmt.Fprintln(os.Stderr, "basefreqs -- A:", basefreqs[0], " C:", basefreqs[1], " G:", basefreqs[2], " T:", basefreqs[3])
	fmt.Fprintln(os.Stderr, "modelparams --")
	fmt.Fprintln(os.Stderr, " - ", modelparams[0], modelparams[1], modelparams[2])
	fmt.Fprintln(os.Stderr, modelparams[0], " - ", modelparams[3], modelparams[4])
	fmt.Fprintln(os.Stderr, modelparams[1], modelparams[3], " - ", 1.0)
	fmt.Fprintln(os.Stderr, modelparams[2], modelparams[4], 1.0, "-")
}

func main() {
	rand.Seed(uint64(time.Now().UTC().UnixNano()))
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
	md := flag.String("m", "1.0,1.0,1.0,1.0,1.0", "five params for GTR")
	wks := flag.Int("w", 4, "number of threads")
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")
	flag.Parse()
	if len(*tfn) == 0 {
		flag.PrintDefaults()
		os.Exit(1)
	}
	if len(*afn) == 0 {
		flag.PrintDefaults()
		os.Exit(1)
	}
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}
	x := gophy.NewDNAModel()
	//x.SetupQJC()
	x.SetMap()

	//read a tree file
	f, err := os.Open(*tfn)
	if err != nil {
		fmt.Println(err)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	var rt *gophy.Node
	t := gophy.NewTree()
	for scanner.Scan() {
		ln := scanner.Text()
		if len(ln) < 2 {
			continue
		}
		rt = gophy.ReadNewickString(ln)
		t.Instantiate(rt)
	}

	fmt.Println(t.Rt.Newick(true))

	//read a seq file
	nsites := 0
	seqs := map[string]string{}
	seqnames := make([]string, 0)
	for _, i := range gophy.ReadSeqsFromFile(*afn) {
		seqs[i.NM] = i.SQ
		seqnames = append(seqnames, i.NM)
		nsites = len(i.SQ)
	}
	bf := gophy.GetEmpiricalBaseFreqs(seqs)
	x.SetBaseFreqs(bf)
	// get the site patternas
	patterns, patternsint, gapsites, constant, uninformative, _ := gophy.GetSitePatterns(seqs, nsites, seqnames)
	patternval, _ := gophy.PreparePatternVecs(t, patternsint, seqs)
	//list of sites
	fmt.Fprintln(os.Stderr, "nsites:", nsites)
	fmt.Fprintln(os.Stderr, "patterns:", len(patterns), len(patternsint))
	fmt.Fprintln(os.Stderr, "onlygaps:", len(gapsites))
	fmt.Fprintln(os.Stderr, "constant:", len(constant))
	fmt.Fprintln(os.Stderr, "uninformative:", len(uninformative))
	// model things
	mds := strings.Split(*md, ",")
	modelparams := make([]float64, 5)
	if len(mds) != 5 {
		fmt.Fprintln(os.Stderr, "your model contains ", len(mds), " params, not 5")
		os.Exit(1)
	} else {
		for i, j := range mds {
			f, err := strconv.ParseFloat(j, 64)
			if err != nil {
				fmt.Fprintln(os.Stderr, "problem parsing ", j, " as float in model specs")
				os.Exit(1)
			}
			modelparams[i] = f
		}
	}
	printModel(modelparams, bf)
	x.SetRateMatrix(modelparams)
	x.SetupQGTR()
	fmt.Println(x.Q)
	//
	l := gophy.PCalcLikePatterns(t, x, patternval, *wks)
	fmt.Println("starting lnL:", l)
	fmt.Println("getting starting branch lengths (parsimony)")
	s := gophy.PCalcSankParsPatterns(t, patternval, *wks)
	fmt.Println("sank:", s)
	gophy.EstParsBL(t, patternval, nsites)
	fmt.Println("start:\n" + t.Rt.Newick(true) + ";")
	//os.Exit(0)

	// calc likelihood
	start := time.Now()
	l = gophy.PCalcLikePatterns(t, x, patternval, *wks)
	fmt.Println("lnL:", l)

	/*
		//optimize model
		gophy.OptimizeGTR(t, x, patternval, 10)
		//MCMC(t, x, patternval, 3, "temp.mcmc.tre")
		//print the matrix
		for i := 0; i < 4; i++ {
			for j := 0; j < 4; j++ {
				fmt.Print(x.R.At(i, j), " ")
			}
			fmt.Print("\n")
		}
		//print the basefreqs
		fmt.Println(x.BF)
		fmt.Println(t.Rt.Newick(true))

		//fmt.Println(gophy.GetGammaCats(10, 5, false))
		l = gophy.PCalcLikePatterns(t, x, patternval, *wks)
		fmt.Println("lnL:", l)
	*/

	gophy.OptimizeBLNR(t, x, patternval, 10)
	l = gophy.PCalcLikePatterns(t, x, patternval, *wks)
	fmt.Println("final:\n" + t.Rt.Newick(true) + ";")
	fmt.Println("lnL:", l)
	end := time.Now()
	fmt.Fprintln(os.Stderr, end.Sub(start))
}
