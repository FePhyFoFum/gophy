// calcmsbl will calculate the branch lengths for multistate characters
//
package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"time"

	"github.com/FePhyFoFum/gophy"

	"golang.org/x/exp/rand"
)

func main() {
	rand.Seed(uint64(time.Now().UTC().UnixNano()))
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
	wks := flag.Int("w", 4, "number of threads")
	flag.Parse()
	if len(*tfn) == 0 {
		fmt.Fprintln(os.Stderr, "need a tree filename (-t)")
		os.Exit(1)
	}
	if len(*afn) == 0 {
		fmt.Fprintln(os.Stderr, "need a seq filename (-s)")
		os.Exit(1)
	}

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
	fmt.Println(t.Rt.Newick(false))

	//read a seq file
	nsites := 0
	seqs := map[string][]string{}
	mseqs, numstates := gophy.ReadMSeqsFromFile(*afn)
	seqnames := make([]string, 0)
	for _, i := range mseqs {
		seqs[i.NM] = i.SQs
		seqnames = append(seqnames, i.NM)
		nsites = len(i.SQ)
	}
	x := gophy.NewMultStateModel()
	x.NumStates = numstates
	x.SetMap()
	bf := gophy.GetEmpiricalBaseFreqsMS(mseqs, x.NumStates)
	x.SetBaseFreqs(bf)
	x.EBF = x.BF
	// get the site patternas
	patterns, patternsint, gapsites, constant, uninformative, _ := gophy.GetSitePatternsMS(mseqs, x)
	patternval, _ := gophy.PreparePatternVecsMS(t, patternsint, seqs, x)
	//list of sites
	fmt.Fprintln(os.Stderr, "nsites:", nsites)
	fmt.Fprintln(os.Stderr, "patterns:", len(patterns), len(patternsint))
	fmt.Fprintln(os.Stderr, "onlygaps:", len(gapsites))
	fmt.Fprintln(os.Stderr, "constant:", len(constant))
	fmt.Fprintln(os.Stderr, "uninformative:", len(uninformative))

	//start := time.Now()
	// calc likelihood
	w := 10
	if nsites < w {
		w = nsites
	}
	x.SetupQJC()
	l := gophy.PCalcLikePatternsMS(t, x, patternval, *wks)
	fmt.Println("starting lnL:", l)

	//optimize branch lengths
	fmt.Println("start:\n" + t.Rt.Newick(true) + ";")
	gophy.OptimizeBLNRMS(t, x, patternval, 10)
	l = gophy.PCalcLikePatternsMS(t, x, patternval, *wks)
	fmt.Println("ln:", l)
	fmt.Println("end:\n" + t.Rt.Newick(true) + ";")
}
