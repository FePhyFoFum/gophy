package main

import (
	"flag"
	"fmt"
	"os"
	"sort"

	"github.com/FePhyFoFum/gophy"
)

func main() {
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
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
	trees := gophy.ReadTreesFromFile(*tfn)
	fmt.Println(len(trees), "trees read")

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
	for _, t := range trees {
		patternval, patternvec := gophy.PreparePatternVecsMS(t, patternsint, seqs, x)
		//this is necessary to get order of the patters in the patternvec since they have no order
		// this will be used with fullpattern to reconstruct the sequences
		sv := gophy.NewSortedIdxSlice(patternvec)
		sort.Sort(sv)
		//fmt.Println(sv.IntSlice, sv.Idx)
		//list of sites
		fmt.Fprintln(os.Stderr, "nsites:", nsites)
		fmt.Fprintln(os.Stderr, "patterns:", len(patterns), len(patternsint))
		fmt.Fprintln(os.Stderr, "onlygaps:", len(gapsites))
		fmt.Fprintln(os.Stderr, "constant:", len(constant))
		fmt.Fprintln(os.Stderr, "uninformative:", len(uninformative))
		// model things
		gophy.PCalcSankParsPatternsMultState(t, x, patternval, 1)
		gophy.EstParsBLMultState(t, x, patternval, nsites)
		fmt.Println(t.Rt.Newick(true) + ";")
	}
}
