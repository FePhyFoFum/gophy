package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"runtime/pprof"
	"sort"
	"time"

	"github.com/FePhyFoFum/gophy"
	"gonum.org/v1/gonum/mat"
)

func main() {
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
	//md := flag.Bool("m", false, "model params free")
	//ebf := flag.Bool("b", true, "use empirical base freqs (alt is estimate)")
	wks := flag.Int("w", 4, "number of threads")
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")
	flag.Parse()
	if len(*tfn) == 0 {
		fmt.Fprintln(os.Stderr, "need a tree filename (-t)")
		os.Exit(1)
	}
	if len(*afn) == 0 {
		fmt.Fprintln(os.Stderr, "need a seq filename (-s)")
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
	x := gophy.NewMULTModel()
	x.NumStates = 4
	x.SetMap()
	fmt.Println(x.CharMap)

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
	seqs := map[string][]string{}
	mseqs := gophy.ReadMSeqsFromFile(*afn)
	seqnames := make([]string, 0)
	for _, i := range mseqs {
		seqs[i.NM] = i.SQs
		seqnames = append(seqnames, i.NM)
		nsites = len(i.SQ)
	}
	bf := gophy.GetEmpiricalBaseFreqsMS(mseqs, x.NumStates)
	x.SetBaseFreqs(bf)
	fmt.Println(x.BF)
	// get the site patternas
	patterns, patternsint, gapsites, constant, uninformative, fullpattern := gophy.GetSitePatternsMS(mseqs, x)
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
	x.SetupQJC()
	//x.SetupQGTR()
	fmt.Println(x.Q)
	//os.Exit(0)

	start := time.Now()
	l := gophy.PCalcLikePatternsMS(t, x, patternval, *wks)
	fmt.Println("starting lnL:", l)
	gophy.OptimizeMULT1R(t, x, patternval, *wks)
	l = gophy.PCalcLikePatternsMS(t, x, patternval, *wks)
	fmt.Println("optimized lnL:", l)

	//ancestral states
	fmt.Println("--------------------------------")
	fmt.Println("------------anc states----------")
	fmt.Println("--------------------------------")
	rets := gophy.CalcAncStatesMS(x, t, patternval)
	for i := range rets {
		fmt.Println(i)
		//fmt.Println(rets[i])
		for _, j := range fullpattern {
			//fmt.Println(j, sv.Idx[j])
			actj := sv.Idx[j]
			fmt.Println(" ", j, rets[i][actj])
		}
		fmt.Print("\n")
	}
	end := time.Now()
	fmt.Fprintln(os.Stderr, end.Sub(start))

	//stochastic time
	fmt.Println("--------------------------------")
	fmt.Println("------------stoch time----------")
	fmt.Println("--------------------------------")
	patternloglikes := make([]float64, len(patternval))
	for i := 0; i < len(patternloglikes); i++ {
		patternloglikes[i] = gophy.CalcLogLikeOneSiteMS(t, x, i)
	}
	sttimes := make(map[*gophy.Node][][]float64)
	for _, nd := range t.Post {
		if nd != t.Rt {
			sttimes[nd] = make([][]float64, len(patternval))
			for j := range patternval {
				sttimes[nd][j] = make([]float64, x.NumStates)
			}
		}
	}
	for st := 0; st < x.NumStates; st++ {
		retsS := gophy.CalcStochMapMS(x, t, patternval, true, st, st)
		for i := range retsS {
			if i == t.Rt {
				continue
			}
			for j := range retsS[i] {
				sl := patternloglikes[j]
				s := 0.
				for _, m := range retsS[i][j] {
					s += m
				}
				sttimes[i][j][st] = s / math.Exp(sl)
			}
		}
	}
	for nd := range sttimes {
		if nd == t.Rt {
			continue
		}
		fmt.Println(nd)
		for _, j := range fullpattern {
			actj := sv.Idx[j]
			fmt.Println("  ", j, " ", sttimes[nd][actj])
		}
	}
	//stochatic number
	/*
	 *
	 */
	fmt.Println("--------------------------------")
	fmt.Println("----------stoch number----------")
	fmt.Println("--------------------------------")
	stnum := make(map[*gophy.Node][]*mat.Dense)
	for _, nd := range t.Post {
		if nd != t.Rt {
			stnum[nd] = make([]*mat.Dense, len(patternval))
			for j := range patternval {
				stnum[nd][j] = mat.NewDense(x.NumStates, x.NumStates, nil)
				stnum[nd][j].Zero()
			}
		}
	}
	for st := 0; st < x.NumStates; st++ {
		for st2 := 0; st2 < x.NumStates; st2++ {
			if st == st2 {
				continue
			}
			retsS := gophy.CalcStochMapMS(x, t, patternval, false, st2, st) //this must be reversed, don't get confused
			for i := range retsS {
				if i == t.Rt {
					continue
				}
				for j := range retsS[i] {
					sl := patternloglikes[j]
					s := 0.
					for _, m := range retsS[i][j] {
						s += m
					}
					stnum[i][j].Set(st, st2, s/math.Exp(sl))
				}
			}
		}
	}
	for nd := range stnum {
		if nd == t.Rt {
			continue
		}
		fmt.Println(nd)
		for _, j := range fullpattern {
			actj := sv.Idx[j]
			gophy.PrintMatrix(stnum[nd][actj], true)
		}
	}
}
