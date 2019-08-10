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
	f, err := os.Create("staterec.log")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	lg := bufio.NewWriter(f)
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
	anc := flag.Bool("a", false, "calc anc states")
	stt := flag.Bool("i", false, "calc stochastic time (states will also be calculated)")
	stn := flag.Bool("n", false, "calc stochastic number (states will also be calculated)")
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
	fmt.Fprint(lg, x.BF)
	// get the site patternas
	patterns, patternsint, gapsites, constant, uninformative, fullpattern := gophy.GetSitePatternsMS(mseqs, x.GetCharMap(), x.GetNumStates())

	for _, t := range trees {
		patternval, patternvec := gophy.PreparePatternVecsMS(t, patternsint, seqs, x.GetCharMap(), x.GetNumStates())
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

		optimizeThings(t, x, patternval, *wks)

		//ancestral states
		if *anc || *stt || *stn {
			fmt.Println("--------------------------------")
			fmt.Println("------------anc states----------")
			fmt.Println("--------------------------------")
			ancState(t, x, patternval, sv, fullpattern)
		}

		patternloglikes := make([]float64, len(patternval))
		if *stt || *stn {
			for i := 0; i < len(patternloglikes); i++ {
				patternloglikes[i] = gophy.CalcLogLikeOneSiteMS(t, x, i)
			}
		}
		//stochastic time
		if *stt {
			fmt.Println("--------------------------------")
			fmt.Println("------------stoch time----------")
			fmt.Println("--------------------------------")
			stochTime(t, x, patternval, sv, fullpattern, patternloglikes)
		}

		//stochatic number
		if *stn {
			fmt.Println("--------------------------------")
			fmt.Println("----------stoch number----------")
			fmt.Println("--------------------------------")
			stochNumber(t, x, patternval, sv, fullpattern, patternloglikes)
		}
	}
	//close log
	err = lg.Flush()
	if err != nil {
		log.Fatal(err)
	}
}

func optimizeThings(t *gophy.Tree, x *gophy.MultStateModel, patternval []float64, wks int) {
	x.SetupQJC()
	l := gophy.PCalcLikePatternsMS(t, x, patternval, wks)
	fmt.Println("starting lnL:", l)
	gophy.OptimizeMS1R(t, x, patternval, wks)
	gophy.OptimizeMKMS(t, x, x.Q.At(0, 1), patternval, false, wks)
	gophy.PrintMatrix(x.Q, false)
	l = gophy.PCalcLikePatternsMS(t, x, patternval, wks)
	fmt.Println("optimized lnL:", l)
}

func ancState(t *gophy.Tree, x *gophy.MultStateModel, patternval []float64, sv *gophy.SortedIntIdxSlice, fullpattern []int) {
	start := time.Now()
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
}

func stochTime(t *gophy.Tree, x *gophy.MultStateModel, patternval []float64, sv *gophy.SortedIntIdxSlice, fullpattern []int, patternloglikes []float64) {

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
}

func stochNumber(t *gophy.Tree, x *gophy.MultStateModel, patternval []float64, sv *gophy.SortedIntIdxSlice, fullpattern []int, patternloglikes []float64) {
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
