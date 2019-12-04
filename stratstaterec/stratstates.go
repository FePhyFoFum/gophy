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
	"strconv"
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
	tfn := flag.String("t", "", "file containing a tree calibrated to time")
	afn := flag.String("s", "", "seq filename")
	rfn := flag.String("r", "", "file containing temporal ranges")
	//anc := flag.Bool("a", false, "calc anc states")
	//stn := flag.Bool("n", false, "calc stochastic number (states will also be calculated)")
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
	if len(*rfn) == 0 {
		fmt.Fprintln(os.Stderr, "need a stratigraphic filename (-r)")
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
		ancSeqs := make(map[*gophy.Node][]string)
		for _, n := range t.Pre {
			if len(n.Chs) == 0 {
				ancSeqs[n] = seqs[n.Nam]
			}
		}
		gophy.ReadStrat(*rfn, t)
		gophy.MakeStratHeights(t)
		fmt.Println(t.Rt.NewickFloatBL("TimeLen"))
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
		fmt.Println("--------------------------------")
		fmt.Println("-----estimating anc states------")
		fmt.Println("--------------------------------")
		//patternProbs := ancState(t, x, patternval, sv, fullpattern)
		//fmt.Println(patternProbs)

		rets := gophy.CalcAncStatesMS(x, t, patternval)
		for i := range rets {
			//fmt.Println(rets[i])
			var states []string
			for _, j := range fullpattern {
				//fmt.Println(j, sv.Idx[j])
				actj := sv.Idx[j]
				maxm := 0
				maxv := rets[i][actj][0]
				for m, v := range rets[i][actj] {
					if v > maxv {
						maxm = m
						maxv = v
					}
				}
				states = append(states, strconv.Itoa(maxm))
			}
			ancSeqs[i] = states
		}
		fmt.Println("-------------------------------------")
		fmt.Println("------NEW VS REDERIVED STATES--------")
		fmt.Println("-------------------------------------")
		//nchanges, newstates, nochange := calcNewVsRederived(t, ancSeqs)
		//entropyThruTime(t, ancSeqs, nchanges, newstates, nochange)
		entropyThruTime(t, ancSeqs)
	}
	//close log
	err = lg.Flush()
	if err != nil {
		log.Fatal(err)
	}
}

func calcNewVsRederived(t *gophy.Tree, ancSeqs map[*gophy.Node][]string) (nchanges, newstates, nochange float64) {
	sorted := timeSortNodes(t)
	seenSites := make(map[int]map[string]bool)
	nchanges = 0
	newstates = 0
	nochange = 0
	//fmt.Println("time\tnchanges\tnew_states\tnochange")
	for _, n := range sorted {
		curseq := ancSeqs[n]
		parseq := ancSeqs[n.Par]
		for i, s := range curseq {
			if n == t.Rt {
				if _, ok := seenSites[i]; !ok {
					seenSites[i] = make(map[string]bool)
					seenSites[i][s] = true
					seenSites[i]["-"] = true
					seenSites[i]["?"] = true
				}
				continue
			}
			if s != parseq[i] && s != "-" && s != "?" {
				nchanges++
			} else {
				nochange++
			}
			if _, ok := seenSites[i][s]; !ok {
				seenSites[i][s] = true
				newstates++
			}
		}
		/*if n.Height != 0 {
			fmt.Println(-n.Height, nchanges, newstates, nochange)
		} else {
			fmt.Println(n.Height, nchanges, newstates, nochange)
		}*/
	}
	return
}

/*
func entropyThruTime(t *gophy.Tree, ancSeqs map[*gophy.Node][]string, totalChanges, totalNewStates, totalNoChange float64) {
	pNewState := totalNewStates / (totalChanges + totalNoChange)
	pReversal := (totalChanges - totalNewStates) / (totalChanges + totalNoChange)
	pNoChange := (totalNoChange / (totalChanges + totalNoChange))
	samp := (totalChanges + totalNoChange)
	fmt.Println(pNewState, pReversal, pNoChange, samp)
	sorted := timeSortNodes(t)
	seenSites := make(map[int]map[string]bool)
	nchanges := 0
	newstates := 0
	nochange := 0
	entropy := 0.0
	information := 0.0
	fmt.Println("time\tnchanges\tnew_states\tnochange\tentropy")
	for _, n := range sorted {
		curseq := ancSeqs[n]
		parseq := ancSeqs[n.Par]
		nodenew := 0.0
		nodechanges := 0.0
		nodenochange := 0.0
		for i, s := range curseq {
			if n == t.Rt {
				if _, ok := seenSites[i]; !ok {
					seenSites[i] = make(map[string]bool)
					seenSites[i][s] = true
					seenSites[i]["-"] = true
					seenSites[i]["?"] = true
				}
				continue
			}
			if s != parseq[i] && s != "-" && s != "?" {
				nchanges++
				nodechanges++
				if _, ok := seenSites[i][s]; !ok {
					seenSites[i][s] = true
					newstates++
					nodenew++
					entropy += math.Log2(1.0 / pNewState) // float64(nchanges+nochange))
				} else {
					entropy += math.Log2(1.0 / pReversal) // float64(nchanges+nochange))
				}
			} else {
				nochange++
				nodenochange++
				entropy += math.Log2(1.0 / pNoChange) // float64(nchanges+nochange))
			}
		}
		if nchanges+nochange != 0 {
			if nchanges != 0 {
				information += (nodenew * math.Log2(1.0/(float64(newstates)/float64(nchanges+nochange))))
			}
			reder := nodechanges - nodenew
			fmt.Println(reder, nodenew, newstates, nchanges, nochange, nodenochange, nodechanges, information)
			if reder != 0.0 {
				information += (reder * math.Log2(1.0/(float64(nchanges-newstates)/float64(nchanges+nochange))))
			}
			information += (nodenochange * math.Log2(1.0/(float64(nochange)/(float64(nchanges+nochange)))))
		}
		fmt.Println(information / float64(nchanges+nochange))

		/*
			if n.Height != 0 {
				fmt.Println(-n.Height, nchanges, newstates, nochange, entropy/float64(nchanges+nochange))
			} else {
				fmt.Println(n.Height, nchanges, newstates, nochange, entropy/float64(nchanges+nochange))
			}
	}
}
*/
func entropyThruTime(t *gophy.Tree, ancSeqs map[*gophy.Node][]string) {
	samp := 0.0
	//fmt.Println(pNewState, pReversal, pNoChange, samp)
	sorted := timeSortNodes(t)
	seenSites := make(map[int]map[string]bool)
	nchanges := 0.0
	newstates := 0.0
	nochange := 0.0
	entropy := 0.0
	information := 0.0
	fmt.Println("time\tnchanges\tnew_states\tnochange\tentropy")
	for _, n := range sorted {
		curseq := ancSeqs[n]
		parseq := ancSeqs[n.Par]
		nodenew := 0.0
		nodechanges := 0.0
		nodenochange := 0.0
		for i, s := range curseq {
			if n == t.Rt {
				if _, ok := seenSites[i]; !ok {
					seenSites[i] = make(map[string]bool)
					seenSites[i][s] = true
					seenSites[i]["-"] = true
					seenSites[i]["?"] = true
				}
				continue
			}
			if s != parseq[i] && s != "-" && s != "?" {
				nchanges++
				nodechanges++
				if _, ok := seenSites[i][s]; !ok {
					seenSites[i][s] = true
					newstates++
					nodenew++
				}
			} else {
				nochange++
				nodenochange++
			}
		}
		reder := nodechanges - nodenew
		if nchanges+nochange != 0 {
			var pNew, pReder, pNoChange float64
			pNew = (float64(newstates) / float64(nchanges+nochange))
			if nodenew != 0 {
				information += (nodenew * (pNew * math.Log2(1.0/pNew)))
			}
			pReder = (float64(nchanges-newstates) / float64(nchanges+nochange))
			if reder != 0.0 {
				information += (reder * (pReder * math.Log2(1.0/pReder)))
			}
			pNoChange = (float64(nochange) / float64(nchanges+nochange))
			information += (nodenochange * (pNoChange * math.Log2(1.0/pNoChange)))
			//fmt.Println(pNew, pReder, pNoChange)
		}
		//fmt.Println(reder, nodenew, newstates, nchanges, nochange, nodenochange, nodechanges, information)
		samp += (nodechanges + nodenochange)
		if nchanges+nochange > 0 {
			entropy = information / samp
		}
		//fmt.Println(entropy, information/samp)
		if n.Height != 0 {
			fmt.Println(-n.Height, nchanges, newstates, nochange, entropy)
		} else {
			fmt.Println(n.Height, nchanges, newstates, nochange, entropy)
		}
	}
}

func timeSortNodes(t *gophy.Tree) []*gophy.Node {
	for _, n := range t.Pre {
		if len(n.Chs) == 0 {
			n.Height = n.FData["FAD"]
		}
	}
	sorted := gophy.TimeTraverse(t.Pre, false)
	return sorted
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

func ancState(t *gophy.Tree, x *gophy.MultStateModel, patternval []float64, sv *gophy.SortedIntIdxSlice, fullpattern []int) map[*gophy.Node][][]float64 {
	start := time.Now()
	rets := gophy.CalcAncStatesMS(x, t, patternval)
	ancSeqs := make(map[*gophy.Node][][]float64)
	for i := range rets {
		//fmt.Println(i)
		//fmt.Println(rets[i])
		var patternVecs [][]float64
		for _, j := range fullpattern {
			//fmt.Println(j, sv.Idx[j])
			actj := sv.Idx[j]
			patternVecs = append(patternVecs, rets[i][actj])
			//fmt.Println(" ", j, rets[i][actj])
		}
		ancSeqs[i] = patternVecs
		fmt.Print("\n")
	}
	end := time.Now()
	fmt.Fprintln(os.Stderr, end.Sub(start))
	return ancSeqs
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
