package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/FePhyFoFum/gophy"
)

func printModel(modelparams []float64, basefreqs []float64) {
	fmt.Fprintln(os.Stderr, "basefreqs -- A:", basefreqs[0], " C:", basefreqs[1], " G:", basefreqs[2], " T:", basefreqs[3])
	fmt.Fprintln(os.Stderr, "modelparams --")
	fmt.Fprintln(os.Stderr, " - ", modelparams[0], modelparams[1], modelparams[2])
	fmt.Fprintln(os.Stderr, modelparams[0], " - ", modelparams[3], modelparams[4])
	fmt.Fprintln(os.Stderr, modelparams[1], modelparams[3], " - ", 1.0)
	fmt.Fprintln(os.Stderr, modelparams[2], modelparams[4], 1.0, "-")
}

func main() {
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
	md := flag.String("m", "1.0,1.0,1.0,1.0,1.0", "five params for GTR")
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
	x := gophy.NewDNAModel()
	//x.SetupQJC()
	x.SetNucMap()

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
	patterns, patternsint, gapsites, constant, uninformative, fullpattern := gophy.GetSitePatterns(seqs, nsites, seqnames)
	patternval, patternvec := gophy.PreparePatternVecs(t, patternsint, seqs)
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

	//
	start := time.Now()
	l := gophy.PCalcLikePatterns(t, x, patternval, *wks)
	fmt.Println("starting lnL:", l)
	rets := gophy.CalcAncStates(x, t, patternval)
	charM := gophy.GetRevNucMap()
	for i := range rets {
		fmt.Println(i)
		//fmt.Println(rets[i])
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
			fmt.Print(charM[maxm])
		}
		fmt.Print("\n")
	}
	end := time.Now()
	fmt.Fprintln(os.Stderr, end.Sub(start))
}
