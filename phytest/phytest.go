package main

import (
	"bufio"
	"flag"
	"fmt"
	"gophy"
	"log"
	"math"
	"os"
	"runtime/pprof"
	"strconv"
	"strings"
	"time"

	"golang.org/x/exp/rand"
)

func slidingWindow(x, w float64) float64 {
	return (rand.Float64() * ((x + w/2) - (x - w/2))) + (x - w/2.)
}

// MCMC simple MCMC for branch lengths
func MCMC(t *gophy.Tree, x *gophy.DNAModel, nsites int, wks int, outfilename string) {
	f, err := os.Create(outfilename)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	wr := bufio.NewWriter(f)
	w := 0.5
	printf := 100
	iter := 10000
	curlike := gophy.PCalcLogLike(t, x, nsites, wks)
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
		newlike = gophy.PCalcLogLikeMarked(t, x, nsites, wks)
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

func getPatterns(seqs map[string]string, nsites int, seqnames []string) (patterns map[string][]int,
	patternsint map[int]float64, gapsites []int, constant []int, uninformative []int) {
	patterns = make(map[string][]int)
	for k := 0; k < nsites; k++ {
		tp := ""
		As, Cs, Gs, Ts, gapcount := 0, 0, 0, 0, 0
		for _, j := range seqnames {
			tp += string(seqs[j][k])
			switch c := string(seqs[j][k]); c {
			case "A":
				As++
			case "C":
				Cs++
			case "G":
				Gs++
			case "T":
				Ts++
			case "-":
				gapcount++
			case "N":
				gapcount++
			default:
				//fmt.Println(c)
			}
		}
		efc := len(seqs) - gapcount
		if gapcount == len(seqs) {
			gapsites = append(gapsites, k)
			continue
		} else if As >= efc || Cs >= efc || Gs >= efc || Ts >= efc {
			constant = append(constant, k)
		}
		twocount := 0
		if As >= 2 {
			twocount++
		}
		if Cs >= 2 {
			twocount++
		}
		if Gs >= 2 {
			twocount++
		}
		if Ts >= 2 {
			twocount++
		}
		if twocount < 2 {
			uninformative = append(uninformative, k)
		}
		if _, ok := patterns[tp]; !ok {
			patterns[tp] = make([]int, 0)
		}
		patterns[tp] = append(patterns[tp], k)
	}
	patternsint = make(map[int]float64) // key is first site, value is the number of that one
	for _, j := range patterns {
		patternsint[j[0]] = float64(len(j))
	}
	return
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
		fmt.Fprintln(os.Stderr, "need a filename")
		os.Exit(1)
	}
	if len(*afn) == 0 {
		fmt.Fprintln(os.Stderr, "need a filename")
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
	patterns, patternsint, gapsites, constant, uninformative := getPatterns(seqs, nsites, seqnames)
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
	//each site has a vector of 4 for DNA conditionals
	/*for _, n := range t.Post {
		n.Data = make([][]float64, nsites)
		for i := 0; i < nsites; i++ {
			n.Data[i] = []float64{0.0, 0.0, 0.0, 0.0}
		}
		if len(n.Chs) == 0 {
			for i := 0; i < nsites; i++ {
				for j := range x.CharMap[string(seqs[n.Nam][i])] {
					n.Data[i][j] = 1.0
				}
			}
		}
	}*/
	//TRYING PATTERNS

	patternvec := make([]int, len(patternsint))     //which site
	patternval := make([]float64, len(patternsint)) //log of number of sites
	count := 0
	for i := range patternsint {
		patternvec[count] = i
		patternval[count] = patternsint[i]
		count++
	}
	for _, n := range t.Post {
		n.Data = make([][]float64, len(patternsint))
		for i := 0; i < len(patternsint); i++ {
			n.Data[i] = []float64{0.0, 0.0, 0.0, 0.0}
		}
		if len(n.Chs) == 0 {
			count := 0
			for _, i := range patternvec {
				if _, ok := x.CharMap[string(seqs[n.Nam][i])]; !ok {
					if string(seqs[n.Nam][i]) != "-" && string(seqs[n.Nam][i]) != "N" {
						fmt.Println(string(seqs[n.Nam][i]))
						os.Exit(0)
					}
				}
				for j := range x.CharMap[string(seqs[n.Nam][i])] {
					n.Data[count][j] = 1.0
				}
				count++
			}
		}
	}

	// calc likelihood
	start := time.Now()
	w := 10
	if nsites < w {
		w = nsites
	}
	l := gophy.PCalcLikePatterns(t, x, patternval, *wks)
	fmt.Println("lnL:", l)
	//fmt.Println(t.Rt.Newick(true))
	/*
		for _, n := range t.Post {
			OptimizeBL(n, t, x, nsites, 10)
		}
		for _, n := range t.Post {
			OptimizeBL(n, t, x, nsites, 10)
		}
		for _, n := range t.Post {
			OptimizeBL(n, t, x, nsites, 10)
		}
		for _, n := range t.Post {
			OptimizeBL(n, t, x, nsites, 10)
		}
		for _, n := range t.Post {
			OptimizeBL(n, t, x, nsites, 10)
		}*/
	//MCMC(t, x, nsites, 10, "temp.mcmc.tre")
	end := time.Now()
	fmt.Fprintln(os.Stderr, end.Sub(start))
}
