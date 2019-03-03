package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"runtime/pprof"
	"strconv"
	"strings"
	"time"

	"github.com/FePhyFoFum/gophy"

	"golang.org/x/exp/rand"
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
	rand.Seed(uint64(time.Now().UTC().UnixNano()))
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
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

	maxint := 0
	nodemodels := map[*gophy.Node]int{}
	//get number of models from tree labels
	for _, i := range t.Pre {
		if len(i.Nam) == 0 {
			nodemodels[i] = 0
		} else if len(i.Chs) == 0 { //if parent has the num and it is tip, then it gets the num
			if strings.Contains(i.Nam, "#") {
				spls := strings.Split(i.Nam, "#")
				i.Nam = spls[0]
				num64, _ := strconv.ParseInt(spls[1], 0, 0)
				nodemodels[i] = int(num64)
				fmt.Println(i, num64)
			}
			nodemodels[i] = nodemodels[i.Par]

		} else if i.Nam[0] == '#' {
			num64, _ := strconv.ParseInt(i.Nam[1:], 0, 0)
			num := int(num64)
			if num > maxint {
				maxint = num
			}
			fmt.Println(i, num)
			nodemodels[i] = int(num)
		} else {
			nodemodels[i] = 0
		}
	}
	fmt.Println("number of models:", maxint+1)
	//read a seq file
	nsites := 0
	seqs := map[string]string{}
	seqnames := make([]string, 0)
	for _, i := range gophy.ReadSeqsFromFile(*afn) {
		seqs[i.NM] = i.SQ
		seqnames = append(seqnames, i.NM)
		nsites = len(i.SQ)
	}
	models := make([]*gophy.DNAModel, maxint+1)
	for i := 0; i <= maxint; i++ {
		// start model
		x := gophy.NewDNAModel()
		x.SetNucMap()
		//empirical freqs for just the relevant seqs
		tseqs := map[string]string{}
		for tn := range nodemodels {
			if len(tn.Chs) == 0 && nodemodels[tn] == i {
				tseqs[tn.Nam] = seqs[tn.Nam]
			}
		}
		bf := gophy.GetEmpiricalBaseFreqs(tseqs)
		x.SetBaseFreqs(bf)
		// model things
		modelparams := []float64{1.0, 1.0, 1.0, 1.0, 1.0}
		printModel(modelparams, bf)
		x.SetRateMatrix(modelparams)
		x.SetupQGTR()
		models[i] = x
	}
	//end

	// get the site patternas
	patterns, patternsint, gapsites, constant, uninformative := gophy.GetSitePatterns(seqs, nsites, seqnames)
	//list of sites
	fmt.Fprintln(os.Stderr, "nsites:", nsites)
	fmt.Fprintln(os.Stderr, "patterns:", len(patterns), len(patternsint))
	fmt.Fprintln(os.Stderr, "onlygaps:", len(gapsites))
	fmt.Fprintln(os.Stderr, "constant:", len(constant))
	fmt.Fprintln(os.Stderr, "uninformative:", len(uninformative))

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
				if _, ok := models[0].CharMap[string(seqs[n.Nam][i])]; !ok {
					if string(seqs[n.Nam][i]) != "-" && string(seqs[n.Nam][i]) != "N" {
						fmt.Println(string(seqs[n.Nam][i]))
						os.Exit(0)
					}
				}
				for _, j := range models[0].CharMap[string(seqs[n.Nam][i])] {
					n.Data[count][j] = 1.0
				}
				count++
			}
		}
	}

	start := time.Now()
	// calc likelihood
	w := 10
	if nsites < w {
		w = nsites
	}
	fmt.Println(t.Rt.Newick(true))

	l := gophy.PCalcLikePatternsMul(t, models, nodemodels, patternval, *wks)
	fmt.Println("lnL:", l)

	//optimize branch lengths
	fmt.Println("optimize bl1")
	gophy.OptimizeBLSMul(t, models, nodemodels, patternval, 10)
	//optimize model
	fmt.Println("optimize model")
	gophy.OptimizeGTRMul(t, models, nodemodels, patternval, 10)
	//optimize branch lengths
	fmt.Println("optimize bl2")
	gophy.OptimizeBLSMul(t, models, nodemodels, patternval, 10)
	end := time.Now()
	fmt.Fprintln(os.Stderr, end.Sub(start))
	//print the matrix
	for _, x := range models {
		for i := 0; i < 4; i++ {
			for j := 0; j < 4; j++ {
				if j < i {
					fmt.Print("- ")
				} else {
					fmt.Print(x.R.At(i, j), " ")
				}
			}
			fmt.Print("\n")
		}
		//print the basefreqs
		fmt.Println(x.BF)
	}
	fmt.Println(t.Rt.Newick(true))
	l = gophy.PCalcLikePatternsMul(t, models, nodemodels, patternval, *wks)
	fmt.Println("lnL:", l)
}
