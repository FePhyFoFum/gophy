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
	fmt.Fprintln(os.Stderr, "basefreqs --", basefreqs)
	fmt.Fprintln(os.Stderr, "modelparams --", modelparams)
	/*
		fmt.Fprintln(os.Stderr, " - ", modelparams[0], modelparams[1], modelparams[2])
		fmt.Fprintln(os.Stderr, modelparams[0], " - ", modelparams[3], modelparams[4])
		fmt.Fprintln(os.Stderr, modelparams[1], modelparams[3], " - ", 1.0)
		fmt.Fprintln(os.Stderr, modelparams[2], modelparams[4], 1.0, "-")
	*/
}

func getModels(t *gophy.Tree, single bool) (nodemodels map[*gophy.Node]int, maxint int) {
	nodemodels = map[*gophy.Node]int{}
	for _, i := range t.Pre {
		if len(i.Nam) == 0 {
			nodemodels[i] = 0
		} else if len(i.Chs) == 0 { //if parent has the num and it is tip, then it gets the num
			if strings.Contains(i.Nam, "#") {
				spls := strings.Split(i.Nam, "#")
				i.Nam = spls[0]
				if single == false {
					num64, _ := strconv.ParseInt(spls[1], 0, 0)
					num := int(num64)
					nodemodels[i] = num
					if num > maxint {
						maxint = num
					}
					fmt.Println(i, num)
				} else {
					nodemodels[i] = 0
				}
			} else {
				nodemodels[i] = nodemodels[i.Par]
			}

		} else if i.Nam[0] == '#' {
			if single == false {
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
		} else {
			nodemodels[i] = 0
		}
	}
	return
}

func main() {
	rand.Seed(uint64(time.Now().UTC().UnixNano()))
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
	indc := flag.Bool("c", false, "ind. comp.")
	indr := flag.Bool("r", false, "ind. rate mat.")
	opt := flag.Bool("o", false, "optimized?")
	//gam := flag.Bool("g", false, "gamma rates?")
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
	compfree := false
	if *indc {
		fmt.Println("comp free")
		compfree = true
	}
	ratefree := false
	if *indr {
		fmt.Println("rate")
		ratefree = true
	}
	single := false
	if ratefree == false && compfree == false {
		single = true
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

	//get number of models from tree labels
	nodemodels, maxint := getModels(t, single)
	fmt.Println("number of models:", maxint+1)
	//read a seq file
	nsites := 0
	seqs := map[string][]string{}
	seqnames := make([]string, 0)
	mseqs, numstates := gophy.ReadMSeqsFromFile(*afn)
	for _, i := range mseqs {
		seqs[i.NM] = i.SQs
		seqnames = append(seqnames, i.NM)
		nsites = len(i.SQs)
	}
	models := make([]gophy.StateModel, maxint+1)
	for i := 0; i <= maxint; i++ {
		// start model
		x := gophy.NewMultStateModel()
		x.NumStates = numstates
		x.SetMap()
		//empirical freqs for just the relevant seqs
		tseqs := make([]gophy.MSeq, 0)
		for tn := range nodemodels {
			if len(tn.Chs) == 0 && nodemodels[tn] == i {
				for tnn := range mseqs {
					if mseqs[tnn].NM == tn.Nam {
						tseqs = append(tseqs, mseqs[tnn])
					}
				}
			}
		}
		bf := gophy.GetEmpiricalBaseFreqsMS(tseqs, x.NumStates)
		x.SetBaseFreqs(bf)
		x.EBF = x.BF
		// model things
		modelparams := make([]float64, ((((x.NumStates * x.NumStates) - x.NumStates) / 2.) - 1))
		for m := range modelparams {
			modelparams[m] = 1.0
		}
		printModel(modelparams, bf)
		x.SetScaledRateMatrix(modelparams, true)
		fmt.Println(x.R)
		x.SetupQGTR()
		models[i] = x
	}
	// get the site patternas
	patterns, patternsint, gapsites, constant, uninformative, _ := gophy.GetSitePatternsMS(mseqs, models[0].GetCharMap(), models[0].GetNumStates())
	patternval, _ := gophy.PreparePatternVecsMS(t, patternsint, seqs, models[0].GetCharMap(), models[0].GetNumStates())
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

	l := gophy.PCalcLikePatternsMSMUL(t, models, nodemodels, patternval, *wks)
	fmt.Println("starting lnL:", l)

	if *opt {
		//optimize branch lengths
		fmt.Println("optimize bl")
		//optimize branch lengths
		fmt.Println("getting starting branch lengths (parsimony)")
		s := gophy.PCalcSankParsPatternsMultState(t, numstates, patternval, *wks)
		fmt.Println("sank:", s, nsites)
		gophy.EstParsBLMultState(t, numstates, patternval, nsites)

		fmt.Println("start:\n" + t.Rt.Newick(true) + ";")
		gophy.OptimizeBLNRMSMul(t, models, nodemodels, patternval, *wks)
		l = gophy.PCalcLikePatternsMSMUL(t, models, nodemodels, patternval, *wks)
		fmt.Println("ln:", l)
		//optimize model
		if compfree && ratefree == false {
			fmt.Println("optimize model: comp free , rate shared")
			gophy.OptimizeGTRMSCompSharedRM(t, models, nodemodels, patternval, *wks)
		}
		if compfree && ratefree {
			gophy.OptimizeGTRBPMSMul(t, models, nodemodels, patternval, *wks)
		}
		if compfree == false && ratefree {
			//comp emp
			gophy.OptimizeGTRMSMul(t, models, nodemodels, patternval, *wks)
		}
		if single {
			fmt.Println("optimize model: single")
			//comp emp
			gophy.OptimizeGTRMSMul(t, models, nodemodels, patternval, *wks)
		}
		//optimize branch lengths
		fmt.Println("optimize bl2")
		gophy.OptimizeBLNRMSMul(t, models, nodemodels, patternval, *wks)

		//print the matrix
		fmt.Println(t.Rt.Newick(true))
		l = gophy.PCalcLikePatternsMSMUL(t, models, nodemodels, patternval, *wks)
		fmt.Println("lnL:", l)
	}
}

/*
	gophy.TritomyRoot(t)
	fmt.Println("optimize model")
	gophy.OptimizeGTRMul(t, models, nodemodels, patternval, 4)
	l := gophy.PCalcLikePatternsMul(t, models, nodemodels, patternval, *wks)
	l = -gophy.OptimizeBLSMul(t, models, nodemodels, patternval, 10)
	fmt.Println("lnL:", l)
	os.Exit(0)
	moves := gophy.NNIMoves(t)
	for i, j := range moves {
		gophy.SwapBranch(j[0], j[1])
		fmt.Println(" ", i, t.Rt.Newick(false)+";")
		fmt.Println("optimize bl")
		l = -gophy.OptimizeBLSMul(t, models, nodemodels, patternval, 10)
		fmt.Println("  bllnL:", l)
		fmt.Println("  " + t.Rt.Newick(true) + ";")
		gophy.SwapBranch(j[0], j[1])
	}
	//reroot
	nrt := t.Rt.Chs[0].Chs[1]
	fmt.Println(nrt)
	gophy.Reroot(nrt, t)
	gophy.TritomyRoot(t)
	t.Instantiate(t.Rt)
	t.Rt.Data = make([][]float64, len(patternsint))
	for i := 0; i < len(patternsint); i++ {
		t.Rt.Data[i] = []float64{0.0, 0.0, 0.0, 0.0}
	}
	//end reroot
	l = gophy.PCalcLikePatternsMul(t, models, nodemodels, patternval, *wks)
	fmt.Println("lnL:", l)
	fmt.Println(t.Rt.Newick(true))
	moves = gophy.NNIMoves(t)
	for i, j := range moves {
		gophy.SwapBranch(j[0], j[1])
		fmt.Println(" ", i, t.Rt.Newick(false)+";")
		fmt.Println("optimize bl")
		l = -gophy.OptimizeBLSMul(t, models, nodemodels, patternval, 10)
		fmt.Println("  bllnL:", l)
		fmt.Println("  " + t.Rt.Newick(true) + ";")
		gophy.SwapBranch(j[0], j[1])
	}
	fmt.Println(t.Rt.Newick(true) + ";")
	os.Exit(0)
*/
