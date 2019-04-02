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
	patternval := gophy.PreparePatternVecs(t, patternsint, seqs)
	//list of sites
	fmt.Fprintln(os.Stderr, "nsites:", nsites)
	fmt.Fprintln(os.Stderr, "patterns:", len(patterns), len(patternsint))
	fmt.Fprintln(os.Stderr, "onlygaps:", len(gapsites))
	fmt.Fprintln(os.Stderr, "constant:", len(constant))
	fmt.Fprintln(os.Stderr, "uninformative:", len(uninformative))

	start := time.Now()
	// calc likelihood
	w := 10
	if nsites < w {
		w = nsites
	}

	l := gophy.PCalcLikePatternsMul(t, models, nodemodels, patternval, *wks)
	fmt.Println("starting lnL:", l)

	//optimize branch lengths
	fmt.Println("optimize bl")
	//optimize branch lengths
	fmt.Println("getting starting branch lengths (parsimony)")
	s := gophy.PCalcSankParsPatterns(t, patternval, *wks)
	fmt.Println("sank:", s)
	gophy.EstParsBL(t, patternval, nsites)

	fmt.Println("start:\n" + t.Rt.Newick(true) + ";")
	gophy.OptimizeBLNRMult(t, models, nodemodels, patternval, 10)
	//optimize model

	if compfree && ratefree == false {
		fmt.Println("optimize model: comp free , rate shared")
		gophy.OptimizeGTRCompSharedRM(t, models, nodemodels, patternval, 10)
	}
	if compfree && ratefree {
		gophy.OptimizeGTRBPMul(t, models, nodemodels, patternval, 10)
	}
	if compfree == false && ratefree {
		//comp emp
		gophy.OptimizeGTRMul(t, models, nodemodels, patternval, 10)
	}
	if single {
		//comp emp
		gophy.OptimizeGTRMul(t, models, nodemodels, patternval, 10)
	}
	/*
		//optimize branch lengths
		fmt.Println("optimize bl2")
		gophy.OptimizeBLSMul(t, models, nodemodels, patternval, 10)
		//optimize model
		fmt.Println("optimize model2")
		gophy.OptimizeGTRMul(t, models, nodemodels, patternval, 10)
	*/

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
