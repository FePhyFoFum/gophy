package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"runtime/pprof"
	"sort"
	"strconv"
	"time"

	"github.com/FePhyFoFum/gophy"

	"golang.org/x/exp/rand"
)

// this will determine how many tips are not already covered by another model
func getVisible(nd *gophy.Node) (numvisible int) {
	dontinclude := map[*gophy.Node]bool{}
	for _, i := range nd.PreorderArray() {
		if _, ok := i.IData["shift"]; ok {
			for _, j := range i.GetTips() {
				dontinclude[j] = true
			}
		}
	}
	numvisible = len(nd.GetTips()) - len(dontinclude)
	return
}

type keyVal struct {
	Key   *gophy.Node
	Value float64
}

//returns a list with the smallest values
func sortAicMap(nodevalues map[*gophy.Node]float64) (ss []keyVal) {
	//sort the values
	ss = []keyVal{}
	for k, v := range nodevalues {
		ss = append(ss, keyVal{k, v})
	}

	sort.Slice(ss, func(i, j int) bool {
		return ss[i].Value < ss[j].Value
	})
	return
}

//mark the model
func markModel(nd *gophy.Node, modelnum int) {
	nd.IData["shift"] = modelnum
	nd.Nam = strconv.Itoa(modelnum)
	for _, i := range nd.PreorderArray() {
		if _, ok := i.IData["shift"]; !ok {
			i.IData["shift"] = modelnum
		}
	}
}

//mark the model
func unmarkModel(nd *gophy.Node, modelnum int) {
	delete(nd.IData, "shift")
	nd.Nam = ""
	for _, i := range nd.PreorderArray() {
		if _, ok := i.IData["shift"]; ok {
			if i.IData["shift"] == modelnum {
				delete(i.IData, "shift")
			}
		}
	}
}

func moveMarksToLabels(nd *gophy.Node) {
	for _, i := range nd.PreorderArray() {
		if len(i.Chs) == 0 {
			continue
		}
		if _, ok := i.IData["shift"]; ok {
			i.Nam = strconv.Itoa(i.IData["shift"])
		}
	}
}

func getNodeModels(curmodels []*gophy.DNAModel, curnodemodels map[*gophy.Node]int,
	y *gophy.DNAModel, nd *gophy.Node) (models []*gophy.DNAModel, nodemodels map[*gophy.Node]int, modelint int) {
	modelint = len(curmodels)
	markModel(nd, modelint)
	models = curmodels
	models = append(models, y)
	nodemodels = curnodemodels
	for _, i := range nd.PreorderArray() {
		if i.IData["shift"] == modelint {
			nodemodels[i] = modelint
		}
	}
	return
}

func main() {
	rand.Seed(uint64(time.Now().UTC().UnixNano()))
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
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
	minset := 10 // have to have at least 8 tips in there

	//read tree
	t := gophy.ReadTreeFromFile(*tfn)
	seqs, patternsint, nsites, bf := gophy.ReadPatternsSeqsFromFile(*afn)
	patternval, _ := gophy.PreparePatternVecs(t, patternsint, seqs)

	// model things
	modelparams := make([]float64, 5)
	for i := range modelparams {
		modelparams[i] = 1.0
	}
	//root model
	x := gophy.NewDNAModel()
	x.SetMap()
	fmt.Fprintln(os.Stderr, "emp:", bf)
	x.SetBaseFreqs(bf)
	x.SetRateMatrix(modelparams)
	x.SetupQGTR()
	l := gophy.PCalcLikePatterns(t, x, patternval, *wks)
	gophy.OptimizeGTR(t, x, patternval, *wks)
	gophy.OptimizeBF(t, x, patternval, *wks)
	l = gophy.PCalcLikePatterns(t, x, patternval, *wks)
	fmt.Fprintln(os.Stderr, "lnL:", l)

	//for each node in the tree if the number of tips is > minset
	allmodels := make([]*gophy.DNAModel, 1) //number of clades that have enough taxa and aren't the root
	modelmap := make(map[*gophy.Node]int)
	allmodels[0] = x
	count := 1
	for _, i := range t.Post {
		if i == t.Rt {
			continue
		}
		if len(i.GetTips()) > minset {
			y := gophy.NewDNAModel()
			y.SetMap()
			y.SetBaseFreqs(bf)
			y.SetRateMatrix(modelparams)
			y.SetupQGTR()
			modelmap[i] = count
			allmodels = append(allmodels, y)
			count++
		}
	}
	fmt.Fprintln(os.Stderr, "num of models:", count)

	nodemodels := make(map[*gophy.Node]int)
	numbaseparams := ((2 * float64(len(t.Tips))) - 3.) + 5. + 3. //tree ones and GTR and one BF
	nodevalues := make(map[*gophy.Node]float64)                  // key node, value aicc
	saic := gophy.CalcAICC(l, numbaseparams, nsites)
	fmt.Fprintln(os.Stderr, "starting AICc:", saic)
	cur := 1
	for _, i := range t.Post {
		if i == t.Rt {
			continue
		}
		if getVisible(i) > minset {
			fmt.Fprint(os.Stderr, "\rOn "+strconv.Itoa(cur)+"/"+strconv.Itoa(count))
			y := allmodels[modelmap[i]]
			models := []*gophy.DNAModel{allmodels[0], y}
			gophy.OptimizeBFRMSubClade(t, i, false, y, patternval, *wks)
			for _, j := range t.Post {
				nodemodels[j] = 0
			}
			for _, j := range i.PreorderArray() {
				nodemodels[j] = 1
			}
			lm := gophy.PCalcLikePatternsMul(t, models, nodemodels, patternval, *wks)
			naic := gophy.CalcAICC(lm, numbaseparams+3., nsites)
			nodevalues[i] = naic
			cur++
		}
	}
	keys := sortAicMap(nodevalues)
	currentaic := saic
	curmodels := []*gophy.DNAModel{allmodels[0]}
	curnodemodels := make(map[*gophy.Node]int)
	//initialize things
	for _, j := range t.Post {
		curnodemodels[j] = 0
	}
	fmt.Fprintln(os.Stderr, "\rfinal est")
	curparams := numbaseparams
	cur = 0
	for _, k := range keys {
		fmt.Fprintln(os.Stderr, "On "+strconv.Itoa(cur)+"/"+strconv.Itoa(count))
		if getVisible(k.Key) > minset {
			fmt.Println(k.Key, k.Value)
			y := allmodels[modelmap[k.Key]]
			testmodels, testnodemodels, modelint := getNodeModels(curmodels, curnodemodels, y, k.Key)
			//need to be able to send better starting points so it doesn't take as long
			gophy.OptimizeGTRBPMul(t, testmodels, testnodemodels, true, patternval, *wks)
			lm := gophy.PCalcLikePatternsMul(t, testmodels, testnodemodels, patternval, *wks)
			naic := gophy.CalcAICC(lm, curparams+8., nsites)
			fmt.Fprintln(os.Stderr, " -- ", naic, lm, currentaic)
			if naic < currentaic {
				fmt.Fprintln(os.Stderr, naic, "<", currentaic, curparams+8)
				//check for nested shifts in this one.
				//use OptimizeGTRCompSharedRMSubClade
				currentaic = naic
				curparams += 8
				curmodels = testmodels
				curnodemodels = testnodemodels
			} else {
				unmarkModel(k.Key, modelint)
			}
		}
		cur++
	}
	moveMarksToLabels(t.Rt)
	fmt.Fprintln(os.Stderr, "Final models")
	fmt.Fprintln(os.Stderr, "-------")
	for i, j := range curmodels {
		fmt.Fprintln(os.Stderr, i, j.BF)
		fmt.Fprintln(os.Stderr, i, j.R)
	}
	fmt.Println(t.Rt.Newick(true) + ";")
}
