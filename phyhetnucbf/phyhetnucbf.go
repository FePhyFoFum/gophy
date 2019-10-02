package main

import (
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
	"gonum.org/v1/gonum/floats"

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
	models = []*gophy.DNAModel{}
	for _, i := range curmodels {
		models = append(models, i)
	}
	models = append(models, y)
	nodemodels = make(map[*gophy.Node]int)
	for i, j := range curnodemodels {
		nodemodels[i] = j
	}
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
	uc1 := flag.Bool("ue", false, "uncertainty existence")
	uc2 := flag.Bool("ul", false, "uncertainty location")
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
	fmt.Fprintln(os.Stderr, "lnL:", l)
	useLog := false //slower
	if math.IsInf(l, -1) {
		fmt.Println("initial value problem, switch to big")
		useLog = true
		l = gophy.PCalcLogLikePatterns(t, x, patternval, *wks)
		fmt.Fprintln(os.Stderr, "lnL:", l)
	}
	gophy.OptimizeGTR(t, x, patternval, useLog, *wks)
	gophy.OptimizeBF(t, x, patternval, useLog, *wks)
	if useLog {
		l = gophy.PCalcLogLikePatterns(t, x, patternval, *wks)
	} else {
		l = gophy.PCalcLikePatterns(t, x, patternval, *wks)
	}
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
	saic := gophy.CalcBIC(l, numbaseparams, nsites)
	currentaic := saic
	//start with the sorted set
	fmt.Fprintln(os.Stderr, "starting BIC:", saic)
	cur := 1
	for _, i := range t.Post {
		if i == t.Rt {
			continue
		}
		if getVisible(i) > minset {
			fmt.Fprint(os.Stderr, "\rOn "+strconv.Itoa(cur)+"/"+strconv.Itoa(count))
			y := allmodels[modelmap[i]]
			models := []*gophy.DNAModel{allmodels[0], y}
			gophy.OptimizeBFSubClade(t, i, false, y, patternval, useLog, *wks)
			for _, j := range t.Post {
				nodemodels[j] = 0
			}
			for _, j := range i.PreorderArray() {
				nodemodels[j] = 1
			}
			lm := 1.0
			if useLog {
				lm = gophy.PCalcLogLikePatternsMul(t, models, nodemodels, patternval, *wks)
			} else {
				lm = gophy.PCalcLikePatternsMul(t, models, nodemodels, patternval, *wks)
			}
			naic := gophy.CalcBIC(lm, numbaseparams+3., nsites)
			nodevalues[i] = naic
			cur++
		}
	}
	keys := sortAicMap(nodevalues)
	curmodels := []*gophy.DNAModel{allmodels[0]}
	curnodemodels := make(map[*gophy.Node]int)
	//get the final model configuration
	for _, j := range t.Post {
		curnodemodels[j] = 0
	}
	fmt.Fprintln(os.Stderr, "\rfinal est")
	curparams := numbaseparams
	modelnodes := []*gophy.Node{}
	cur = 0
	for _, k := range keys {
		fmt.Fprintln(os.Stderr, "On "+strconv.Itoa(cur)+"/"+strconv.Itoa(count))
		if getVisible(k.Key) > minset {
			//fmt.Println(k.Key, k.Value)
			//shortcut, if k.Value is worse than the best by > 10 don't consider it
			if k.Value-saic > 10 { // 10 is arbitrary. pick another measure
				fmt.Println(k.Value)
				cur++
				continue
			}
			y := allmodels[modelmap[k.Key]]
			testmodels, testnodemodels, modelint := getNodeModels(curmodels, curnodemodels, y, k.Key)
			//need to be able to send better starting points so it doesn't take as long
			lm := 1.0
			if useLog {
				gophy.OptimizeGTRCompSharedRM(t, testmodels, testnodemodels, true, patternval, true, *wks)
				lm = gophy.PCalcLogLikePatternsMul(t, testmodels, testnodemodels, patternval, *wks)
			} else {
				gophy.OptimizeGTRCompSharedRM(t, testmodels, testnodemodels, true, patternval, false, *wks)
				lm = gophy.PCalcLikePatternsMul(t, testmodels, testnodemodels, patternval, *wks)
			}
			naic := gophy.CalcBIC(lm, curparams+3., nsites)
			fmt.Fprintln(os.Stderr, " -- ", naic, lm, currentaic)
			if naic < currentaic {
				fmt.Fprintln(os.Stderr, naic, "<", currentaic, curparams+3)
				currentaic = naic
				curparams += 3
				curmodels = testmodels
				curnodemodels = testnodemodels
				modelnodes = append(modelnodes, k.Key)
			} else {
				unmarkModel(k.Key, modelint)
			}
		}
		cur++
	}
	//at this point we should have the good configuration
	//uncertainty
	// existence, check with and without the
	if *uc1 {
		mainaic := 1.0 //same as math.Exp((currentaic - currentaic) / 2)
		for _, i := range modelnodes {
			fmt.Fprintln(os.Stderr, i, i.IData["shift"])
			//unmark just the modelnodes
			testnodemodels := make(map[*gophy.Node]int)
			testmodels := []*gophy.DNAModel{}
			for _, j := range curmodels {
				testmodels = append(testmodels, j.DeepCopyDNAModel())
			}
			for k, j := range curnodemodels {
				testnodemodels[k] = j
			}
			for _, k := range i.PreorderArray() {
				if i.IData["shift"] == k.IData["shift"] {
					testnodemodels[k] = 0
				}
			}
			//get aic, could reestimate the model but not right now, if this happens it needs to be copied so that we can
			//   just reestimate one model
			//   go back to the original as well
			tlm := 0.0
			if useLog {
				gophy.OptimizeGTRCompSharedRMSingleModel(t, testmodels,
					testnodemodels, true, 0, patternval, true, *wks)
				tlm = gophy.PCalcLogLikePatternsMul(t, testmodels, testnodemodels, patternval, *wks)
			} else {
				gophy.OptimizeGTRCompSharedRMSingleModel(t, testmodels,
					testnodemodels, true, 0, patternval, false, *wks)
				tlm = gophy.PCalcLikePatternsMul(t, testmodels, testnodemodels, patternval, *wks)
			}
			taic := gophy.CalcBIC(tlm, curparams-3., nsites)
			tstat := math.Exp((currentaic - taic) / 2)
			i.FData["uc1"] = mainaic / (mainaic + tstat)
			fmt.Fprintln(os.Stderr, currentaic, taic, i.FData["uc1"])
		}
	}
	//uncertainty
	// location
	if *uc2 {
		mainaic := 1.0 //same as math.Exp((currentaic - currentaic) / 2)
		for _, i := range modelnodes {
			fmt.Fprintln(os.Stderr, i, i.IData["shift"])
			testnodes := []*gophy.Node{}
			if _, ok := i.Par.IData["shift"]; !ok {
				testnodes = append(testnodes, i.Par)
			}
			for _, j := range i.Chs {
				if j.IData["shift"] == i.IData["shift"] {
					testnodes = append(testnodes, j)
				}
			}
			tevals := []float64{mainaic}
			for _, ts := range testnodes {
				//unmark just the modelnodes
				testnodemodels := make(map[*gophy.Node]int)
				testmodels := []*gophy.DNAModel{}
				for _, j := range curmodels {
					testmodels = append(testmodels, j.DeepCopyDNAModel())
				}
				for k, j := range curnodemodels {
					testnodemodels[k] = j
				}
				if ts != i.Par {
					if _, ok := i.Par.IData["shift"]; !ok {
						testnodemodels[i] = 0
					} else {
						testnodemodels[i] = testnodemodels[i.Par]
					}
					//get other child
					for _, k := range ts.GetSib().PreorderArray() {
						testnodemodels[k] = 0
					}
				} else { // the node is the parent
					testnodemodels[ts] = i.IData["shift"]
					for _, k := range ts.PreorderArray() {
						if _, ok := k.IData["shift"]; !ok {
							testnodemodels[k] = i.IData["shift"]
						}
					}
				}
				//get aic, could reestimate the model but not right now, if this happens it needs to be copied so that we can
				//   go back to the original as well
				tlm := 0.0
				if useLog {
					gophy.OptimizeGTRCompSharedRMSingleModel(t, testmodels,
						testnodemodels, true, i.IData["shift"], patternval, true, *wks)
					tlm = gophy.PCalcLogLikePatternsMul(t, testmodels, testnodemodels, patternval, *wks)
				} else {
					gophy.OptimizeGTRCompSharedRMSingleModel(t, testmodels,
						testnodemodels, true, i.IData["shift"], patternval, false, *wks)
					tlm = gophy.PCalcLikePatternsMul(t, testmodels, testnodemodels, patternval, *wks)
				}
				taic := gophy.CalcBIC(tlm, curparams, nsites)
				tstat := math.Exp((currentaic - taic) / 2)
				tevals = append(tevals, tstat)
			}
			waic := mainaic / floats.Sum(tevals)
			fmt.Fprintln(os.Stderr, currentaic, waic)
		}
	}
	moveMarksToLabels(t.Rt)
	fmt.Fprintln(os.Stderr, "Final models")
	fmt.Fprintln(os.Stderr, "-------")
	for i, j := range curmodels {
		fmt.Fprintln(os.Stderr, i, j.BF)
	}
	//make a nexus so it is easier to read in the fig-tree
	fmt.Println(t.Rt.Newick(true) + ";")
}
