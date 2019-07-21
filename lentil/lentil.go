package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"

	"github.com/FePhyFoFum/gophy"
	"gonum.org/v1/gonum/stat"
)

func main() {
	mfn := flag.String("m", "", "map filename")
	tfn := flag.String("t", "", "tree filename")
	wks := flag.Int("w", 2, "how many workers?")
	flag.Parse()
	if len(*mfn) == 0 {
		fmt.Fprintln(os.Stderr, "need a map filename (-m)")
		os.Exit(1)
	}
	if len(*tfn) == 0 {
		fmt.Fprintln(os.Stderr, "need a tree filename (-t)")
		os.Exit(1)
	}
	maptips := make(map[string]int)
	mapints := make(map[int]string)
	numtips := 0
	//read map tree
	fmt.Fprint(os.Stderr, "reading map tree\n")
	fc, err := os.Open(*mfn)
	if err != nil {
		fmt.Println(err)
	}
	defer fc.Close()
	csc := bufio.NewScanner(fc)
	var mapt gophy.Tree
	for csc.Scan() {
		ln := csc.Text()
		if len(ln) < 2 {
			continue
		}
		rt := gophy.ReadNewickString(ln)
		mapt.Instantiate(rt)
		for _, n := range mapt.Tips {
			if _, ok := maptips[n.Nam]; !ok {
				maptips[n.Nam] = numtips
				mapints[numtips] = n.Nam
				numtips++
			}
		}
	}
	//read tree set (all the tips have to be in the map)
	trees := make([]gophy.Tree, 0)
	ntrees := 0
	fc, err = os.Open(*tfn)
	if err != nil {
		fmt.Println(err)
	}
	defer fc.Close()
	scanner := bufio.NewReader(fc)
	for {
		ln, err := scanner.ReadString('\n')
		if len(ln) > 0 {
			fmt.Fprintf(os.Stderr, "\rreading tree %d", ntrees+1)
			rt := gophy.ReadNewickString(ln)
			var t gophy.Tree
			t.Instantiate(rt)
			trees = append(trees, t)
			ntrees++
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Printf("read %s bytes: %v", ln, err)
			break
		}
	}
	fmt.Fprintf(os.Stderr, "\n")

	// do analyses
	spQuartsD, qToN, bigQuartetArr, tipLens, tipTrees := prepareSpTree(mapt, maptips)
	//processGeneTrees(mapt, trees, maptips, mapints, spQuartsD, bigQuartetArr, tipLens, tipTrees)
	pprocessGeneTrees(*wks, mapt, trees, maptips, mapints, spQuartsD, bigQuartetArr, tipLens, tipTrees)
	calcValues(mapt, qToN, spQuartsD, tipLens)
	writeOutput(mapt)
}

func prepareSpTree(specTree gophy.Tree, maptips map[string]int) (spQuartsD map[int][]float64,
	qToN map[int]*gophy.Node, bigQuartetArr []gophy.Quartet, tipLens map[int][]float64,
	tipTrees map[int]map[int]bool) {
	spQuartsD = map[int][]float64{}
	qToN = map[int]*gophy.Node{}
	bigQuartetArr = []gophy.Quartet{}
	qcount := 0
	tipLens = map[int][]float64{}     //node int
	tipTrees = map[int]map[int]bool{} //node int, tree int map
	for _, i := range specTree.Pre {
		if i == specTree.Rt {
			continue
		}
		t, err := gophy.GetQuartet(i, specTree, maptips)
		if err == nil {
			t.Index = qcount
			spQuartsD[t.Index] = []float64{}
			qToN[t.Index] = i
			bigQuartetArr = append(bigQuartetArr, t)
			qcount++
		} else { //tips don't have quartets so just store it in the tree
			i.Num = qcount
			tipLens[qcount] = []float64{}
			tipTrees[qcount] = map[int]bool{}
			qToN[qcount] = i
			bigQuartetArr = append(bigQuartetArr, gophy.Quartet{})
			qcount++
		}
	}
	return
}

func processGeneTrees(specTree gophy.Tree, trees []gophy.Tree, maptips map[string]int,
	mapints map[int]string, spQuartsD map[int][]float64, bigQuartetArr []gophy.Quartet,
	tipLens map[int][]float64, tipTrees map[int]map[int]bool) {
	for count, tree1 := range trees {
		//fmt.Fprintln(os.Stderr, count)
		fmt.Fprintf(os.Stderr, "\rprocessing tree %d/%d", count+1, len(trees))
		//for each species tree quart
		for i := range spQuartsD {
			// for each node in the gene tree
			for _, j := range tree1.Pre {
				gqu, err := gophy.GetQuartet(j, tree1, maptips)
				if err == nil {
					// if the gene tree quartet matches the species tree quart
					if bigQuartetArr[i].Match(gqu) {
						//add the length to the species tree quart
						spQuartsD[i] = append(spQuartsD[i], j.Len)
						//special check for tip lengths
						keep := checkTips(gqu, bigQuartetArr[i])
						for k := range keep {
							ln, err1 := tree1.GetTipByName(mapints[k])
							if err1 != nil {
								fmt.Println(err1)
								os.Exit(0)
							}
							st, err2 := specTree.GetTipByName(mapints[k])
							if err2 != nil {
								fmt.Println(err2)
								os.Exit(0)
							}
							if _, ok := tipTrees[st.Num][count]; !ok {
								tipLens[st.Num] = append(tipLens[st.Num], ln.Len)
								tipTrees[st.Num][count] = true
							}
						}
						//end tip check
					}
				}
			}
		}
	}
	fmt.Fprintf(os.Stderr, "\n")
}

//pprocessGeneTrees parallel process
func pprocessGeneTrees(workers int, specTree gophy.Tree, trees []gophy.Tree, maptips map[string]int,
	mapints map[int]string, spQuartsD map[int][]float64, bigQuartetArr []gophy.Quartet,
	tipLens map[int][]float64, tipTrees map[int]map[int]bool) {
	for count, tree1 := range trees {
		//fmt.Fprintln(os.Stderr, count)
		fmt.Fprintf(os.Stderr, "\rprocessing tree %d/%d", count+1, len(trees))
		//for each species tree quart
		//job is the i_th spQuartsD and j_th tree1.Pre
		jobs := make(chan []int, len(spQuartsD)*len(tree1.Tips)*2)
		results := make(chan GeneWorkerReturn, len(spQuartsD)*len(tree1.Tips)*2)
		for w := 1; w <= workers; w++ {
			go pprocessGeneTreesWorker(trees, specTree, bigQuartetArr, maptips,
				mapints, jobs, results)
		}
		njobs := 0
		for i := range spQuartsD {
			for j := range tree1.Pre {
				if tree1.Pre[j] != tree1.Rt && len(tree1.Pre[j].Chs) > 0 {
					jobs <- []int{i, j, count}
					njobs++
				}
			}
		}
		close(jobs)
		for i := 0; i < njobs; i++ {
			val := <-results
			if val.match == true {
				spQuartsD[val.spi] = append(spQuartsD[val.spi], val.jlen)
				for j := range val.stNum {
					if _, ok := tipTrees[val.stNum[j]][val.treen[j]]; !ok {
						tipLens[val.stNum[j]] = append(tipLens[val.stNum[j]], val.lnlen[j])
						tipTrees[val.stNum[j]][val.treen[j]] = true
					}
				}
			}
		}
	}
	fmt.Fprintf(os.Stderr, "\n")
}

//GeneWorkerReturn just returning the necessary values
type GeneWorkerReturn struct {
	match bool
	spi   int
	jlen  float64
	stNum []int
	lnlen []float64
	treen []int
}

//pprocessGeneTreesWorker parallel process worker
func pprocessGeneTreesWorker(trees []gophy.Tree, specTree gophy.Tree,
	bigQuartetArr []gophy.Quartet, maptips map[string]int, mapints map[int]string,
	jobs <-chan []int, results chan<- GeneWorkerReturn) {
	for j := range jobs {
		gwr := GeneWorkerReturn{match: false}
		spi, noden, treen := j[0], j[1], j[2]
		gqu, _ := gophy.GetQuartet(trees[treen].Pre[noden], trees[treen], maptips)
		// if the gene tree quartet matches the species tree quart
		if bigQuartetArr[spi].Match(gqu) {
			//			fmt.Println(bigQuartetArr[spi].StringWithNames(mapints), gqu.StringWithNames(mapints))
			gwr.match = true
			//add the length to the species tree quart
			gwr.spi = spi
			gwr.jlen = trees[treen].Pre[noden].Len
			//special check for tip lengths
			keep := checkTips(gqu, bigQuartetArr[spi])
			gwr.stNum = make([]int, len(keep))
			gwr.treen = make([]int, len(keep))
			gwr.lnlen = make([]float64, len(keep))
			count := 0
			for k := range keep {
				ln, _ := trees[treen].GetTipByName(mapints[k])
				st, _ := specTree.GetTipByName(mapints[k])
				gwr.stNum[count] = st.Num
				gwr.treen[count] = treen
				gwr.lnlen[count] = ln.Len
				count++
			}
			//end tip check
		}
		results <- gwr
	}
}

func checkTips(gqu gophy.Quartet, bqi gophy.Quartet) map[int]bool {
	check := map[int]bool{} // tip
	keep := map[int]bool{}  // keepers need to change to strings later
	for _, k := range gqu.Lts {
		if len(k) == 1 { //potential tip info
			for lk := range k {
				check[lk] = true
			}
		}
	}
	for _, k := range gqu.Rts {
		if len(k) == 1 { //potential tip info
			for lk := range k {
				check[lk] = true
			}
		}
	}
	for _, k := range bqi.Lts {
		if len(k) == 1 && gophy.IntMapIntersects(check, k) {
			for lk := range k {
				keep[lk] = true
			}
		}
	}
	for _, k := range bqi.Rts {
		if len(k) == 1 && gophy.IntMapIntersects(check, k) {
			for lk := range k {
				keep[lk] = true
			}
		}
	}
	return keep
}

func setNodeVals(nd *gophy.Node, vals []float64) {
	nd.FData["mean"] = stat.Mean(vals, nil)
	nd.FData["median"] = gophy.MedianF(vals)
	nd.FData["min"] = gophy.MinF(vals)
	nd.FData["max"] = gophy.MaxF(vals)
	nd.FData["supp"] = float64(len(vals))
	lc, hc := vals[0], vals[0]
	if len(vals) > 1 {
		lc, hc = gophy.ConfInt95NormalF(vals)
	}
	nd.FData["CONFH"] = hc
	nd.FData["CONFL"] = math.Max(0, lc)
}

func calcValues(specTree gophy.Tree, qToN map[int]*gophy.Node, spQuartsD map[int][]float64,
	tipLens map[int][]float64) {
	for i, j := range qToN {
		var a []float64
		if len(j.Chs) == 0 {
			a = tipLens[i]
		} else {
			a = spQuartsD[i]
		}
		if len(a) > 0 {
			setNodeVals(j, a)
		}
	}
}

func writeOutput(specTree gophy.Tree) {
	calcsF := []string{"mean", "median", "min", "max", "CONFH", "CONFL", "supp"}
	for _, i := range calcsF {
		fmt.Print(specTree.Rt.NewickFloatBL(i) + ";\n")
	}
}
