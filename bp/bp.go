package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/FePhyFoFum/gophy"
)

// RunParams options for the run
type RunParams struct {
	BlCut          float64
	SupCut         float64
	TIgnore        []string
	IncludeTips    bool
	IncludeNodeMap bool
}

/*
 This will calculate many things things necessary for comparing
 bipartitions.
*/

func main() {
	rp := RunParams{BlCut: 0, SupCut: 0, TIgnore: nil}
	wks := flag.Int("w", 2, "how many workers?")
	cut := flag.Float64("scut", 0.0, " support cutoff (if support is present in the trees)")
	blcut := flag.Float64("bcut", 0.0, "branch length cutoff")
	oed := flag.Bool("e", false, "output edges?")
	oconf := flag.String("oconf", "", "run conflict by giving an output filename")
	rf := flag.Bool("rf", false, "run rf?")
	rfp := flag.Bool("rfp", false, "run rf (partial overlap)?")
	rfw := flag.Bool("rfw", false, "run rfw (weighted)?")
	rfwp := flag.Bool("rfwp", false, "run rfwp (weighted w/ missing taxa penalty)?")
	comp := flag.String("c", "", "compare biparts to those in this file")
	pca := flag.Bool("pca", false, "pairwise compare all trees (with a pool of biparts)")
	pc := flag.Bool("pc", false, "pairwise compare each tree (with each tree)")
	fn := flag.String("t", "", "tree filename")
	ig := flag.String("ig", "", "ignore these taxa (comma not space separated)")
	v := flag.Bool("v", false, "verbose results?")
	tv := flag.Bool("tv", false, "for the -c option. print the results on the comp tree")
	rng := flag.String("rng", "", "range of trees to check in a large tree file like -rng 0-100 for the first hundred")
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")
	memprofile := flag.String("memprofile", "", "write mem profile to file")
	flag.Parse()
	if len(os.Args) < 2 {
		flag.PrintDefaults()
		os.Exit(1)
	}
	//filename things
	if *wks > 1 {
		fmt.Fprintln(os.Stderr, "workers:", *wks)
	}
	if *cut > 0.0 {
		fmt.Fprintln(os.Stderr, "cutoff set to:", *cut)
		rp.SupCut = *cut
	}
	if *blcut > 0.0 {
		fmt.Fprintln(os.Stderr, "blcutoff set to:", *blcut)
		rp.BlCut = *blcut
	}
	if *rfw || *rfwp {
		rp.IncludeNodeMap = true
		rp.IncludeTips = true
	}
	if len(*fn) == 0 {
		fmt.Fprintln(os.Stderr, "need a filename")
		flag.PrintDefaults()
		os.Exit(1)
	}
	if *v == true {
		fmt.Fprintln(os.Stderr, "verbose")
	}
	if *tv == true {
		fmt.Fprintln(os.Stderr, "tree verbose")
	}

	//cpu profile code
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}
	//end profile code

	rngstart := 0
	rngstop := 0
	rngcheck := false
	if len(*rng) > 0 {
		rngs := strings.Split(*rng, "-")
		rv, err := strconv.Atoi(rngs[0])
		if err != nil {
			fmt.Fprintln(os.Stderr, "error processing ", rngs, " as range")
			os.Exit(0)
		}
		rngstart = rv
		rv, err = strconv.Atoi(rngs[1])
		if err != nil {
			fmt.Fprintln(os.Stderr, "error processing ", rngs, " as range")
			os.Exit(0)
		}
		rngstop = rv
		fmt.Fprintln(os.Stderr, "only looking at trees from", rngs[0], "-", rngs[1])
		rngcheck = true
	}
	ignore := []string{}
	if len(*ig) > 0 {
		ignore = strings.Split(*ig, ",")
		fmt.Fprintln(os.Stderr, "ignoring:", ignore)
	}
	rp.TIgnore = ignore
	f, err := os.Open(*fn)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	defer f.Close()
	scanner := bufio.NewReader(f)
	bps := make([]gophy.Bipart, 0) // list of all biparts
	bpts := make(map[int][]int)    // key tree index, value bipart index list
	ntrees := 0
	readtrees := 0
	trees := make([]gophy.Tree, 0)
	numtips := 0
	skipped := 0
	maptips := make(map[string]int)
	mapints := make(map[int]string)

	//start a timer
	start := time.Now()

	// reading the trees
	fmt.Fprint(os.Stderr, "reading trees\n")
	for {
		ln, err := scanner.ReadString('\n')
		if len(ln) > 0 {
			if rngcheck {
				if ntrees < rngstart {
					ntrees++
					continue
				}
				if ntrees >= rngstop {
					continue
				}
			}
			rt := gophy.ReadNewickString(ln)
			var t gophy.Tree
			t.Index = ntrees
			t.Instantiate(rt)
			trees = append(trees, t)
			for _, n := range t.Tips {
				if gophy.StringSliceContains(ignore, n.Nam) {
					continue
				}
				if _, ok := maptips[n.Nam]; !ok {
					maptips[n.Nam] = numtips
					mapints[numtips] = n.Nam
					numtips++
				}
			}
			ntrees++
			readtrees++
			if ntrees%10 == 0 {
				fmt.Fprint(os.Stderr, ".")
			}
			if ntrees%100 == 0 {
				fmt.Fprint(os.Stderr, "\n")
			}
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Printf("read %s bytes: %v", ln, err)
			break
		}
	}
	fmt.Fprint(os.Stderr, "\n")

	//read edges
	bps = PReadTrees(trees, *wks, rp, mapints, maptips)
	count := 0
	for i := range bps {
		bps[i].Index = count
		count++
	}
	//end read edges

	end := time.Now()
	fmt.Fprintln(os.Stderr, "trees read:", readtrees)
	fmt.Fprintln(os.Stderr, "edges skipped:", skipped)
	fmt.Fprintln(os.Stderr, "edges read:", len(bps), end.Sub(start))

	//-----------------
	// now doing things with these edges
	//---------
	// get the tree indices for the biparts
	// used for rf, rfp, and pc
	if *rf || *rfp || *rfw || *rfwp || *pc {
		bpts = make(map[int][]int, ntrees)
		for i, b := range bps {
			for _, j := range b.TreeIndices {
				bpts[j] = append(bpts[j], i) // key is tree and value is bipart
			}
		}
	}
	// output edges
	if *oed {
		fmt.Println("--edges--")
		gophy.OutputEdges(mapints, bps, ntrees, *v)
	}
	// compare to some other tree or bipart
	if len(*comp) > 0 {
		runCompare(rp, ignore, *comp, *wks, mapints, maptips, bps, readtrees, *v, *tv)
	}

	// output verbose conflict information
	if len(*oconf) > 0 {
		runConflict(*oconf, *wks, bps, mapints)
	}
	// concordance

	//pairwise comparisons
	if *pca {
		fmt.Println("--biparts compared to those in the pool of biparts")
		gophy.CompareTreeToBiparts(bps, bps, *wks, mapints, *v, *tv)
		/* This is the old comparison to a pool. I don't think we wnat this
		for j, k := range bpts { // j is tree index, k is list of biparts in bps
			comptreebps := make([]gophy.Bipart, 0)
			for _, m := range k {
				comptreebps = append(comptreebps, bps[m])
			}
			fmt.Println("tree", j, ":", len(comptreebps), "biparts from compare tree")
			start := time.Now()
			gophy.CompareTreeToBiparts(bps, comptreebps, *wks, mapints, *v)
			end := time.Now()
			fmt.Fprintln(os.Stderr, "comp done:", end.Sub(start))
		}*/
	}

	// pairwise comparison each
	if *pc {
		fmt.Println("--bipart compared to those in each tree--")
		for j := 0; j < ntrees; j++ {
			comptreebps1 := make([]gophy.Bipart, 0)
			k := bpts[j]
			for _, m := range k {
				comptreebps1 = append(comptreebps1, bps[m])
			}
			for i := 0; i < ntrees; i++ {
				if j < i {
					k = bpts[i]
					comptreebps2 := make([]gophy.Bipart, 0)
					for _, m := range k {
						comptreebps2 = append(comptreebps2, bps[m])
					}
					fmt.Println("comparing", j, i)
					start := time.Now()
					gophy.CompareTreeToBiparts(comptreebps2, comptreebps1, *wks, mapints, *v, *tv)
					end := time.Now()
					fmt.Fprintln(os.Stderr, "comp done:", end.Sub(start))
				}
			}
		}
	}

	//calculate rf

	if *rf {
		runRf(ntrees, *wks, bpts, bps)
	}

	// calculate rf (partial overlap)
	if *rfp {
		runRfp(ntrees, *wks, bpts, bps)
	}
	// calculate rfw (weighted partial overlap)
	if *rfw || *rfwp {
		if *rfwp {
			runRfw(ntrees, true, *wks, bpts, bps, *v)
		} else {
			runRfw(ntrees, false, *wks, bpts, bps, *v)
		}
	}

	//memprofile
	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.WriteHeapProfile(f)
		f.Close()
	}
}

//
//
//   END MAIN
//
//

func runRf(ntrees int, workers int, bpts map[int][]int, bps []gophy.Bipart) {
	jobs := make(chan []int, ntrees*ntrees)
	results := make(chan []int, ntrees*ntrees)
	start := time.Now()
	for w := 1; w <= workers; w++ {
		go gophy.PCalcSliceIntDifferenceInt(bpts, jobs, results)
	}

	njobs := 0
	for i := 0; i < ntrees; i++ {
		for j := 0; j < ntrees; j++ {
			if i < j {
				jobs <- []int{i, j}
				njobs++
			}
		}
	}
	close(jobs)
	for i := 0; i < njobs; i++ {
		val := <-results
		fmt.Println(val[0], val[1], ":", val[2]*2)
	}
	end := time.Now()
	fmt.Println(end.Sub(start))
}

func runRfp(ntrees int, workers int, bpts map[int][]int, bps []gophy.Bipart) {
	jobs := make(chan []int, ntrees*ntrees)
	results := make(chan []int, ntrees*ntrees)
	start := time.Now()
	for w := 1; w <= workers; w++ {
		go gophy.PCalcRFDistancesPartial(bpts, bps, jobs, results)
	}

	njobs := 0
	for i := 0; i < ntrees; i++ {
		for j := 0; j < ntrees; j++ {
			if i < j {
				jobs <- []int{i, j}
				njobs++
			}
		}
	}
	close(jobs)
	for i := 0; i < njobs; i++ {
		val := <-results
		fmt.Println(strconv.Itoa(val[0]) + " " + strconv.Itoa(val[1]) + ": " + strconv.Itoa(val[2]))
	}
	end := time.Now()
	fmt.Fprintln(os.Stderr, end.Sub(start))
}

func runRfw(ntrees int, tippenalty bool, workers int, bpts map[int][]int, bps []gophy.Bipart, verbose bool) {
	jobs := make(chan []int, ntrees*ntrees)
	results := make(chan gophy.Rfwresult, ntrees*ntrees)
	start := time.Now()
	for w := 1; w <= workers; w++ {
		go gophy.PCalcRFDistancesPartialWeighted(bpts, tippenalty, bps, jobs, results)
	}

	njobs := 0
	for i := 0; i < ntrees; i++ {
		for j := 0; j < ntrees; j++ {
			if i < j {
				jobs <- []int{i, j}
				njobs++
			}
		}
	}
	close(jobs)
	for i := 0; i < njobs; i++ {
		val := <-results
		if verbose {
			fmt.Println(strconv.Itoa(val.Tree1) + " " + strconv.Itoa(val.Tree2) + ": " + strconv.FormatFloat(val.Weight, 'f', -1, 64) +
				" " + strconv.FormatFloat(val.MaxDev, 'f', -1, 64))
		} else {
			fmt.Println(strconv.Itoa(val.Tree1) + " " + strconv.Itoa(val.Tree2) + ": " + strconv.FormatFloat(val.Weight, 'f', -1, 64))
		}
	}
	end := time.Now()
	fmt.Fprintln(os.Stderr, end.Sub(start))
}

func runConflict(outfile string, workers int, bps []gophy.Bipart, mapints map[int]string) {
	fmt.Println("--general conflict--")
	f, err := os.Create(outfile)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	jobs := make(chan []int, len(bps)*len(bps))
	results := make(chan []int, len(bps)*len(bps))
	start := time.Now()

	for w := 1; w <= workers; w++ {
		go gophy.PConflicts(bps, jobs, results)
	}

	for i := range bps {
		for j := range bps {
			if i < j {
				jobs <- []int{i, j}
			}
		}
	}
	close(jobs)
	confs := make(map[int][]int) // key is bipart and value are the conflicts
	for i := range bps {
		for j := range bps {
			if i < j {
				x := <-results
				if x[2] == 1 {
					confs[x[0]] = append(confs[x[0]], x[1])
				}
			}
		}
	}
	// printing this takes a long time. wonder how we can speed this up
	//presenting results by sorting
	count2 := 0
	maxcount2 := 20
	maxcount := 20
	//first sort the confs
	nn := map[int][]int{}
	var sortedCounts []int
	for v := range confs {
		nn[bps[v].Ct] = append(nn[bps[v].Ct], v)
	}
	for k := range nn {
		sortedCounts = append(sortedCounts, k)
	}
	sort.Sort(sort.Reverse(sort.IntSlice(sortedCounts)))
	var sortedConfs []int
	for _, m := range sortedCounts {
		for _, k := range nn[m] {
			sortedConfs = append(sortedConfs, k)
		}
	}
	//then sort the list within
	//for x, y := range confs {
	for _, x := range sortedConfs {
		y := confs[x]
		fmt.Fprint(w, bps[x].Ct, " ", bps[x].NewickWithNames(mapints)+"\n")
		n := map[int][]int{}
		var a []int
		for _, v := range y {
			n[bps[v].Ct] = append(n[bps[v].Ct], v)
		}
		for k := range n {
			a = append(a, k)
		}
		sort.Sort(sort.Reverse(sort.IntSlice(a)))
		count := 0
		for _, k := range a {
			for _, s := range n[k] {
				//s is the bps index, k is the count
				fmt.Fprint(w, "  ", bps[s].Ct, " "+bps[s].NewickWithNames(mapints)+"\n")
				if count >= maxcount {
					break
				}
				count++
			}
			if count >= maxcount {
				break
			}
		}
		if count2 > maxcount2 {
			break
		}
		count2++
	}
	end := time.Now()
	err = w.Flush()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("conf done:", end.Sub(start))
}

func runCompare(rp RunParams, ignore []string, compfile string, workers int, mapints map[int]string,
	maptips map[string]int, bps []gophy.Bipart, numtrees int, verbose bool, treeverbose bool) {
	fmt.Fprintln(os.Stderr, "--biparts compared to those in", compfile, "--")
	fc, err := os.Open(compfile)
	if err != nil {
		fmt.Println(err)
	}
	defer fc.Close()
	csc := bufio.NewScanner(fc)
	comptreebps := make([]gophy.Bipart, 0)
	/*
	   read tree and get biparts
	*/
	var t gophy.Tree
	for csc.Scan() {
		ln := csc.Text()
		if len(ln) < 2 {
			continue
		}
		rt := gophy.ReadNewickString(ln)
		t.Instantiate(rt)
		for _, n := range t.Tips {
			if _, ok := maptips[n.Nam]; !ok {
				fmt.Fprintln(os.Stderr, "need to figure out what to do when these tips don't map on the compare tree", n.Nam)
				ignore = append(ignore, n.Nam)
			}
		}
		for _, n := range t.Post {
			if len(n.Chs) > 1 && n != t.Rt {
				if rp.SupCut > 0.0 && len(n.Nam) > 0 {
					if s, err := strconv.ParseFloat(n.Nam, 32); err == nil {
						if s < rp.SupCut {
							continue
						}
					}
				}
				lt := make(map[int]bool)
				rt := make(map[int]bool)
				for _, t := range t.Tips {
					if gophy.StringSliceContains(ignore, t.Nam) == false {
						rt[maptips[t.Nam]] = true
					}
				}
				for _, t := range n.GetTips() {
					if gophy.StringSliceContains(ignore, t.Nam) == false {
						lt[maptips[t.Nam]] = true
						delete(rt, maptips[t.Nam])
					}
				}
				if len(rt) < 2 {
					continue
				}
				tbp := gophy.Bipart{Lt: lt, Rt: rt, Nds: []*gophy.Node{n}}
				comptreebps = append(comptreebps, tbp)
			}
		}
		break
	}
	fmt.Fprintln(os.Stderr, "read", len(comptreebps), "biparts from compare tree")
	start := time.Now()
	//make it a for for better memory things
	// this is just going to run it on each edge independently
	for i := range comptreebps {
		tc := []gophy.Bipart{comptreebps[i]}
		gophy.CompareTreeToBiparts(bps, tc, workers, mapints, verbose, treeverbose)
	}
	//
	end := time.Now()
	if treeverbose {
		fmt.Println("TREES WITH CONFLICT (FIRST), CONCORDANCE (SECOND), UNSUPPORTED (THIRD), PROPS AFTER")
		t.Rt.Nam = ""
		for _, n := range t.Post {
			if len(n.Chs) > 1 && n != t.Rt {
				n.Nam = n.SData["conf"]
			}
		}
		fmt.Println(t.Rt.Newick(true) + ";")
		for _, n := range t.Post {
			if len(n.Chs) > 1 && n != t.Rt {
				n.Nam = n.SData["conc"]
			}
		}
		fmt.Println(t.Rt.Newick(true) + ";")
		for _, n := range t.Post {
			if len(n.Chs) > 1 && n != t.Rt {
				n.Nam = strconv.Itoa(numtrees - (int(n.FData["conc"]) + int(n.FData["conf"])))
			}
		}
		fmt.Println(t.Rt.Newick(true) + ";")
		//
		for _, n := range t.Post {
			if len(n.Chs) > 1 && n != t.Rt {
				n.Nam = strconv.FormatFloat(n.FData["conf"]/float64(numtrees), 'f', 3, 64)
			}
		}
		fmt.Println(t.Rt.Newick(true) + ";")
		for _, n := range t.Post {
			if len(n.Chs) > 1 && n != t.Rt {
				n.Nam = strconv.FormatFloat(n.FData["conc"]/float64(numtrees), 'f', 3, 64)
			}
		}
		fmt.Println(t.Rt.Newick(true) + ";")
		for _, n := range t.Post {
			if len(n.Chs) > 1 && n != t.Rt {
				n.Nam = strconv.FormatFloat(float64(numtrees-(int(n.FData["conc"])+int(n.FData["conf"])))/float64(numtrees), 'f', 3, 64)
			}
		}
		fmt.Println(t.Rt.Newick(true) + ";")

	}
	fmt.Fprintln(os.Stderr, "comp done:", end.Sub(start))
}

// PmergeBps merge the biparts from two sets for parallel use
func PmergeBps(jobs <-chan [][]gophy.Bipart, results chan<- []gophy.Bipart) {
	for j := range jobs {
		bpst1, bpst2 := j[0], j[1]
		if j[1] != nil {
			var toadd []gophy.Bipart
			for i := range bpst2 {
				match := false
				for k := range bpst1 {
					if bpst1[k].Equals(bpst2[i]) == true {
						match = true
						bpst1[k].Ct = bpst1[k].Ct + bpst2[i].Ct
						bpst1[k].TreeIndices = append(bpst1[k].TreeIndices, bpst2[i].TreeIndices...)
						bpst1[k].Nds = append(bpst1[k].Nds, bpst2[i].Nds...)
						for ke, va := range bpst2[i].NdsM {
							bpst1[k].NdsM[ke] = va
						}
						break
					}
				}
				if match == false {
					toadd = append(toadd, bpst2[i])
				}
			}
			for _, x := range toadd {
				bpst1 = append(bpst1, x)
			}
		}
		results <- bpst1
	}
}

// PDeconstructTrees deconstruct tree for parallel use
func PDeconstructTrees(rp RunParams, maptips map[string]int, mapints map[int]string, jobs <-chan gophy.Tree, results chan<- []gophy.Bipart) {
	for t := range jobs {
		bps := make([]gophy.Bipart, 0)
		for _, n := range t.Post {
			if rp.BlCut > 0 && n.Len < rp.BlCut {
				continue
			}
			//tips, only used for rfw
			if rp.IncludeTips && len(n.Chs) == 0 {
				lt := make(map[int]bool)
				rt := make(map[int]bool)
				for _, t := range n.GetTips() {
					if gophy.StringSliceContains(rp.TIgnore, t.Nam) == false {
						lt[maptips[t.Nam]] = true
					}
				}
				if len(lt) < 1 {
					continue
				}
				nm := make(map[int]*gophy.Node)
				nm[t.Index] = n
				tbp := gophy.Bipart{Lt: lt, Rt: rt, Ct: 1, NdsM: nm}
				tbp.TreeIndices = append(tbp.TreeIndices, t.Index)
				tbp.Nds = append(tbp.Nds, n)
				bps = append(bps, tbp)
			}
			//internal nodes
			if len(n.Chs) > 1 && n != t.Rt {
				if rp.SupCut > 0.0 && len(n.Nam) > 0 {
					if s, err := strconv.ParseFloat(n.Nam, 32); err == nil {
						if s < rp.SupCut {
							continue
						}
					}
				}
				lt := make(map[int]bool)
				rt := make(map[int]bool)
				for _, t := range t.Tips {
					if gophy.StringSliceContains(rp.TIgnore, t.Nam) == false {
						rt[maptips[t.Nam]] = true
					}
				}
				for _, t := range n.GetTips() {
					if gophy.StringSliceContains(rp.TIgnore, t.Nam) == false {
						lt[maptips[t.Nam]] = true
						delete(rt, maptips[t.Nam])
					}
				}
				if len(rt) < 2 {
					continue
				}
				nm := make(map[int]*gophy.Node)
				if rp.IncludeNodeMap {
					nm[t.Index] = n
				}
				tbp := gophy.Bipart{Lt: lt, Rt: rt, Ct: 1, NdsM: nm}
				tbp.TreeIndices = append(tbp.TreeIndices, t.Index)
				tbp.Nds = append(tbp.Nds, n)
				//checks just the root case where there can be dups given how we get things
				if n.Par == t.Rt {
					index := gophy.BipartSliceContains(bps, tbp)
					if index == -1 {
						//	bpind := len(bps)
						bps = append(bps, tbp)
					}
				} else {
					bps = append(bps, tbp)
				}
			}
		}
		results <- bps
	}
}

// PReadTrees read the trees and return
func PReadTrees(trees []gophy.Tree, workers int, rp RunParams, mapints map[int]string, maptips map[string]int) (bps []gophy.Bipart) {
	//
	//parallel edge reading
	//
	jobs := make(chan gophy.Tree, len(trees))
	results := make(chan []gophy.Bipart, len(trees))
	for w := 1; w <= workers; w++ {
		go PDeconstructTrees(rp, maptips, mapints, jobs, results)
	}
	for _, i := range trees {
		jobs <- i
	}
	close(jobs)
	jobs2 := make(chan [][]gophy.Bipart, len(trees))
	results2 := make(chan []gophy.Bipart, len(trees))
	for w := 1; w < workers; w++ {
		go PmergeBps(jobs2, results2)
	}
	var x1 []gophy.Bipart
	var x2 []gophy.Bipart
	njobs := 0
	second := false
	for i := range trees {
		if x1 == nil {
			x1 = <-results
			second = false
		} else {
			x2 = <-results
			x := [][]gophy.Bipart{x1, x2}
			jobs2 <- x
			njobs++
			x1 = nil
			second = true
		}
		if i%50 == 0 {
			fmt.Fprintln(os.Stderr, i)
		}
	}
	// if there is an odd number
	if second == false {
		x := [][]gophy.Bipart{x1, nil}
		jobs2 <- x
		njobs++
	}
	// merge these now
	going := true
	x1 = nil
	x2 = nil
	njobs2 := 0
	for going {
		if njobs2 > 0 {
			njobs = njobs2
		}
		njobs2 = 0
		second = false
		for i := 0; i < njobs; i++ {
			if x1 == nil {
				x1 = <-results2
				second = false
			} else {
				x2 = <-results2
				x := [][]gophy.Bipart{x1, x2}
				jobs2 <- x
				njobs2++
				x1 = nil
				second = true
			}
			if i%10 == 0 {
				fmt.Fprint(os.Stderr, "-")
			}
		}
		if njobs == 1 && x1 != nil && second == false {
			if njobs == 1 {
				bps = x1
				break
			} else {
				x := [][]gophy.Bipart{x1, nil}
				jobs2 <- x
				njobs2++
			}
		}
		//break
		if njobs2 == 0 {
			going = false
			break
		}
	}
	fmt.Fprint(os.Stderr, "\n")
	close(jobs2)
	//
	//end parallel edge reading
	//
	return
}
