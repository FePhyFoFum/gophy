package main

import (
	"bufio"
	"flag"
	"fmt"
	"gophy"
	"log"
	"os"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"time"
)

// RunParams options for the run
type RunParams struct {
	BlCut   float64
	SupCut  float64
	TIgnore []string
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
					if gophy.SliceStringContains(rp.TIgnore, t.Nam) == false {
						rt[maptips[t.Nam]] = true
					}
				}
				for _, t := range n.GetTips() {
					if gophy.SliceStringContains(rp.TIgnore, t.Nam) == false {
						lt[maptips[t.Nam]] = true
						delete(rt, maptips[t.Nam])
					}
				}
				if len(rt) < 2 {
					continue
				}
				tbp := gophy.Bipart{Lt: lt, Rt: rt, Ct: 1}
				tbp.TreeIndices = append(tbp.TreeIndices, t.Index)
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

/*
 This will calculate many things things necessary for comparing
 bipartitions.
*/

func main() {
	if len(os.Args) < 2 {
		fmt.Println("bp -t treefile")
		os.Exit(1)
	}
	rp := RunParams{BlCut: 0, SupCut: 0, TIgnore: nil}
	wks := flag.Int("wks", 2, "how many workers?")
	cut := flag.Float64("cutoff", 0.0, "cutoff (if support is present in the trees)")
	blcut := flag.Float64("blcut", 0.0, "branch length cutoff")
	//oed := flag.Bool(p"ed", false, "output edges?")
	oconf := flag.String("oconf", "", "run conflict by giving an output filename")
	rf := flag.Bool("rf", false, "run rf?")
	rfp := flag.Bool("rfp", false, "run rf (partial overlap)?")
	comp := flag.String("comp", "", "compare biparts to those in this file")
	fn := flag.String("t", "", "tree filename")
	ig := flag.String("ig", "", "ignore these taxa (comma not space separated)")
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")
	flag.Parse()
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
	if len(*fn) == 0 {
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
	scanner := bufio.NewScanner(f)
	bps := make([]gophy.Bipart, 0) // list of all biparts
	bpts := make(map[int][]int)    // key tree index, value bipart index list
	ntrees := 0
	trees := make([]gophy.Tree, 0)
	numtips := 0
	skipped := 0
	maptips := make(map[string]int)
	mapints := make(map[int]string)
	start := time.Now()
	// reading the trees
	fmt.Fprint(os.Stderr, "reading trees\n")
	for scanner.Scan() {
		ln := scanner.Text()
		if len(ln) < 2 {
			continue
		}
		rt := gophy.ReadNewickString(ln)
		var t gophy.Tree
		t.Index = ntrees
		t.Instantiate(rt)
		trees = append(trees, t)
		for _, n := range t.Tips {
			if gophy.SliceStringContains(ignore, n.Nam) {
				continue
			}
			if _, ok := maptips[n.Nam]; !ok {
				maptips[n.Nam] = numtips
				mapints[numtips] = n.Nam
				numtips++
			}
		}
		ntrees++
		if ntrees%10 == 0 {
			fmt.Fprint(os.Stderr, ".")
		}
		if ntrees%100 == 0 {
			fmt.Fprint(os.Stderr, "\n")
		}
	}
	fmt.Fprint(os.Stderr, "\n")

	//
	//parallel edge reading
	//
	jobs := make(chan gophy.Tree, len(trees))
	results := make(chan []gophy.Bipart, len(trees))
	for w := 1; w <= *wks; w++ {
		go PDeconstructTrees(rp, maptips, mapints, jobs, results)
	}
	for _, i := range trees {
		jobs <- i
	}
	close(jobs)
	jobs2 := make(chan [][]gophy.Bipart, len(trees))
	results2 := make(chan []gophy.Bipart, len(trees))
	for w := 1; w < *wks; w++ {
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
	end := time.Now()
	//
	//end parallel edge reading
	//
	fmt.Fprintln(os.Stderr, "trees read:", ntrees)
	fmt.Fprintln(os.Stderr, "edges skipped:", skipped)
	fmt.Fprintln(os.Stderr, "edges read:", len(bps), end.Sub(start))

	if len(*comp) > 0 {
		fmt.Println("--biparts compared to those in", *comp, "--")
		fc, err := os.Open(*comp)
		if err != nil {
			fmt.Println(err)
		}
		defer fc.Close()
		csc := bufio.NewScanner(fc)
		comptreebps := make([]gophy.Bipart, 0)
		/*
		   read tree and get biparts
		*/
		for csc.Scan() {
			ln := csc.Text()
			if len(ln) < 2 {
				continue
			}
			rt := gophy.ReadNewickString(ln)
			var t gophy.Tree
			t.Instantiate(rt)
			for _, n := range t.Tips {
				if _, ok := maptips[n.Nam]; !ok {
					fmt.Println("need to figure out what to do when these tips don't map on the compare tree", n.Nam)
					ignore = append(ignore, n.Nam)
				}
			}
			for _, n := range t.Post {
				if len(n.Chs) > 1 && n != t.Rt {
					lt := make(map[int]bool)
					rt := make(map[int]bool)
					for _, t := range t.Tips {
						if gophy.SliceStringContains(ignore, t.Nam) == false {
							rt[maptips[t.Nam]] = true
						}
					}
					for _, t := range n.GetTips() {
						if gophy.SliceStringContains(ignore, t.Nam) == false {
							lt[maptips[t.Nam]] = true
							delete(rt, maptips[t.Nam])
						}
					}
					if len(rt) < 2 {
						continue
					}
					tbp := gophy.Bipart{Lt: lt, Rt: rt}
					comptreebps = append(comptreebps, tbp)
				}
			}
			break
		}
		fmt.Println("read", len(comptreebps), "biparts from compare tree")
		start := time.Now()
		gophy.CompareTreeToBiparts(bps, comptreebps, *wks, mapints)
		end := time.Now()
		fmt.Fprintln(os.Stderr, "comp done:", end.Sub(start))
	}
	if len(*oconf) > 0 {
		fmt.Println("--general conflict--")
		f, err := os.Create(*oconf)
		if err != nil {
			log.Fatal(err)
		}
		defer f.Close()
		w := bufio.NewWriter(f)
		jobs := make(chan []int, len(bps)*len(bps))
		results := make(chan []int, len(bps)*len(bps))
		start := time.Now()

		for w := 1; w <= *wks; w++ {
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
	/*
		calculate rf
	*/
	if *rf || *rfp {
		bpts = make(map[int][]int, ntrees)
		for i, b := range bps {
			for _, j := range b.TreeIndices {
				bpts[j] = append(bpts[j], i)
			}
		}
	}
	if *rf {
		jobs := make(chan []int, ntrees*ntrees)
		results := make(chan []int, ntrees*ntrees)
		start := time.Now()
		for w := 1; w <= *wks; w++ {
			go gophy.PCalcSliceIntDifferenceInt(bpts, jobs, results)
		}

		for i := 0; i < ntrees; i++ {
			for j := 0; j < ntrees; j++ {
				if i < j {
					jobs <- []int{i, j}
				}
			}
		}
		close(jobs)
		for i := 0; i < ntrees; i++ {
			for j := 0; j < ntrees; j++ {
				if i < j {
					val := <-results
					fmt.Println(val[0], val[1], ":", val[2]*2)
				}
			}
		}
		end := time.Now()
		fmt.Println(end.Sub(start))
	}
	/*
		calculate rf (partial overlap)
	*/
	if *rfp {
		jobs := make(chan []int, ntrees*ntrees)
		results := make(chan []int, ntrees*ntrees)
		start := time.Now()
		for w := 1; w <= *wks; w++ {
			go gophy.PCalcRFDistancesPartial(bpts, bps, jobs, results)
		}

		for i := 0; i < ntrees; i++ {
			for j := 0; j < ntrees; j++ {
				if i < j {
					jobs <- []int{i, j}
				}
			}
		}
		close(jobs)
		for i := 0; i < ntrees; i++ {
			for j := 0; j < ntrees; j++ {
				if i < j {
					val := <-results
					fmt.Println(strconv.Itoa(val[0]) + " " + strconv.Itoa(val[1]) + ": " + strconv.Itoa(val[2]))
				}
			}
		}
		end := time.Now()
		fmt.Fprintln(os.Stderr, end.Sub(start))
	}
}
