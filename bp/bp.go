package main

import (
	"bufio"
	"flag"
	"fmt"
	"gophy"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
	"time"
)

/*
 This will calculate many things things necessary for comparing
 bipartitions.
*/

func main() {
	if len(os.Args) < 2 {
		fmt.Println("bp -fn filename")
		os.Exit(1)
	}
	wks := flag.Int("wks", 1, "how many workers?")
	cut := flag.Float64("cutoff", 0.0, "cutoff (if support is present in the trees)")
	blcut := flag.Float64("blcut", 0.0, "branch length cutoff")
	comp := flag.String("comp", "", "compare biparts to those in this file")
	oconf := flag.String("oconf", "", "run conflict by giving an output filename")
	rf := flag.Bool("rf", false, "run rf?")
	rfp := flag.Bool("rfp", false, "run rf (partial overlap)?")
	oed := flag.Bool("ed", false, "output edges?")
	fn := flag.String("t", "", "tree filename")
	ig := flag.String("ig", "", "ignore these taxa (comma not space separated)")
	flag.Parse()
	//filename things
	if *cut > 0.0 {
		fmt.Fprintln(os.Stderr, "cutoff set to:", *cut)
	}
	if *blcut > 0.0 {
		fmt.Fprintln(os.Stderr, "blcutoff set to:", *blcut)
	}
	if len(*fn) == 0 {
		fmt.Fprintln(os.Stderr, "need a filename")
		os.Exit(1)
	}
	ignore := []string{}
	if len(*ig) > 0 {
		ignore = strings.Split(*ig, ",")
		fmt.Fprintln(os.Stderr, "ignoring:", ignore)
	}
	f, err := os.Open(*fn)
	if err != nil {
		fmt.Println(err)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	bps := make([]gophy.Bipart, 0) // list of all biparts
	bpts := make(map[int][]int)    // key tree index, value bipart index list
	bpbpts := make(map[int][]int)  // key bipart index, value tree index list
	bpsCounts := make(map[int]int) // key bipart index, value count
	ntrees := 0
	numtips := 0
	skipped := 0
	maptips := make(map[string]int)
	mapints := make(map[int]string)
	start := time.Now()
	for scanner.Scan() {
		ln := scanner.Text()
		if len(ln) < 2 {
			continue
		}
		rt := gophy.ReadNewickString(ln)
		var t gophy.Tree
		t.Instantiate(rt)
		//trees = append(trees, t)
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
		for _, n := range t.Post {
			if *blcut > 0 && n.Len < *blcut {
				skipped++
				continue
			}
			if len(n.Chs) > 1 && n != t.Rt {
				if *cut > 0.0 && len(n.Nam) > 0 {
					if s, err := strconv.ParseFloat(n.Nam, 32); err == nil {
						if s < *cut {
							skipped++
							continue
						}
					}
				}
				lt := make(map[int]bool)
				rt := make(map[int]bool)
				for _, t := range t.Tips {
					if gophy.SliceStringContains(ignore, t.Nam) == false {
						rt[maptips[t.Nam]] = true
					}
				}
				for _, t := range n.Tips() {
					if gophy.SliceStringContains(ignore, t.Nam) == false {
						lt[maptips[t.Nam]] = true
						delete(rt, maptips[t.Nam])
					}
				}
				if len(rt) < 2 {
					continue
				}
				tbp := gophy.Bipart{Lt: lt, Rt: rt}
				//could check to merge here and if we already have one, use the one we have
				index := gophy.BipartSliceContains(bps, tbp)
				if index == -1 {
					bpind := len(bps)
					bps = append(bps, tbp)
					bpsCounts[len(bps)-1]++
					bpts[ntrees] = append(bpts[ntrees], bpind)
					bpbpts[bpind] = append(bpbpts[bpind], ntrees)
				} else {
					if gophy.IntSliceContains(bpts[ntrees], index) == false {
						bpts[ntrees] = append(bpts[ntrees], index)
						bpbpts[index] = append(bpbpts[index], ntrees)
						bpsCounts[index]++
					}
				}
			}
		}
		ntrees++
	}
	end := time.Now()
	fmt.Fprintln(os.Stderr, "trees read:", ntrees)
	fmt.Fprintln(os.Stderr, "edges skipped:", skipped)
	fmt.Fprintln(os.Stderr, "edges read:", len(bps), end.Sub(start))

	/*
	   output edges
	*/
	if *oed {
		fmt.Println("--edges--")
		for i, b := range bps {
			fmt.Println(b.StringWithNames(mapints), len(bpbpts[i]), float64(len(bpbpts[i]))/float64(ntrees))
		}
	}
	/*
			   checking conflicts between biparts and those in this tree
		       TODO: add how many trees for each and maybe sort them
	*/
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
					for _, t := range n.Tips() {
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
		jobs := make(chan []int, len(bps)*len(comptreebps))
		results := make(chan []int, len(bps)*len(comptreebps))
		start := time.Now()
		for w := 1; w <= *wks; w++ {
			go gophy.PConflictsCompTree(bps, comptreebps, jobs, results)
		}
		for j := range comptreebps {
			for i := range bps {
				jobs <- []int{i, j}
			}
		}
		close(jobs)
		confs := make(map[int][]int)   // key is bipart and value are the conflicts
		allconfs := make(map[int]bool) // list of all the conflicting biparts
		for range comptreebps {
			for range bps {
				x := <-results
				if x[2] == 1 {
					confs[x[1]] = append(confs[x[1]], x[0])
					allconfs[x[0]] = true
				}
			}
		}
		/*
		 going to make a set of concordance biparts of the set of conflicting biparts
		*/
		jobs = make(chan []int, len(allconfs)*len(allconfs))
		results = make(chan []int, len(allconfs)*len(allconfs))
		start = time.Now()

		for w := 1; w <= *wks; w++ {
			go gophy.PConcordance(bps, jobs, results)
		}

		for i := range allconfs {
			for j := range allconfs {
				jobs <- []int{i, j}
			}
		}
		close(jobs)
		bpsConcCounts := make(map[int]int) // key bipart index, value number of concordant bps
		for i := range allconfs {
			for j := range allconfs {
				if i < j {
					x := <-results
					if x[2] == 1 {
						if _, ok := bpsConcCounts[x[0]]; ok {
							bpsConcCounts[x[0]]++
						} else {
							bpsConcCounts[x[0]] = 0
						}
						if _, ok := bpsConcCounts[x[1]]; ok {
							bpsConcCounts[x[1]]++
						} else {
							bpsConcCounts[x[1]] = 0
						}
					}
				}
			}
		}
		/*
		   sorting the results so that the larger bps are listed first. stop printing after a few.
		   add a sys command for listing all the results
		*/
		for x, y := range confs {
			fmt.Print(comptreebps[x].NewickWithNames(mapints) + "\n")
			n := map[int][]int{}
			var a []int
			for _, v := range y {
				//n[bpsCounts[v]] = append(n[bpsCounts[v]], v)
				n[bpsConcCounts[v]] = append(n[bpsConcCounts[v]], v)
			}
			for k := range n {
				a = append(a, k)
			}
			sort.Sort(sort.Reverse(sort.IntSlice(a)))
			count := 0
			for _, k := range a {
				for _, s := range n[k] {
					//s is the bps index, k is the count
					//fmt.Print(" ", bpsCounts[s], " "+bps[s].NewickWithNames(mapints)+"\n")
					fmt.Print(" ", bpsConcCounts[s], " "+bps[s].NewickWithNames(mapints)+"\n")
					// fmt.Println(s, k)
					if count >= 10 {
						break
					}
					count++
				}
				if count >= 5 {
					break
				}
			}
		}
		end := time.Now()
		fmt.Fprintln(os.Stderr, "comp done:", end.Sub(start))
	}
	/*
	   checking conflicts between all the biparts
	*/
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

		for x, y := range confs {
			fmt.Fprint(w, bps[x].StringWithNames(mapints)+"\n")
			//fmt.Println(bps[x].StringWithNames(mapints))
			for _, n := range y {
				fmt.Fprint(w, " "+bps[n].StringWithNames(mapints)+"\n")
				//fmt.Println(" ",bps[n].StringWithNames(mapints))
				//fmt.Println(" ",bps[n])
			}
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
