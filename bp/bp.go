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
	"time"
)

/*
 This will calculate many things things necessary for comparing
 bipartitions.
*/

func main() {
	if len(os.Args) < 2 {
		fmt.Println("bipart -fn filename")
		os.Exit(1)
	}
	wks := flag.Int("wks", 1, "how many workers?")
	cut := flag.Float64("cutoff", 0.0, "cutoff (if support is present in the trees)")
	comp := flag.String("comp", "", "compare biparts to those in this file")
	oconf := flag.String("oconf", "", "run conflict by giving an output filename")
	rf := flag.Bool("rf", false, "run rf?")
	oed := flag.Bool("oed", false, "output edges?")
	fn := flag.String("fn", "", "filename")
	flag.Parse()
	//filename things
	if *cut > 0.0 {
		fmt.Println("cutoff set to:", *cut)
	}
	if len(*fn) == 0 {
		fmt.Println("need a filename")
		os.Exit(1)
	}
	f, err := os.Open(*fn)
	if err != nil {
		fmt.Println(err)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	bps := make([]gophy.Bipart, 0)  // list of all biparts
	bpts := make(map[int][]int)     // key tree index, value bipart index list
	bpbpts := make(map[int][]int)   // key bipart index, value tree index list
	bps_counts := make(map[int]int) // key bipart index, value count
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
			if _, ok := maptips[n.Nam]; !ok {
				maptips[n.Nam] = numtips
				mapints[numtips] = n.Nam
				numtips += 1
			}
		}
		for _, n := range t.Post {
			if len(n.Chs) > 1 && n != t.Rt {
				if *cut > 0.0 && len(n.Nam) > 0 {
					if s, err := strconv.ParseFloat(n.Nam, 32); err == nil {
						if s < *cut {
							skipped += 1
							continue
						}
					}
				}
				lt := make(map[int]bool)
				rt := make(map[int]bool)
				for _, t := range t.Tips {
					rt[maptips[t.Nam]] = true
				}
				for _, t := range n.Tips() {
					lt[maptips[t.Nam]] = true
					delete(rt, maptips[t.Nam])
				}
				if len(rt) < 2 {
					continue
				}
				tbp := gophy.Bipart{lt, rt}
				//could check to merge here and if we already have one, use the one we have
				index := gophy.BipartSliceContains(bps, tbp)
				if index == -1 {
					bpind := len(bps)
					bps = append(bps, tbp)
					bps_counts[len(bps)-1] = 1
					bpts[ntrees] = append(bpts[ntrees], bpind)
					bpbpts[bpind] = append(bpbpts[bpind], ntrees)
				} else {
					if gophy.IntSliceContains(bpts[ntrees], index) == false {
						bpts[ntrees] = append(bpts[ntrees], index)
						bpbpts[index] = append(bpbpts[index], ntrees)
						bps_counts[index] += 1
					}
				}
			}
		}
		ntrees += 1
	}
	end := time.Now()
	fmt.Println("trees read:", ntrees)
	fmt.Println("edges skipped:", skipped)
	fmt.Println("edges read:", len(bps), end.Sub(start))

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
					fmt.Println("need to figure out what to do when these tips don't map on the compare tree")
				}
			}
			for _, n := range t.Post {
				if len(n.Chs) > 1 && n != t.Rt {
					lt := make(map[int]bool)
					rt := make(map[int]bool)
					for _, t := range t.Tips {
						rt[maptips[t.Nam]] = true
					}
					for _, t := range n.Tips() {
						lt[maptips[t.Nam]] = true
						delete(rt, maptips[t.Nam])
					}
					if len(rt) < 2 {
						continue
					}
					tbp := gophy.Bipart{lt, rt}
					fmt.Println(lt)
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
		for j, _ := range comptreebps {
			for i, _ := range bps {
				jobs <- []int{i, j}
			}
		}
		close(jobs)
		confs := make(map[int][]int) // key is bipart and value are the conflicts
		for range comptreebps {
			for range bps {
				x := <-results
				if x[2] == 1 {
					confs[x[1]] = append(confs[x[1]], x[0])
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
				n[bps_counts[v]] = append(n[bps_counts[v]], v)
			}
			for k := range n {
				a = append(a, k)
			}
			sort.Sort(sort.Reverse(sort.IntSlice(a)))
			count := 0
			for _, k := range a {
				for _, s := range n[k] {
					//s is the bps index, k is the count
					fmt.Print(" ", bps_counts[s], " "+bps[s].NewickWithNames(mapints)+"\n")
					fmt.Println(s, k)
					if count >= 5 {
						break
					}
					count += 1
				}
				if count >= 5 {
					break
				}
			}
		}
		end := time.Now()
		fmt.Println("comp done:", end.Sub(start))
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

		for i, _ := range bps {
			for j, _ := range bps {
				if i < j {
					jobs <- []int{i, j}
				}
			}
		}
		close(jobs)
		confs := make(map[int][]int) // key is bipart and value are the conflicts
		for i, _ := range bps {
			for j, _ := range bps {
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

		       need to add ignoring taxa if they are missing
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
					fmt.Println(val[0], val[1], ":", val[2])
				}
			}
		}
		end := time.Now()
		fmt.Println(end.Sub(start))
	}
}
