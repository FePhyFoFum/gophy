package main

import (
	"flag"
	"fmt"
	"gophy"
	"io/ioutil"
	"os"
	"strings"
	"time"
    //"log"
    //"bufio"
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
	oconf := flag.String("oconf", "", "run conflict by giving an output filename")
	rf := flag.Bool("rf", false, "run rf?")
	oed := flag.Bool("oed", false, "output edges?")
	fn := flag.String("fn", "", "filename")
	flag.Parse()
	//filename things
	if len(*fn) == 0 {
		fmt.Println("need a filename")
		os.Exit(1)
	}
	b, err := ioutil.ReadFile(*fn)
	if err != nil {
		fmt.Println(err)
	}
	ss := string(b)
	bps := make([]gophy.Bipart, 0) // list of all biparts
	bpts := make(map[int][]int)    // key tree index, value bipart index list
	bpbpts := make(map[int][]int)  // key bipart index, value tree index list
	ntrees := 0
	numtips := 0
	maptips := make(map[string]int)
	mapints := make(map[int]string)
	start := time.Now()
	for _, ln := range strings.Split(ss, "\n") {
		if len(ln) < 2 {
			continue
		}
		rt := gophy.ReadNewickString(ln)
		var t gophy.Tree
		t.Instantiate(rt)
		//trees = append(trees, t)
		for i, n := range t.Tips {
			if _, ok := maptips[n.Nam]; !ok {
				maptips[n.Nam] = i
				mapints[i] = n.Nam
				numtips += 1
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
				//could check to merge here and if we already have one, use the one we have
				index := gophy.BipartSliceContains(bps, tbp)
				if index == -1 {
					bpind := len(bps)
					bps = append(bps, tbp)
					bpts[ntrees] = append(bpts[ntrees], bpind)
					bpbpts[bpind] = append(bpbpts[bpind], ntrees)
				} else {
					if gophy.IntSliceContains(bpts[ntrees], index) == false {
						bpts[ntrees] = append(bpts[ntrees], index)
						bpbpts[index] = append(bpbpts[index], ntrees)
					}
				}
			}
		}
		ntrees += 1
	}
	end := time.Now()
	fmt.Println("trees read:", ntrees)
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
	   checking conflicts between all the biparts
	*/
	if len(*oconf) > 0 {
        fmt.Println("--conflict--")
        /*f, err := os.Create(*oconf)
        if err != nil {
            log.Fatal(err)
        }
        defer f.Close()
        w := bufio.NewWriter(f)*/
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
                        confs[x[0]] = append(confs[x[0]],x[1])
                    }
				}
			}
		}
        // printing this takes a long time. wonder how we can speed this up
        /*
        for x, y := range confs {
            fmt.Fprint(w,bps[x].StringWithNames(mapints)+"\n")
            //fmt.Println(bps[x].StringWithNames(mapints))
            for _, n := range y {
                fmt.Fprint(w," "+bps[n].StringWithNames(mapints)+"\n")
                //fmt.Println(" ",bps[n].StringWithNames(mapints))
                //fmt.Println(" ",bps[n])
            }
        }*/
		end := time.Now()
        /*err = w.Flush()
        if err != nil{
            log.Fatal(err)
        }*/
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
