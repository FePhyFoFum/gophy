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
)

//SitePart just a charmap but might add more later
type SitePart struct {
	CharMap    map[string][]int
	Columns    []int
	PatternStr string
	Biparts    []gophy.Bipart
}

//getSitePartNameString change sitepart to named site part
func (s SitePart) getSitePartNameString(names []string) string {
	st := "("
	for _, j := range s.CharMap {
		st += "("
		for i, m := range j {
			st += names[m]
			if i+1 < len(j) {
				st += ","
			}
		}
		st += ")"
	}
	st += ")"
	return st
}

func (s SitePart) getSitePartCount() (c int) {
	c = 0
	for _, j := range s.CharMap {
		if len(j) > 1 {
			c++
		}
	}
	return
}

func constructBiparts(s SitePart) (bps []gophy.Bipart) {
	bps = make([]gophy.Bipart, 0)
	for i, j := range s.CharMap {
		if len(j) > 1 {
			bp := gophy.Bipart{}
			bp.Lt = make(map[int]bool)
			bp.Rt = make(map[int]bool)
			for _, m := range j {
				bp.Lt[m] = true
			}
			for k, l := range s.CharMap {
				if k != i {
					for _, n := range l {
						bp.Rt[n] = true
					}
				}
			}
			eq := false
			for j := range bps {
				if bps[j].Equals(bp) {
					eq = true
					break
				}
			}
			if eq == false {
				bps = append(bps, bp)
			}
		}
	}
	return
}

func combineBiparts(bps []gophy.Bipart, namesmap map[int]string) {
	rt := gophy.Node{nil, nil, "root", map[string]string{}, map[string]float64{}, map[string]int{}, 0, 0., nil, false, 0., map[float64]bool{}, nil, nil, nil, nil}
	nodesmap := make(map[int]*gophy.Node)
	for i := range namesmap {
		nd := gophy.Node{&rt, nil, strconv.Itoa(i), map[string]string{}, map[string]float64{}, map[string]int{}, 0, 0., nil, false, 0., map[float64]bool{}, nil, nil, nil, nil}
		nodesmap[i] = &nd
		rt.Chs = append(rt.Chs, &nd)
	}
	finalsets := make([]gophy.Bipart, 0)
	c := 0
	for _, i := range bps {
		nds := make([]*gophy.Node, 0)
		for j := range i.Lt {
			nds = append(nds, nodesmap[j])
		}
		m := gophy.GetMrca(nds, &rt)
		if m == &rt && len(m.Chs) <= 3 {
			continue
		} else if m != &rt && len(m.Chs) <= 2 {
			continue
		} else {
			chs := make([]*gophy.Node, 0)
			ochs := make([]*gophy.Node, 0)
			for _, j := range m.Chs {
				t := j.GetTips()
				if gophy.NodeNamesSliceIntersects(nds, t) {
					chs = append(chs, j)
				} else {
					ochs = append(ochs, j)
				}
			}
			nd := gophy.Node{m, chs, "c" + strconv.Itoa(c), map[string]string{}, map[string]float64{}, map[string]int{}, 0, 0., nil, false, 0., map[float64]bool{}, nil, nil, nil, nil}
			nd.Chs = chs
			for _, j := range chs {
				j.Par = &nd
			}
			m.Chs = ochs
			m.Chs = append(m.Chs, &nd)
			finalsets = append(finalsets, i)
		}
		c++
	}

	//fmt.Println("finalsets:", finalsets)
	t := gophy.NewTree()
	t.Instantiate(&rt)
	for _, i := range rt.GetTips() {
		n, _ := strconv.Atoi(i.Nam)
		i.Nam = namesmap[n]
	}
	fmt.Fprintln(os.Stderr, rt.Newick(false)+";")
}

func main() {
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "fasta aln filename")
	//wks := flag.Int("w", 4, "number of threads")
	v := flag.Bool("v", false, "verbose")
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")
	flag.Parse()
	if len(os.Args) < 2 {
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

	//fmt.Fprintln(os.Stderr, "threads:", *wks)

	//read a seq file
	nsites := 0
	seqs := map[string]string{}
	seqnames := make([]string, 0)
	namesmap := make(map[int]string)
	count := 0
	for _, i := range gophy.ReadSeqsFromFile(*afn) {
		seqs[i.NM] = i.SQ
		seqnames = append(seqnames, i.NM)
		nsites = len(i.SQ)
		namesmap[count] = i.NM
		count++
	}
	bf := gophy.GetEmpiricalBaseFreqs(seqs)
	//end read a seq file

	start := time.Now()

	fmt.Println("BF", bf)
	// get the site patternas
	patterns, patternsint, gapsites, constant, uninformative, _ := gophy.GetSitePatterns(seqs, nsites, seqnames)

	fmt.Println(" patterns:", len(patternsint), " gaps:", len(gapsites),
		" constant:", len(constant), " uninformative:", len(uninformative))

	//construct siteparts
	siteparts := make([]SitePart, len(patterns))
	c := 0
	for j, m := range patterns {
		sp := SitePart{}
		sp.CharMap = make(map[string][]int)
		sp.Columns = m
		sp.PatternStr = j
		for i := range seqnames {
			sji := string(j[i])
			sp.CharMap[sji] = append(sp.CharMap[sji], i)
		}
		sp.Biparts = constructBiparts(sp)
		siteparts[c] = sp
		c++
	}
	//end construct siteparts

	//read a tree file
	if len(*tfn) > 0 {
		t := gophy.ReadTreeFromFile(*tfn)
		fmt.Fprintln(os.Stderr, "read tree:", len(t.Tips))
	}
	//end read tree file

	bps := make([]gophy.Bipart, 0)
	bpsc := make([]int, 0)
	//summary
	for _, sp := range siteparts {
		if len(sp.CharMap) == 1 {
			//fmt.Println(" ", sp.getSitePartNameString(seqnames))
		} else {
			if sp.getSitePartCount() > 1 {
				if *v {
					fmt.Println(sp.PatternStr, sp.getSitePartNameString(seqnames), "x", len(sp.Columns), " #bps:", len(sp.Biparts))
				}
				for _, j := range sp.Biparts {
					eq := false
					for ci, i := range bps {
						if j.Equals(i) {
							eq = true
							bpsc[ci] += len(sp.Columns)
							break
						}
					}
					if eq == false {
						bps = append(bps, j)
						bpsc = append(bpsc, len(sp.Columns))
					}
				}
			}
		}
	}

	fmt.Println("\n--biparts--")
	tmp := make([]int, len(bpsc))
	copy(tmp, bpsc)
	ss := gophy.NewSortedIdxSliceD(tmp...)
	sort.Sort(ss)
	sbps := make([]gophy.Bipart, 0)
	for i := len(ss.Idx) - 1; i >= 0; i-- {
		if bpsc[ss.Idx[i]] == 1 && *v {
			fmt.Println(bps[ss.Idx[i]].NewickWithNames(namesmap), bpsc[ss.Idx[i]])
		} else if bpsc[ss.Idx[i]] > 1 {
			fmt.Println(bps[ss.Idx[i]].NewickWithNames(namesmap), bpsc[ss.Idx[i]])
			sbps = append(sbps, bps[ss.Idx[i]])
		}
	}
	combineBiparts(sbps, namesmap)
	//end summary
	end := time.Now()
	fmt.Println(end.Sub(start))
}
