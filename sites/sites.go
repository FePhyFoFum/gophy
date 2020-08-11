// sites is a simple utility to describe things about a dataset
//
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

//TODO: need to add the functionality for gaps
func combineBiparts(bps []gophy.Bipart, namesmap map[int]string) {
	//need to add the support from the smaller ones to the larger ones and go from there
	rt := gophy.Node{Nam: "root"}
	nodesmap := make(map[int]*gophy.Node)
	for i := range namesmap {
		nd := gophy.Node{Par: &rt, Nam: strconv.Itoa(i)}
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
			nd := gophy.Node{Par: m, Chs: chs, Nam: "c" + strconv.Itoa(c)}
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
	cbp := flag.Bool("c", false, "do you want to combine biparts?")
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

	//read a seq file
	seqs := map[string][]string{}
	seqnames := make([]string, 0)
	namesmap := make(map[int]string)
	namesrevmap := make(map[string]int)
	count := 0
	mseqs, numstates := gophy.ReadMSeqsFromFile(*afn)
	for _, i := range mseqs {
		seqs[i.NM] = i.SQs
		seqnames = append(seqnames, i.NM)
		namesmap[count] = i.NM
		namesrevmap[i.NM] = count
		count++
	}
	x := gophy.NewMultStateModel()
	x.NumStates = numstates
	x.SetMap()
	bf := gophy.GetEmpiricalBaseFreqsMS(mseqs, x.NumStates)
	//end read a seq file

	start := time.Now()

	fmt.Println("BF", bf)
	// get the site patternas
	patterns, patternsint, gapsites, constant, uninformative, _ :=
		gophy.GetSitePatternsMS(mseqs, x.GetCharMap(), x.GetNumStates())

	fmt.Println(" patterns:", len(patternsint), " gaps:", len(gapsites),
		" constant:", len(constant), " uninformative:", len(uninformative))
	//construct siteparts
	siteparts := make([]SitePart, len(patterns))
	colmap := make(map[int]int)
	c := 0
	for j, m := range patterns {
		sp := SitePart{}
		sp.CharMap = make(map[string][]int)
		sp.Columns = m
		sp.PatternStr = j
		for i := range seqnames {
			sji := string(j[i])
			if sji == "-" || sji == "?" {
				continue
			}
			sp.CharMap[sji] = append(sp.CharMap[sji], i)
		}
		sp.Biparts = constructBiparts(sp)
		siteparts[c] = sp
		for _, k := range sp.Columns {
			colmap[k] = c
		}
		c++
	}
	//end construct siteparts
	var t gophy.Tree
	//read a tree file
	if len(*tfn) > 0 {
		t = *gophy.ReadTreeFromFile(*tfn)
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
							bps[ci].Ct += len(sp.Columns)
							for _, k := range sp.Columns {
								bps[ci].TreeIndices = append(bps[ci].TreeIndices, k)
							}
							break
						}
					}
					if eq == false {
						j.Ct += len(sp.Columns)
						for _, k := range sp.Columns {
							j.TreeIndices = append(j.TreeIndices, k)
						}
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
	if *cbp {
		combineBiparts(sbps, namesmap)
	}
	//compare to the tree
	if *tfn != "" {
		fmt.Println("\n--treemap--")
		for i := range t.Pre {
			n := t.Pre[i]
			if len(n.Chs) < 2 {
				continue
			}
			fmt.Println(n)
			ignore := []string{}
			lt := make(map[int]bool)
			rt := make(map[int]bool)
			for _, t := range t.Tips {
				if gophy.StringSliceContains(ignore, t.Nam) == false {
					rt[namesrevmap[t.Nam]] = true
				}
			}
			for _, t := range n.GetTips() {
				if gophy.StringSliceContains(ignore, t.Nam) == false {
					lt[namesrevmap[t.Nam]] = true
					delete(rt, namesrevmap[t.Nam])
				}
			}
			if len(rt) < 2 {
				continue
			}
			tbp := gophy.Bipart{Lt: lt, Rt: rt, Nds: []*gophy.Node{n}}
			//bps that conflict (bc0) and are concordant (bc1)
			bc0 := 0
			bc1 := 0
			//sites that conflict (bc0v) and are concordant (bc1v)
			bc0v := make(map[int]bool)
			bc1v := make(map[int]bool)
			for _, t2 := range bps {
				cc := 0
				if t2.ConflictsWith(tbp) {
					cc = 1
				}
				cw := 0
				if t2.ConcordantWith(tbp) {
					cw = 1
				}
				if cc == 1 {
					bc0 += t2.Ct
					for _, j := range t2.TreeIndices {
						bc0v[j] = true
					}
				}
				if cw == 1 {
					bc1 += t2.Ct
					for _, j := range t2.TreeIndices {
						bc1v[j] = true
					}
				}

			}
			fmt.Println(" ", bc0, bc1, len(bc0v), len(bc1v))
			n.Nam = strconv.Itoa(bc0) + "/" + strconv.Itoa(bc1) + "/" +
				strconv.Itoa(len(bc0v)) + "/" + strconv.Itoa(len(bc1v))
			for x := range bc0v {
				fmt.Println(siteparts[colmap[x]].PatternStr, siteparts[colmap[x]].Columns)
			}
		}
		fmt.Println("\n--newick--\n" + t.Rt.Newick(false) + ";")
	}
	//end summary
	end := time.Now()
	fmt.Println(end.Sub(start))
}
