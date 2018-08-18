package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"runtime/pprof"
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

func main() {
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "fasta aln filename")
	wks := flag.Int("w", 4, "number of threads")
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

	fmt.Fprintln(os.Stderr, "threads:", *wks)

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
	patterns, patternsint, gapsites, constant, uninformative := gophy.GetSitePatterns(seqs, nsites, seqnames)

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
				fmt.Println(sp.PatternStr, sp.getSitePartNameString(seqnames), "x", len(sp.Columns), " #bps:", len(sp.Biparts))
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

	for i, j := range bps {
		fmt.Println(j.NewickWithNames(namesmap), bpsc[i])
	}
	//end summary
	end := time.Now()
	fmt.Fprintln(os.Stderr, end.Sub(start))
}
