package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"runtime/pprof"
	"time"

	"github.com/FePhyFoFum/gophy"
)

type SitePart struct {
	CharMap map[string][]int
}

func (s SitePart) GetSitePartNameString(names []string) string {
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
	for _, i := range gophy.ReadSeqsFromFile(*afn) {
		seqs[i.NM] = i.SQ
		seqnames = append(seqnames, i.NM)
		nsites = len(i.SQ)
	}
	bf := gophy.GetEmpiricalBaseFreqs(seqs)

	fmt.Println("BF", bf)
	// get the site patternas
	patterns, patternsint, gapsites, constant, uninformative := gophy.GetSitePatterns(seqs, nsites, seqnames)

	fmt.Println(patterns, patternsint, gapsites, constant, uninformative)
	for j, m := range patterns {
		fmt.Println(j, m)
		sp := SitePart{}
		sp.CharMap = make(map[string][]int)
		for i, _ := range seqnames {
			sp.CharMap[string(j[i])] = append(sp.CharMap[string(j[i])], i)
			//fmt.Print(" ", k, " ", string(j[i]), "\n")
		}
		fmt.Println(" ", sp.GetSitePartNameString(seqnames))
	}
	//read a tree file
	if len(*tfn) > 0 {
		f, err := os.Open(*tfn)
		if err != nil {
			fmt.Println(err)
		}
		defer f.Close()
		scanner := bufio.NewScanner(f)
		var rt *gophy.Node
		t := gophy.NewTree()
		for scanner.Scan() {
			ln := scanner.Text()
			if len(ln) < 2 {
				continue
			}
			rt = gophy.ReadNewickString(ln)
			t.Instantiate(rt)
		}
	}
	//end read tree file

	start := time.Now()
	end := time.Now()
	fmt.Fprintln(os.Stderr, end.Sub(start))
}
