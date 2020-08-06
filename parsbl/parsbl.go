// parsbl is a tool to apply parsimony branch lengths to a tree. It reads
// multistate files and so if you have a nucleotide file, you need to give
// it the -n flag and it will output a multistate file.
//
package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"strings"

	"github.com/FePhyFoFum/gophy"
)

func createMultiStateFile(infile string) {
	f, err := os.Create(infile + ".gophy.ms")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	lg := bufio.NewWriter(f)
	tseqs := gophy.ReadSeqsFromFile(infile)
	for _, i := range tseqs {
		lg.WriteString(">" + i.NM + "\n")
		vals := make([]string, len(i.SQ))
		for j := range i.SQ {
			if string(i.SQ[j]) == "A" {
				vals[j] = "0"
			} else if string(i.SQ[j]) == "C" {
				vals[j] = "1"
			} else if string(i.SQ[j]) == "G" {
				vals[j] = "2"
			} else if string(i.SQ[j]) == "T" {
				vals[j] = "3"
			} else {
				vals[j] = "-"
			}

		}
		lg.WriteString(strings.Join(vals, " ") + "\n")
	}
	lg.Flush()
}

func calcPBL(nsites int, numstates int, mseqs []gophy.MSeq, seqs map[string][]string,
	trees []*gophy.Tree) {
	x := gophy.NewMultStateModel()
	x.NumStates = numstates
	fmt.Println(x.NumStates)
	x.SetMap()
	bf := gophy.GetEmpiricalBaseFreqsMS(mseqs, x.NumStates)
	x.SetBaseFreqs(bf)
	x.EBF = x.BF

	// get the site patternas
	patterns, patternsint, gapsites, constant, uninformative, _ :=
		gophy.GetSitePatternsMS(mseqs, x.GetCharMap(), x.GetNumStates())
	for _, t := range trees {
		patternval, _ := gophy.PreparePatternVecsMS(t, patternsint, seqs, x.GetCharMap(),
			x.GetNumStates())
		//this is necessary to get order of the patters in the patternvec since they have no order
		// this will be used with fullpattern to reconstruct the sequences
		//list of sites
		fmt.Fprintln(os.Stderr, "nsites:", nsites)
		fmt.Fprintln(os.Stderr, "patterns:", len(patterns), len(patternsint))
		fmt.Fprintln(os.Stderr, "onlygaps:", len(gapsites))
		fmt.Fprintln(os.Stderr, "constant:", len(constant))
		fmt.Fprintln(os.Stderr, "uninformative:", len(uninformative))
		// model things
		gophy.PCalcSankParsPatternsMultState(t, x.GetNumStates(), patternval, 1)
		gophy.EstParsBLMultState(t, x.GetNumStates(), patternval, nsites)
		fmt.Println(t.Rt.Newick(true) + ";")
	}

}

func main() {
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
	nuc := flag.Bool("n", false, "nucleotide data?")
	flag.Parse()
	if len(*tfn) == 0 {
		fmt.Fprintln(os.Stderr, "need a tree filename (-t)")
		os.Exit(1)
	}
	if len(*afn) == 0 {
		fmt.Fprintln(os.Stderr, "need a seq filename (-s)")
		os.Exit(1)
	}
	if *nuc {
		fmt.Fprintln(os.Stderr, "nucleotide data. will write multistate file")
		createMultiStateFile(*afn)
		*afn = *afn + ".gophy.ms"
	}

	//read a tree file
	trees := gophy.ReadTreesFromFile(*tfn)
	fmt.Fprintln(os.Stderr, len(trees), "trees read")

	//read a multistate seq file
	nsites := 0
	seqs := map[string][]string{}
	mseqs, numstates := gophy.ReadMSeqsFromFile(*afn)
	seqnames := make([]string, 0)
	for _, i := range mseqs {
		seqs[i.NM] = i.SQs
		seqnames = append(seqnames, i.NM)
		nsites = len(i.SQs)
	}

	//conduct analysis
	calcPBL(nsites, numstates, mseqs, seqs, trees)
}
