package gophy

import (
	"bufio"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

// MSeq minimal seq struct
type MSeq struct {
	NM  string
	SQ  string
	SQs []string //SQ seperated by spaces
}

// GetEmpiricalBaseFreqsMS get the empirical base freqs from the seqs
func GetEmpiricalBaseFreqsMS(seqs []MSeq, numstates int) (bf []float64) {
	bf = make([]float64, numstates)
	statecounts := make([]int, numstates)
	total := 0
	for _, j := range seqs {
		for _, m := range j.SQs {
			v, _ := strconv.Atoi(m)
			statecounts[v]++
		}
	}
	for _, i := range statecounts {
		total += i
	}
	for i := range statecounts {
		bf[i] = float64(statecounts[i]) / float64(total)
	}
	return
}

//ReadMSeqsFromFile obvious
// file should be fasta in format with starts seperated by spaces and starting at 0 going to whatever
func ReadMSeqsFromFile(filen string) (seqs []MSeq, numstates int) {
	file, err := os.Open(filen)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	var csq string
	var cnm string
	first := true
	reader := bufio.NewReader(file)
	numstates = 0
	for {
		st, err := reader.ReadString('\n')
		if len(st) > 0 {
			if string(st[0]) == ">" {
				if first == true {
					first = false
					cnm = strings.TrimRight(st[1:], "\n")
				} else {
					cs := MSeq{cnm, strings.ToUpper(csq), strings.Split(strings.ToUpper(csq), " ")}
					for _, i := range cs.SQs {
						sa, _ := strconv.Atoi(i)
						if (sa + 1) > numstates {
							numstates = sa + 1
						}
					}
					seqs = append(seqs, cs)
					csq = ""
					cnm = strings.TrimRight(st[1:], "\n")
				}
			} else {
				csq += strings.TrimRight(st, "\n")
			}
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Printf("read %v bytes: %v", st, err)
			break
		}
	}
	// get the last one
	cs := MSeq{cnm, strings.ToUpper(csq), strings.Split(strings.ToUpper(csq), " ")}
	seqs = append(seqs, cs)
	return
}

// ReadPatternsMSeqsFromFile return the seqs and patternsint
func ReadPatternsMSeqsFromFile(sfn string) (seqs map[string][]string,
	patternsint map[int]float64, nsites int, bf []float64, numstates int) {
	nsites = 0
	seqs = map[string][]string{}
	seqnames := make([]string, 0)
	mseqs, numstates := ReadMSeqsFromFile(sfn)
	for _, i := range mseqs {
		seqs[i.NM] = i.SQs
		seqnames = append(seqnames, i.NM)
		nsites = len(i.SQs)
	}
	// get the site patternas
	bf = GetEmpiricalBaseFreqsMS(mseqs, numstates)
	//patterns, patternsint, gapsites, constant, uninformative, _ := GetSitePatterns(seqs, nsites, seqnames)
	_, patternsint, _, _, _, _ = GetSitePatternsMS(mseqs, GetMap(numstates), nsites)

	//list of sites
	//fmt.Fprintln(os.Stderr, "nsites:", nsites)
	//fmt.Fprintln(os.Stderr, "patterns:", len(patterns), len(patternsint))
	//fmt.Fprintln(os.Stderr, "onlygaps:", len(gapsites))
	//fmt.Fprintln(os.Stderr, "constant:", len(constant))
	//fmt.Fprintln(os.Stderr, "uninformative:", len(uninformative))
	return
}
