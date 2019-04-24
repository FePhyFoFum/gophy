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
	numsites := len(seqs[0].SQs)
	for _, j := range seqs {
		for m := 0; m < numsites; m++ {
			v, _ := strconv.Atoi(strings.Trim(j.SQs[m], " "))
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
			log.Printf("read %d bytes: %v", st, err)
			break
		}
	}
	// get the last one
	cs := MSeq{cnm, strings.ToUpper(csq), strings.Split(strings.ToUpper(csq), " ")}
	seqs = append(seqs, cs)
	return
}
