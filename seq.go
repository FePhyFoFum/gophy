package gophy

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"strings"
)

// Seq minimal seq struct
type Seq struct {
	NM string
	SQ string
}

// GetFasta return a fasta string
func (s Seq) GetFasta() (ret string) {
	ret = ">" + s.NM + "\n" + s.SQ + "\n"
	return
}

// GetEmpiricalBaseFreqs get the empirical base freqs from the seqs
func GetEmpiricalBaseFreqs(seqs map[string]string) (bf []float64) {
	bf = make([]float64, 4)
	AC := 0
	CC := 0
	GC := 0
	TC := 0
	total := 0
	for _, j := range seqs {
		AC += strings.Count(j, "A")
		CC += strings.Count(j, "C")
		GC += strings.Count(j, "G")
		TC += strings.Count(j, "T")
	}
	total = AC + CC + GC + TC
	bf[0] = float64(AC) / float64(total)
	bf[1] = float64(CC) / float64(total)
	bf[2] = float64(GC) / float64(total)
	bf[3] = float64(TC) / float64(total)
	return
}

// GetEmpiricalBaseFreqsProt get the empirical base freqs from the seqs
func GetEmpiricalBaseFreqsProt(seqs map[string]string) (bf []float64) {
	bf = make([]float64, 20)
	// A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V
	AC := 0
	RC := 0
	NC := 0
	DC := 0
	CC := 0
	QC := 0
	EC := 0
	GC := 0
	HC := 0
	IC := 0
	LC := 0
	KC := 0
	MC := 0
	FC := 0
	PC := 0
	SC := 0
	TC := 0
	WC := 0
	YC := 0
	VC := 0
	total := 0
	for _, j := range seqs {
		AC += strings.Count(j, "A")
		RC += strings.Count(j, "R")
		NC += strings.Count(j, "N")
		DC += strings.Count(j, "D")
		CC += strings.Count(j, "C")
		QC += strings.Count(j, "Q")
		EC += strings.Count(j, "E")
		GC += strings.Count(j, "G")
		HC += strings.Count(j, "H")
		IC += strings.Count(j, "I")
		LC += strings.Count(j, "L")
		KC += strings.Count(j, "K")
		MC += strings.Count(j, "M")
		FC += strings.Count(j, "F")
		PC += strings.Count(j, "P")
		SC += strings.Count(j, "S")
		TC += strings.Count(j, "T")
		WC += strings.Count(j, "W")
		YC += strings.Count(j, "Y")
		VC += strings.Count(j, "V")
	}
	total = AC + RC + NC + DC + CC + QC + EC + GC + HC + IC + LC + KC + MC + FC + PC + SC + TC + WC + YC + VC
	bf[0] = float64(AC) / float64(total)
	bf[1] = float64(RC) / float64(total)
	bf[2] = float64(NC) / float64(total)
	bf[3] = float64(DC) / float64(total)
	bf[4] = float64(CC) / float64(total)
	bf[5] = float64(QC) / float64(total)
	bf[6] = float64(EC) / float64(total)
	bf[7] = float64(GC) / float64(total)
	bf[8] = float64(HC) / float64(total)
	bf[9] = float64(IC) / float64(total)
	bf[10] = float64(LC) / float64(total)
	bf[11] = float64(KC) / float64(total)
	bf[12] = float64(MC) / float64(total)
	bf[13] = float64(FC) / float64(total)
	bf[14] = float64(PC) / float64(total)
	bf[15] = float64(SC) / float64(total)
	bf[16] = float64(TC) / float64(total)
	bf[17] = float64(WC) / float64(total)
	bf[18] = float64(YC) / float64(total)
	bf[19] = float64(VC) / float64(total)
	return
}

func maxv(nums ...float32) (ret float32) {
	ret = -99999.9
	for _, n := range nums {
		if n > ret {
			ret = n
		}
	}
	return
}

func maxF(F [][]float32) (ini int, inj int, ret float32) {
	ret = -99999.9
	ini = -1
	inj = -1
	for i := range F {
		for j := range F[i] {
			if F[i][j] > ret {
				ret = F[i][j]
				ini = i
				inj = j
			}
		}
	}
	return
}

// NW toy example, scores are all default
func NW(seqs []Seq, in1 int, in2 int) {
	F := make([][]float32, len(seqs[in1].SQ))
	for i := range F {
		F[i] = make([]float32, len(seqs[in2].SQ))
	}
	for i := 0; i < len(seqs[in1].SQ); i++ {
		F[i][0] = float32(-1. * i)
	}
	for i := 0; i < len(seqs[in2].SQ); i++ {
		F[0][i] = float32(-1. * i)
	}
	for i := 1; i < len(seqs[in1].SQ); i++ {
		for j := 1; j < len(seqs[in2].SQ); j++ {
			sc := float32(0.0)
			if seqs[in1].SQ[i] == seqs[in2].SQ[j] {
				sc += 1.0
			} else {
				sc -= 1.0
			}
			match := F[i-1][j-1] + sc
			insert := F[i-1][j] + -1
			del := F[i][j-1] + -1
			F[i][j] = maxv(match, insert, del)
		}
	}
	besti, bestj, bestscore := maxF(F)
	fmt.Println(besti, bestj, bestscore)
}

// PNW NW but parallel
func PNW(seqs []Seq, jobs <-chan []int, results chan<- float32) {
	for j := range jobs {
		in1, in2 := j[0], j[1]
		//fmt.Println(in1, in2)
		F := make([][]float32, len(seqs[in1].SQ))
		for i := range F {
			F[i] = make([]float32, len(seqs[in2].SQ))
		}
		for i := 0; i < len(seqs[in1].SQ); i++ {
			F[i][0] = float32(-1. * i)
		}
		for i := 0; i < len(seqs[in2].SQ); i++ {
			F[0][i] = float32(-1. * i)
		}
		for i := 1; i < len(seqs[in1].SQ); i++ {
			for j := 1; j < len(seqs[in2].SQ); j++ {
				sc := float32(0.0)
				if seqs[in1].SQ[i] == seqs[in2].SQ[j] {
					sc += 1.0
				} else {
					sc -= 1.0
				}
				match := F[i-1][j-1] + sc
				insert := F[i-1][j] + -1
				del := F[i][j-1] + -1
				F[i][j] = maxv(match, insert, del)
			}
		}
		_, _, bestscore := maxF(F)
		//fmt.Println(besti, bestj, bestscore)
		results <- bestscore
	}
}

//ReadSeqsFromFile obvious
func ReadSeqsFromFile(filen string) (seqs []Seq) {
	file, err := os.Open(filen)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	var csq string
	var cnm string
	first := true
	reader := bufio.NewReader(file)
	for {
		st, err := reader.ReadString('\n')
		if len(st) > 0 {
			if string(st[0]) == ">" {
				if first == true {
					first = false
					cnm = strings.TrimRight(st[1:], "\n")
				} else {
					cs := Seq{cnm, strings.ToUpper(csq)}
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
	cs := Seq{cnm, strings.ToUpper(csq)}
	seqs = append(seqs, cs)
	return
}

// ReadPatternsSeqsFromFile return the seqs and patternsint
func ReadPatternsSeqsFromFile(sfn string, nucleotide bool) (seqs map[string]string, patternsint map[int]float64, nsites int, bf []float64) {
	nsites = 0
	seqs = map[string]string{}
	seqnames := make([]string, 0)
	for _, i := range ReadSeqsFromFile(sfn) {
		seqs[i.NM] = i.SQ
		seqnames = append(seqnames, i.NM)
		nsites = len(i.SQ)
	}
	// get the site patternas
	if nucleotide {
		bf = GetEmpiricalBaseFreqs(seqs)
	} else {
		bf = GetEmpiricalBaseFreqsProt(seqs)
	}
	//patterns, patternsint, gapsites, constant, uninformative, _ := GetSitePatterns(seqs, nsites, seqnames)
	_, patternsint, _, _, _, _ = GetSitePatterns(seqs, nsites, seqnames)

	//list of sites
	//fmt.Fprintln(os.Stderr, "nsites:", nsites)
	//fmt.Fprintln(os.Stderr, "patterns:", len(patterns), len(patternsint))
	//fmt.Fprintln(os.Stderr, "onlygaps:", len(gapsites))
	//fmt.Fprintln(os.Stderr, "constant:", len(constant))
	//fmt.Fprintln(os.Stderr, "uninformative:", len(uninformative))
	return
}
