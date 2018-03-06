package gophy

import (
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

// ReadSeqsFromFile give filename seq slice
func ReadSeqsFromFile(filen string) (seqs []Seq) {
	file, err := os.Open(filen)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	buf := make([]byte, 1)
	var cc string
	var csq string
	var cnm string
	first := true
	onname := false
	for {
		n, err := file.Read(buf)
		if n > 0 {
			cc = string(buf[:n])
			if cc == ">" {
				onname = true
				if first == true {
					first = false
				} else {
					cs := Seq{cnm, strings.ToUpper(csq)}
					csq = ""
					cnm = ""
					seqs = append(seqs, cs)
				}
			} else if onname == true {
				if cc == "\n" {
					onname = false
				} else {
					cnm += cc
				}
			} else {
				if cc != "\n" {
					csq += cc
				}
			}
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Printf("read %d bytes: %v", n, err)
			break
		}
	}
	//get the last one
	cs := Seq{cnm, strings.ToUpper(csq)}
	seqs = append(seqs, cs)
	return
}
