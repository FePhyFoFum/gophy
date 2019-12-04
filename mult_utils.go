package gophy

import (
	"fmt"
	"os"
)

//GetSitePatternsMS return site pattens when the datatype for the alignment is a map[string]string
func GetSitePatternsMS(seqs []MSeq, charMap map[string][]int, numstates int) (patterns map[string][]int,
	patternsint map[int]float64, gapsites []int, constant []int, uninformative []int, fullpattern []int) {
	nsites := len(seqs[0].SQs)
	patterns = make(map[string][]int)
	fullpattern = make([]int, nsites)
	np := 0
	npmap := make(map[string]int)
	pamap := make(map[int]string)
	for k := 0; k < nsites; k++ {
		tp := ""
		stats := make([]int, numstates)
		for i := range stats {
			stats[i] = 0
		}
		gapcount := 0
		for _, j := range seqs {
			tp += string(j.SQs[k])
			switch c := string(j.SQs[k]); c {
			case "-":
				gapcount++
			case "N":
				gapcount++
			default:
				stats[charMap[c][0]]++ //
			}
		}
		efc := len(seqs) - gapcount
		befc := false
		for i := range stats {
			if stats[i] >= efc {
				befc = true
			}
		}
		if gapcount == len(seqs) {
			gapsites = append(gapsites, k)
			continue
		} else if befc {
			constant = append(constant, k)
		}
		twocount := 0
		for i := range stats {
			if stats[i] >= 2 {
				twocount++
			}
		}
		if twocount < 2 {
			uninformative = append(uninformative, k)
		}
		if _, ok := patterns[tp]; !ok {
			patterns[tp] = make([]int, 0)
			npmap[tp] = np
			pamap[np] = tp
			np++
		}
		patterns[tp] = append(patterns[tp], k)
		fullpattern[k] = npmap[tp]
	}
	patternsint = make(map[int]float64) // key is first site, value is the number of that one
	patternsintmap := make(map[int]int) // key is first site, value is pattern pamap int
	for m, j := range patterns {
		patternsint[j[0]] = float64(len(j))
		patternsintmap[j[0]] = npmap[m]
	}
	return
}

// PreparePatternVecsMS for tree calculations
func PreparePatternVecsMS(t *Tree, patternsint map[int]float64, seqs map[string][]string,
	charMap map[string][]int, numstates int) (patternval []float64, patternvec []int) {
	patternvec = make([]int, len(patternsint))     //which site
	patternval = make([]float64, len(patternsint)) //log of number of sites
	count := 0
	for i := range patternsint {
		patternvec[count] = i
		patternval[count] = patternsint[i]
		count++
	}
	for _, n := range t.Post {
		n.Data = make([][]float64, len(patternsint))
		n.TpConds = make([][]float64, len(patternval))
		for i := 0; i < len(patternsint); i++ {
			n.Data[i] = make([]float64, numstates)
			n.TpConds[i] = make([]float64, numstates)
			for j := 0; j < numstates; j++ {
				n.Data[i][j] = 0.0
				n.TpConds[i][j] = 0.0
			}
		}
		if len(n.Chs) == 0 {
			count := 0
			for _, i := range patternvec {
				if _, ok := charMap[string(seqs[n.Nam][i])]; !ok {
					if string(seqs[n.Nam][i]) != "-" || string(seqs[n.Nam][i]) != "N" {
						fmt.Println(n.Nam, string(seqs[n.Nam][i]))
						os.Exit(0)
					}
				}
				for _, j := range charMap[string(seqs[n.Nam][i])] {
					n.Data[count][j] = 1.0
					n.TpConds[count][j] = 1.0
				}
				count++
			}
		}
	}
	return
}
