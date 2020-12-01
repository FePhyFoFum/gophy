package gophy

import (
	"fmt"
	"os"
	"strconv"
)

//GetSitePatterns return site pattens when the datatype for the alignment is a map[string]string
func GetSitePatterns(seqs map[string]string, nsites int, seqnames []string) (patterns map[string][]int,
	patternsint map[int]float64, gapsites []int, constant []int, uninformative []int, fullpattern []int) {
	patterns = make(map[string][]int)
	fullpattern = make([]int, nsites)
	np := 0
	npmap := make(map[string]int)
	pamap := make(map[int]string)
	for k := 0; k < nsites; k++ {
		tp := ""
		As, Cs, Gs, Ts, gapcount := 0, 0, 0, 0, 0
		for _, j := range seqnames {
			tp += string(seqs[j][k])
			switch c := string(seqs[j][k]); c {
			case "A":
				As++
			case "C":
				Cs++
			case "G":
				Gs++
			case "T":
				Ts++
			case "-":
				gapcount++
			case "N":
				gapcount++
			default:
				//fmt.Println(c)
			}
		}
		efc := len(seqs) - gapcount
		if gapcount == len(seqs) {
			gapsites = append(gapsites, k)
			continue
		} else if As >= efc || Cs >= efc || Gs >= efc || Ts >= efc {
			constant = append(constant, k)
		}
		twocount := 0
		if As >= 2 {
			twocount++
		}
		if Cs >= 2 {
			twocount++
		}
		if Gs >= 2 {
			twocount++
		}
		if Ts >= 2 {
			twocount++
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
	//fmt.Println("fp", fullpattern)
	//fmt.Println("npmap", npmap)
	//fmt.Println("pamap", pamap)
	//fmt.Println("pim", patternsintmap)
	return
}

//GetSitePatternsProt return site pattens when the datatype for the alignment is a map[string]string
func GetSitePatternsProt(seqs map[string]string, nsites int, seqnames []string) (patterns map[string][]int,
	patternsint map[int]float64, gapsites []int, constant []int, uninformative []int, fullpattern []int) {
	patterns = make(map[string][]int)
	fullpattern = make([]int, nsites)
	np := 0
	npmap := make(map[string]int)
	pamap := make(map[int]string)
	for k := 0; k < nsites; k++ {
		tp := ""
		As, Rs, Ns, Ds, Cs, Qs, Es, Gs, Hs, Is, Ls, Ks, Ms, Fs, Ps, Ss, Ts, Ws, Ys, Vs, gapcount := 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
		for _, j := range seqnames {
			tp += string(seqs[j][k])
			switch c := string(seqs[j][k]); c {
			case "A":
				As++
			case "R":
				Rs++
			case "N":
				Ns++
			case "D":
				Ds++
			case "C":
				Cs++
			case "Q":
				Qs++
			case "E":
				Es++
			case "G":
				Gs++
			case "H":
				Hs++
			case "I":
				Is++
			case "L":
				Ls++
			case "K":
				Ks++
			case "M":
				Ms++
			case "F":
				Fs++
			case "P":
				Ps++
			case "S":
				Ss++
			case "T":
				Ts++
			case "W":
				Ws++
			case "Y":
				Ys++
			case "V":
				Vs++
			case "-":
				gapcount++
			case "X":
				gapcount++
			default:
				//fmt.Println(c)
			}
		}
		efc := len(seqs) - gapcount
		if gapcount == len(seqs) {
			gapsites = append(gapsites, k)
			continue
		} else if As >= efc || Rs >= efc || Ns >= efc || Ds >= efc || Cs >= efc || Qs >= efc || Es >= efc || Gs >= efc || Hs >= efc || Is >= efc || Ls >= efc || Ks >= efc || Ms >= efc || Fs >= efc || Ps >= efc || Ss >= efc || Ts >= efc || Ws >= efc || Ys >= efc || Vs >= efc {
			constant = append(constant, k)
		}
		twocount := 0
		if As >= 2 {
			twocount++
		}
		if Rs >= 2 {
			twocount++
		}
		if Ns >= 2 {
			twocount++
		}
		if Ds >= 2 {
			twocount++
		}
		if Cs >= 2 {
			twocount++
		}
		if Qs >= 2 {
			twocount++
		}
		if Es >= 2 {
			twocount++
		}
		if Gs >= 2 {
			twocount++
		}
		if Hs >= 2 {
			twocount++
		}
		if Is >= 2 {
			twocount++
		}
		if Ls >= 2 {
			twocount++
		}
		if Ks >= 2 {
			twocount++
		}
		if Ms >= 2 {
			twocount++
		}
		if Fs >= 2 {
			twocount++
		}
		if Ps >= 2 {
			twocount++
		}
		if Ss >= 2 {
			twocount++
		}
		if Ts >= 2 {
			twocount++
		}
		if Ws >= 2 {
			twocount++
		}
		if Ys >= 2 {
			twocount++
		}
		if Vs >= 2 {
			twocount++
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
	//fmt.Println("fp", fullpattern)
	//fmt.Println("npmap", npmap)
	//fmt.Println("pamap", pamap)
	//fmt.Println("pim", patternsintmap)
	return
}

// PreparePatternVecs for tree calculations
func PreparePatternVecs(t *Tree, patternsint map[int]float64, seqs map[string]string) (patternval []float64, patternvec []int) {
	patternvec = make([]int, len(patternsint))     //which site
	patternval = make([]float64, len(patternsint)) //log of number of sites
	count := 0
	for i := range patternsint {
		patternvec[count] = i
		patternval[count] = patternsint[i]
		count++
	}
	charMap := GetNucMap()
	for _, n := range t.Post {
		n.Data = make([][]float64, len(patternsint))
		n.BData = make([][]*SupFlo, len(patternsint))
		n.TpConds = make([][]float64, len(patternval))
		for i := 0; i < len(patternsint); i++ {
			n.Data[i] = []float64{0.0, 0.0, 0.0, 0.0}
			n.BData[i] = make([]*SupFlo, 4)
			for k := 0; k < 4; k++ {
				n.BData[i][k] = NewSupFlo(0.0, 0)
			}
			n.TpConds[i] = []float64{0.0, 0.0, 0.0, 0.0}
		}
		if len(n.Chs) == 0 {
			count := 0
			for _, i := range patternvec {
				if _, ok := charMap[string(seqs[n.Nam][i])]; !ok {
					if string(seqs[n.Nam][i]) != "-" && string(seqs[n.Nam][i]) != "N" {
						fmt.Println(string(seqs[n.Nam][i]))
						os.Exit(0)
					}
				}
				for _, j := range charMap[string(seqs[n.Nam][i])] {
					n.Data[count][j] = 1.0
					n.BData[count][j].SetFloat64(1.0)
					n.TpConds[count][j] = 1.0
				}
				count++
			}
		}
	}
	return
}

// PreparePatternVecsProt for tree calculations
func PreparePatternVecsProt(t *Tree, patternsint map[int]float64, seqs map[string]string) (patternval []float64, patternvec []int) {
	patternvec = make([]int, len(patternsint))     // which site
	patternval = make([]float64, len(patternsint)) // log of number of sites
	count := 0
	for i := range patternsint {
		patternvec[count] = i
		patternval[count] = patternsint[i]
		count++
	}
	charMap := GetProtMap()
	for _, n := range t.Post {
		n.Data = make([][]float64, len(patternsint))
		n.TpConds = make([][]float64, len(patternval))
		for i := 0; i < len(patternsint); i++ {
			n.Data[i] = []float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			n.TpConds[i] = []float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
		}
		if len(n.Chs) == 0 {
			count := 0
			for _, i := range patternvec {
				if _, ok := charMap[string(seqs[n.Nam][i])]; !ok {
					if string(seqs[n.Nam][i]) != "-" && string(seqs[n.Nam][i]) != "X" {
						fmt.Println(string(seqs[n.Nam][i]))
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

// GetNucMap get the int map for DNA with ambiguities
func GetNucMap() (charMap map[string][]int) {
	charMap = make(map[string][]int)
	charMap["A"] = []int{0}
	charMap["C"] = []int{1}
	charMap["G"] = []int{2}
	charMap["T"] = []int{3}
	charMap["-"] = []int{0, 1, 2, 3}
	charMap["N"] = []int{0, 1, 2, 3}
	charMap["R"] = []int{0, 2}
	charMap["Y"] = []int{1, 3}
	charMap["M"] = []int{0, 1}
	charMap["K"] = []int{2, 3}
	charMap["S"] = []int{1, 2}
	charMap["W"] = []int{0, 3}
	charMap["H"] = []int{0, 1, 3}
	charMap["B"] = []int{1, 2, 3}
	charMap["V"] = []int{0, 1, 2}
	charMap["D"] = []int{0, 2, 3}
	return
}

// GetProtMap get the int map for DNA with ambiguities
func GetProtMap() (charMap map[string][]int) {
	charMap = make(map[string][]int)
	// As, Rs, Ns, Ds, Cs, Qs, Es, Gs, Hs, Is, Ls, Ks, Ms, Fs, Ps, Ss, Ts, Ws, Ys, Vs
	charMap["A"] = []int{0}
	charMap["R"] = []int{1}
	charMap["N"] = []int{2}
	charMap["D"] = []int{3}
	charMap["C"] = []int{4}
	charMap["Q"] = []int{5}
	charMap["E"] = []int{6}
	charMap["G"] = []int{7}
	charMap["H"] = []int{8}
	charMap["I"] = []int{9}
	charMap["L"] = []int{10}
	charMap["K"] = []int{11}
	charMap["M"] = []int{12}
	charMap["F"] = []int{13}
	charMap["P"] = []int{14}
	charMap["S"] = []int{15}
	charMap["T"] = []int{16}
	charMap["W"] = []int{17}
	charMap["Y"] = []int{18}
	charMap["V"] = []int{19}
	charMap["-"] = []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}
	charMap["X"] = []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}
	charMap["B"] = []int{2, 3}
	charMap["Z"] = []int{5, 6}
	return
}

//GetMultMap based on states without the MultStateModel struct
func GetMultMap(numStates int) (charMap map[string][]int) {
	charMap = make(map[string][]int)
	charMap["-"] = make([]int, numStates)
	charMap["N"] = make([]int, numStates)
	for i := 0; i < numStates; i++ {
		charMap[strconv.Itoa(i)] = []int{i}
		charMap["-"][i] = i
		charMap["N"][i] = i
	}
	return
}

// GetRevNucMap ...
func GetRevNucMap() (charMap map[int]string) {
	charMap = make(map[int]string)
	charMap[0] = "A"
	charMap[1] = "C"
	charMap[2] = "G"
	charMap[3] = "T"
	return
}
