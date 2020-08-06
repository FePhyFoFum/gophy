package gophy

import (
	"math"
)

// ParsResult gives the parsimony score (value) for a site
type ParsResult struct {
	value float64
	site  int
}

//PCalcSankParsPatternsMultState parallel caclulation of fitch parsimony costs with patterns
func PCalcSankParsPatternsMultState(t *Tree, numstates int, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan ParsResult, nsites)
	for i := 0; i < wks; i++ {
		go CalcSankParsWorkMultState(t, numstates, jobs, results)
	}
	for i := 0; i < nsites; i++ {
		jobs <- i
	}
	close(jobs)
	rr := ParsResult{}
	for i := 0; i < nsites; i++ {
		rr = <-results
		fl += (rr.value * patternval[rr.site])
	}
	return
}

//CalcSankParsWorkMultState the worker for parallel jobs sent to CalcSankParsNodeMultState
func CalcSankParsWorkMultState(t *Tree, numstates int, jobs <-chan int, results chan<- ParsResult) { //results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				CalcSankParsNodeMultState(n, numstates, j)
			}
			if t.Rt == n {
				sl = math.MaxFloat64
				for i := 0; i < numstates; i++ {
					if n.Data[j][i] < sl {
						sl = n.Data[j][i]
					}
				}
			}
		}
		results <- ParsResult{value: sl, site: j}
	}
}

// CalcSankParsNodeMultState calculates Sank parsimony for a node
func CalcSankParsNodeMultState(nd *Node, numstates int, site int) {
	for i := 0; i < numstates; i++ {
		nd.Data[site][i] = 0.
	}
	for _, c := range nd.Chs {
		// tip , assume that the data are 1.0 and not data are 0.0
		if len(c.Chs) == 0 {
			for i := 0; i < numstates; i++ {
				minh := math.MaxFloat64
				for j := 0; j < numstates; j++ {
					tempv := 0.0
					if c.Data[site][j] == 0.0 {
						tempv += math.MaxFloat64
					} else {
						if i == j {
							tempv = 0.0
						} else {
							tempv = 1.0
						}
					}
					if tempv < minh {
						minh = tempv
					}
				}
				nd.Data[site][i] += minh
			}
		} else {
			for i := 0; i < numstates; i++ {
				minh := math.MaxFloat64
				for j := 0; j < numstates; j++ {
					tempv := 0.0
					if i != j {
						tempv += 1.0
					}
					tempv += c.Data[site][j]
					if tempv < minh {
						minh = tempv
					}
				}
				nd.Data[site][i] += minh
			}
		}
	}
}

// EstParsBLMultState estimates the parsimony branch lengths for multistates
//  this will just pick one when there are equivalent
func EstParsBLMultState(t *Tree, numstates int, patternval []float64, totalsites int) {
	nsites := len(patternval)
	for _, n := range t.Post {
		n.FData["parsbl"] = 0.0
	}
	for i := 0; i < nsites; i++ {
		for _, n := range t.Pre {
			if n == t.Rt {
				minj := 0
				minv := n.Data[i][0]
				for j := 1; j < numstates; j++ {
					if n.Data[i][j] < minv {
						minj = j
						minv = n.Data[i][j]
					}
				}
				n.IData["anc"] = minj
				//fmt.Println(n.Newick(false), n.IData["anc"], n.Data[i])
			} else {
				from := n.Par.IData["anc"]
				if len(n.Chs) > 0 {
					minj := math.MaxInt64
					minv := math.MaxFloat64
					for j := 0; j < numstates; j++ {
						add := 0.
						if j != from {
							add += 1.
						}
						if (n.Data[i][j] + add) < minv {
							minj = j
							minv = (n.Data[i][j] + add)
						}
					}
					if minj != from {
						n.FData["parsbl"] += 1. * patternval[i]
					}
					n.IData["anc"] = minj
				} else {
					minj := 0
					for j := 0; j < numstates; j++ {
						if n.Data[i][j] == 1.0 {
							minj = j
							break
						}
					}
					if minj != from {
						n.FData["parsbl"] += 1. * patternval[i]
					}
					//fmt.Println(n.Newick(false), n.IData["anc"], n.Data[i])
				}
			}
		}
	}
	for _, n := range t.Post {
		//prop or count
		n.Len = math.Max(10e-10, n.FData["parsbl"]/(float64(totalsites)))
		//n.Len = math.Max(0.0, n.FData["parsbl"])
	}
}
