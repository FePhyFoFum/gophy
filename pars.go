package gophy

import (
	"math"
	"strconv"
)

// ParsResult gives the parsimony score (value) for a site
type ParsResult struct {
	value float64
	site  int
}

// PCalcSankParsPatterns parallel caclulation of parsimony costs with patterns
func PCalcSankParsPatterns(t *Tree, numstates int, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan ParsResult, nsites)
	for i := 0; i < wks; i++ {
		go CalcSankParsWork(t, numstates, jobs, results)
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

// CalcSankParsWork ...
func CalcSankParsWork(t *Tree, numstates int, jobs <-chan int, results chan<- ParsResult) { //results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				CalcSankParsNode(n, numstates, j)
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

// CalcSankParsNode ...
func CalcSankParsNode(nd *Node, numstates int, site int) {
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

// CalcSankParsAncStateSingleSite ....
// assumes bifurcating for now
func CalcSankParsAncStateSingleSite(t *Tree, numstates int, site int) {
	ancmaps := make(map[*Node][]int) // ancestral states are stored here
	for _, n := range t.Pre {
		if len(n.Chs) == 0 {
			continue
		}
		if t.Rt == n {
			sl := MinF(n.Data[site])
			states := make([]int, 0)
			//things are stored in n.Data[j][i]
			for i := 0; i < numstates; i++ {
				if n.Data[site][i] == sl {
					states = append(states, i)
				}
			}
			ancmaps[n] = states
		}
		//
		c1v := make(map[int]bool)
		c2v := make(map[int]bool)
		for m := range ancmaps[n] {
			v := n.Data[site][m]
			for i, j := range n.Chs[0].Data[site] {
				c1 := 0.
				if m != i {
					c1 = 1.
				}
				if len(n.Chs[0].Chs) == 0 {
					if j != 1 {
						c1 += math.MaxFloat64
					}
				} else {
					c1 += j
				}
				for k, l := range n.Chs[1].Data[site] {
					c2 := 0.
					if m != k {
						c2 = 1.
					}
					if len(n.Chs[1].Chs) == 0 {
						if l != 1 {
							c2 += math.MaxFloat64
						}
					} else {
						c2 += l
					}
					if c1+c2 == v {
						c1v[i] = true
						c2v[k] = true
					}
				}
			}
		}
		c1vv := make([]int, 0)
		for i := range c1v {
			c1vv = append(c1vv, i)
		}
		c2vv := make([]int, 0)
		for i := range c2v {
			c2vv = append(c2vv, i)
		}
		ancmaps[n.Chs[0]] = c1vv
		ancmaps[n.Chs[1]] = c2vv

		if len(n.Chs) > 0 {
			tnam := ""
			for i, j := range ancmaps[n] {
				tnam += strconv.Itoa(j)
				if i+1 != len(ancmaps[n]) {
					tnam += ","
				}
			}
			n.Nam = "[&values={" + tnam + "}]"
		}
	}
}

// EstParsBL estimate the parsimony branch lengths
func EstParsBL(t *Tree, numstates int, patternval []float64, totalsites int) {
	nsites := len(patternval)
	for _, n := range t.Post {
		n.FData["parsbl"] = 0.0
	}
	for i := 0; i < nsites; i++ {
		for _, n := range t.Pre {
			if t.Rt == n {
				minj := 0
				minv := n.Data[i][0]
				for j := 1; j < numstates; j++ {
					if n.Data[i][j] < minv {
						minj = j
						minv = n.Data[i][j]
					}
				}
				n.IData["anc"] = minj
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

				}
			}
		}
	}
	for _, n := range t.Post {
		n.Len = math.Max(10e-10, n.FData["parsbl"]/(float64(totalsites)))
		//n.Len = math.Max(0.0, n.FData["parsbl"])
	}
}
