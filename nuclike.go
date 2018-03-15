package gophy

import (
	"math"

	"gonum.org/v1/gonum/floats"
)

type LikeResult struct {
	value float64
	site  int
}

/*
 This is for calculating likelihoods for nucleotides
*/

// PCalcLogLike this will calculate log like in parallel
func PCalcLogLike(t *Tree, x *DNAModel, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += CalcLogLikeOneSite(t, x, 0)
	for i := 0; i < wks; i++ {
		go CalcLogLikeWork(t, x, jobs, results)
	}
	for i := 1; i < nsites; i++ {
		jobs <- i
	}
	close(jobs)
	rr := LikeResult{}
	for i := 1; i < nsites; i++ {
		rr = <-results
		fl += rr.value
		//fl += <-results
	}
	return
}

func PCalcLike(t *Tree, x *DNAModel, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += math.Log(CalcLikeOneSite(t, x, 0))
	for i := 0; i < wks; i++ {
		go CalcLikeWork(t, x, jobs, results)
	}
	for i := 1; i < nsites; i++ {
		jobs <- i
	}
	close(jobs)
	rr := LikeResult{}
	for i := 1; i < nsites; i++ {
		rr = <-results
		fl += math.Log(rr.value)
		//fl += <-results
	}
	return
}

func PCalcLikePatterns(t *Tree, x *DNAModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += math.Log(CalcLikeOneSite(t, x, 0))
	for i := 0; i < wks; i++ {
		go CalcLikeWork(t, x, jobs, results)
	}
	for i := 1; i < nsites; i++ {
		jobs <- i
	}
	close(jobs)
	rr := LikeResult{}
	for i := 1; i < nsites; i++ {
		rr = <-results
		fl += math.Log(rr.value) * patternval[rr.site]
		//fl += <-results
	}
	return
}

func PCalcLogLikePatterns(t *Tree, x *DNAModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += CalcLogLikeOneSite(t, x, 0) * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcLogLikeWork(t, x, jobs, results)
	}
	for i := 1; i < nsites; i++ {
		jobs <- i
	}
	close(jobs)
	rr := LikeResult{}
	for i := 1; i < nsites; i++ {
		rr = <-results
		fl += (rr.value * patternval[rr.site])
	}
	return
}

// PCalcLogLikeBack a bit of a shortcut. Could do better, but walks back from the n node to the root
func PCalcLogLikeBack(t *Tree, n *Node, x *DNAModel, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += CalcLogLikeOneSiteBack(t, n, x, 0)
	for i := 0; i < wks; i++ {
		go CalcLogLikeWorkBack(t, n, x, jobs, results)
	}
	for i := 1; i < nsites; i++ {
		jobs <- i
	}
	close(jobs)
	for i := 1; i < nsites; i++ {
		fl += <-results
	}
	return
}

func PCalcLogLikeMarked(t *Tree, x *DNAModel, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	//x.EmptyPDict()
	fl += CalcLogLikeOneSiteMarked(t, x, 0)
	for i := 0; i < wks; i++ {
		go CalcLogLikeWorkMarked(t, x, jobs, results)
	}
	for i := 1; i < nsites; i++ {
		jobs <- i
	}
	close(jobs)
	for i := 1; i < nsites; i++ {
		fl += <-results
	}
	return
}

// CalcLogLikeOneSite just calculate the likelihood of one site
// probably used to populate the PDict in the DNA Model so that we can reuse the calculations
func CalcLogLikeOneSite(t *Tree, x *DNAModel, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			CalcLogLikeNode(n, x, site)
		}
		if t.Rt == n {
			for i := 0; i < 4; i++ {
				t.Rt.Data[site][i] += math.Log(x.BF[i])
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
		}
	}
	return sl
}

func CalcLikeOneSite(t *Tree, x *DNAModel, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			CalcLikeNode(n, x, site)
		}
		if t.Rt == n {
			for i := 0; i < 4; i++ {
				t.Rt.Data[site][i] *= x.BF[i]
			}
			sl = floats.Sum(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLogLikeOneSiteBack like the one above but from nb to the root only
func CalcLogLikeOneSiteBack(t *Tree, nb *Node, x *DNAModel, site int) float64 {
	sl := 0.0
	going := true
	cur := nb
	for going {
		if len(cur.Chs) > 0 {
			CalcLogLikeNode(cur, x, site)
		}
		if cur == t.Rt {
			for i := 0; i < 4; i++ {
				t.Rt.Data[site][i] += math.Log(x.BF[i])
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
			going = false
			break
		}
		cur = cur.Par
	}
	return sl
}

// CalcLogLikeOneSiteMarked this uses the marked machinery to recalculate
func CalcLogLikeOneSiteMarked(t *Tree, x *DNAModel, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			if n.Marked == true {
				CalcLogLikeNode(n, x, site)
				if n != t.Rt {
					n.Par.Marked = true
				}
			}
		}
		if t.Rt == n && n.Marked == true {
			for i := 0; i < 4; i++ {
				t.Rt.Data[site][i] += math.Log(x.BF[i])
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
		} else {
			sl = floats.LogSumExp(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLogLikeWork this is intended for a worker that will be executing this per site
func CalcLogLikeWork(t *Tree, x *DNAModel, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				CalcLogLikeNode(n, x, j)
			}
			if t.Rt == n {
				for i := 0; i < 4; i++ {
					t.Rt.Data[j][i] += math.Log(x.BF[i])
				}
				sl = floats.LogSumExp(t.Rt.Data[j])
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

func CalcLikeWork(t *Tree, x *DNAModel, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				CalcLikeNode(n, x, j)
			}
			if t.Rt == n {
				for i := 0; i < 4; i++ {
					t.Rt.Data[j][i] *= x.BF[i]
				}
				sl = floats.Sum(t.Rt.Data[j])
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

// CalcLogLikeWorkBack this is intended for a worker that will be executing this per site
func CalcLogLikeWorkBack(t *Tree, nb *Node, x *DNAModel, jobs <-chan int, results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		going := true
		cur := nb
		for going {
			if len(cur.Chs) > 0 {
				CalcLogLikeNode(cur, x, j)
			}
			if cur == t.Rt {
				for i := 0; i < 4; i++ {
					t.Rt.Data[j][i] += math.Log(x.BF[i])
				}
				sl = floats.LogSumExp(t.Rt.Data[j])
				going = false
				break
			}
			cur = cur.Par
		}
		results <- sl
	}
}

// CalcLogLikeWorkMarked this is intended to calculate only on the marked nodes back to teh root
func CalcLogLikeWorkMarked(t *Tree, x *DNAModel, jobs <-chan int, results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				if n.Marked == true {
					CalcLogLikeNode(n, x, j)
				}
			}
			if t.Rt == n && n.Marked == true {
				for i := 0; i < 4; i++ {
					t.Rt.Data[j][i] += math.Log(x.BF[i])
				}
				sl = floats.LogSumExp(t.Rt.Data[j])
			} else {
				sl = floats.LogSumExp(t.Rt.Data[j])
			}
		}
		results <- sl
	}
}

// CalcLogLikeNode calculates likelihood for node
func CalcLogLikeNode(nd *Node, model *DNAModel, site int) {
	for i := 0; i < 4; i++ {
		nd.Data[site][i] = 0.
	}
	x1 := 0.0
	x2 := []float64{0.0, 0.0, 0.0, 0.0}
	for _, c := range nd.Chs {
		P := model.GetPMap(c.Len)
		if len(c.Chs) == 0 {
			for i := 0; i < 4; i++ {
				x1 = 0.0
				for j := 0; j < 4; j++ {
					x1 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] += math.Log(x1)
			}
		} else {
			for i := 0; i < 4; i++ {
				for j := 0; j < 4; j++ {
					x2[j] = math.Log(P.At(i, j)) + c.Data[site][j]
				}
				nd.Data[site][i] += floats.LogSumExp(x2)
			}
		}
	}
}

func CalcLikeNode(nd *Node, model *DNAModel, site int) {
	for i := 0; i < 4; i++ {
		nd.Data[site][i] = 1.
	}
	x1 := 0.0
	x2 := 0.0
	for _, c := range nd.Chs {
		P := model.GetPMap(c.Len)
		if len(c.Chs) == 0 {
			for i := 0; i < 4; i++ {
				x1 = 0.0
				for j := 0; j < 4; j++ {
					x1 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] *= x1
			}
		} else {
			for i := 0; i < 4; i++ {
				x2 = 0.0
				for j := 0; j < 4; j++ {
					x2 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] *= x2 //floats.LogSumExp(x2)
			}
		}
	}
}
