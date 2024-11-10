package gophy

import (
	"math"
	"sort"

	"gonum.org/v1/gonum/floats"
)

/**
 * There are several different likelihood calculators here.
 *
 * General outline:
 * PCalc* : parallel likelihood calculators.
 *   LogLike or Like or Sup
 *   Patterns
 *   Marked or Back
 *
 */

/*
 This is for calculating likelihoods for discrete states
*/

// PCalcLogLike this will calculate log like in parallel
func PCalcLogLike(t *Tree, x *DiscreteModel, nsites int, wks int) (fl float64) {
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

// PCalcLike parallel calculate likelihood
func PCalcLike(t *Tree, x *DiscreteModel, nsites int, wks int) (fl float64) {
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

// PCalcLikePatterns parallel caclulation of likelihood with patterns
func PCalcLikePatterns(t *Tree, x *DiscreteModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += math.Log(CalcLikeOneSite(t, x, 0)) * patternval[0]
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

func PCalcSupLikePatterns(t *Tree, x *DiscreteModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeSupResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += CalcSupLikeOneSite(t, x, 0).GetLn().Float64() * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcSupLikeWork(t, x, jobs, results)
	}
	for i := 1; i < nsites; i++ {
		jobs <- i
	}
	close(jobs)
	rr := LikeSupResult{}
	for i := 1; i < nsites; i++ {
		rr = <-results
		fl += rr.value.GetLn().Float64() * patternval[rr.site]
	}
	return
}

// PCalcLikePatternsGamma parallel caclulation of likelihood with patterns with gamma
func PCalcLikePatternsGamma(t *Tree, x *DiscreteModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += math.Log(CalcLikeOneSiteGamma(t, x, 0)) * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcLikeWorkGamma(t, x, jobs, results)
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

// PCalcLikePatternsMarked parallel likelihood caclulation with patterns and just update the values
func PCalcLikePatternsMarked(t *Tree, x *DiscreteModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	//x.EmptyPDict()
	fl += math.Log(CalcLikeOneSiteMarked(t, x, 0)) * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcLikeWorkMarked(t, x, jobs, results)
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

// PCalcLikePatternsMarkedGamma parallel likelihood caclulation with patterns and just update the values
func PCalcLikePatternsMarkedGamma(t *Tree, x *DiscreteModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	//x.EmptyPDict()
	fl += math.Log(CalcLikeOneSiteMarkedGamma(t, x, 0)) * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcLikeWorkMarkedGamma(t, x, jobs, results)
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

// PCalcLogLikePatterns parallel log likeliohood calculation including patterns
func PCalcLogLikePatterns(t *Tree, x *DiscreteModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	x.EmptyPLDict()
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

// PCalcLogLikePatternsGamma parallel log likeliohood calculation including patterns
func PCalcLogLikePatternsGamma(t *Tree, x *DiscreteModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	x.EmptyPLDict()
	fl += CalcLogLikeOneSiteGamma(t, x, 0) * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcLogLikeWorkGamma(t, x, jobs, results)
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
func PCalcLogLikeBack(t *Tree, n *Node, x *DiscreteModel, nsites int, wks int) (fl float64) {
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

// PCalcLogLikeMarked parallel calculation of loglike with just updating
func PCalcLogLikeMarked(t *Tree, x *DiscreteModel, nsites int, wks int) (fl float64) {
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
func CalcLogLikeOneSite(t *Tree, x *DiscreteModel, site int) float64 {
	numstates := x.NumStates
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			CalcLogLikeNode(n, x, site)
		}
		if t.Rt == n {
			for i := 0; i < numstates; i++ {
				t.Rt.Data[site][i] += math.Log(x.BF[i])
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLogLikeOneSiteGamma ...
func CalcLogLikeOneSiteGamma(t *Tree, x *DiscreteModel, site int) float64 {
	numstates := x.NumStates
	tsl := make([]float64, x.GammaNCats)
	for p, g := range x.GammaCats {
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				CalcLogLikeNodeGamma(n, x, site, g)
			}
			if t.Rt == n {
				for i := 0; i < numstates; i++ {
					t.Rt.Data[site][i] += math.Log(x.BF[i])
				}
				tsl[p] = (floats.LogSumExp(t.Rt.Data[site]) + (math.Log(1) - math.Log(float64(x.GammaNCats))))
			}
		}
	}
	return floats.LogSumExp(tsl)
}

// CalcLikeOneSite just one site
func CalcLikeOneSite(t *Tree, x *DiscreteModel, site int) float64 {
	numstates := x.NumStates
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			CalcLikeNode(n, x, site)
		}
		if t.Rt == n {
			for i := 0; i < numstates; i++ {
				t.Rt.Data[site][i] *= x.BF[i]
			}
			sl = floats.Sum(t.Rt.Data[site])
		}
	}
	return sl
}

func CalcSupLikeOneSite(t *Tree, x *DiscreteModel, site int) *SupFlo {
	numstates := x.NumStates
	sl := NewSupFlo(0.0, 0)
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			CalcSupLikeNode(n, x, site)
		}
		if t.Rt == n {
			for i := 0; i < numstates; i++ {
				t.Rt.BData[site][i].MulEqFloat(x.BF[i])
				sl.AddEq(t.Rt.BData[site][i])
			}
		}
	}
	return sl
}

// CalcLikeOneSiteGamma just one site
func CalcLikeOneSiteGamma(t *Tree, x *DiscreteModel, site int) float64 {
	numstates := x.NumStates
	sl := 0.0
	for _, g := range x.GammaCats {
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				CalcLikeNodeGamma(n, x, site, g)
			}
			if t.Rt == n {
				for i := 0; i < numstates; i++ {
					t.Rt.Data[site][i] *= x.BF[i]
				}
				sl += floats.Sum(t.Rt.Data[site]) * (1. / float64(x.GammaNCats))
			}
		}
	}
	return sl
}

// CalcLogLikeOneSiteBack like the one above but from nb to the root only
func CalcLogLikeOneSiteBack(t *Tree, nb *Node, x *DiscreteModel, site int) float64 {
	numstates := x.NumStates
	sl := 0.0
	going := true
	cur := nb
	for going {
		if len(cur.Chs) > 0 {
			CalcLogLikeNode(cur, x, site)
		}
		if cur == t.Rt {
			for i := 0; i < numstates; i++ {
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
func CalcLogLikeOneSiteMarked(t *Tree, x *DiscreteModel, site int) float64 {
	numstates := x.NumStates
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
			for i := 0; i < numstates; i++ {
				t.Rt.Data[site][i] += math.Log(x.BF[i])
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
		} else {
			sl = floats.LogSumExp(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLikeOneSiteMarked this uses the marked machinery to recalculate
func CalcLikeOneSiteMarked(t *Tree, x *DiscreteModel, site int) float64 {
	numstates := x.NumStates
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			if n.Marked == true {
				CalcLikeNode(n, x, site)
				if n != t.Rt {
					n.Par.Marked = true
				}
			}
		}
		if t.Rt == n && n.Marked == true {
			for i := 0; i < numstates; i++ {
				t.Rt.Data[site][i] *= x.BF[i]
			}
			sl = floats.Sum(t.Rt.Data[site])
		} else {
			sl = floats.Sum(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLikeOneSiteMarkedGamma this uses the marked machinery to recalculate
func CalcLikeOneSiteMarkedGamma(t *Tree, x *DiscreteModel, site int) float64 {
	numstates := x.NumStates
	sl := 0.0
	for _, g := range x.GammaCats {
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				if n.Marked == true {
					CalcLikeNodeGamma(n, x, site, g)
					if n != t.Rt {
						n.Par.Marked = true
					}
				}
			}
			if t.Rt == n && n.Marked == true {
				for i := 0; i < numstates; i++ {
					t.Rt.Data[site][i] *= x.BF[i]
				}
				sl += floats.Sum(t.Rt.Data[site]) * (1. / float64(x.GammaNCats))
			} else {
				sl += floats.Sum(t.Rt.Data[site]) * (1. / float64(x.GammaNCats))
			}
		}
	}
	return sl
}

// CalcLogLikeWork this is intended for a worker that will be executing this per site
func CalcLogLikeWork(t *Tree, x *DiscreteModel, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	numstates := x.NumStates
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				CalcLogLikeNode(n, x, j)
			}
			if t.Rt == n {
				for i := 0; i < numstates; i++ {
					t.Rt.Data[j][i] += math.Log(x.BF[i])
				}
				sl = floats.LogSumExp(t.Rt.Data[j])
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

// CalcLogLikeWorkGamma this is intended for a worker that will be executing this per site
func CalcLogLikeWorkGamma(t *Tree, x *DiscreteModel, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	numstates := x.NumStates
	for j := range jobs {
		tsl := make([]float64, x.GammaNCats)
		for p, g := range x.GammaCats {
			for _, n := range t.Post {
				if len(n.Chs) > 0 {
					CalcLogLikeNodeGamma(n, x, j, g)
				}
				if t.Rt == n {
					for i := 0; i < numstates; i++ {
						t.Rt.Data[j][i] += math.Log(x.BF[i])
					}
					tsl[p] = (floats.LogSumExp(t.Rt.Data[j]) + (math.Log(1) - math.Log(float64(x.GammaNCats))))
				}
			}
		}
		results <- LikeResult{value: floats.LogSumExp(tsl), site: j}
	}
}

// CalcLikeWork this is the worker
func CalcLikeWork(t *Tree, x *DiscreteModel, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	numstates := x.NumStates
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				CalcLikeNode(n, x, j)
			}
			if t.Rt == n {
				for i := 0; i < numstates; i++ {
					t.Rt.Data[j][i] *= x.BF[i]
				}
				sl = floats.Sum(t.Rt.Data[j])
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

// CalcLikeWork this is the worker
func CalcSupLikeWork(t *Tree, x *DiscreteModel, jobs <-chan int, results chan<- LikeSupResult) { //results chan<- float64) {
	numstates := x.NumStates
	for j := range jobs {
		sl := NewSupFlo(0.0, 0)
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				CalcSupLikeNode(n, x, j)
			}
			if t.Rt == n {
				for i := 0; i < numstates; i++ {
					t.Rt.BData[j][i].MulEqFloat(x.BF[i])
					sl.AddEq(t.Rt.BData[j][i])
				}
			}
		}
		results <- LikeSupResult{value: sl, site: j}
	}
}

// CalcLikeWorkGamma ...
func CalcLikeWorkGamma(t *Tree, x *DiscreteModel, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	numstates := x.NumStates
	for j := range jobs {
		sl := 0.0
		for _, g := range x.GammaCats {
			for _, n := range t.Post {
				if len(n.Chs) > 0 {
					CalcLikeNodeGamma(n, x, j, g)
				}
				if t.Rt == n {
					for i := 0; i < numstates; i++ {
						t.Rt.Data[j][i] *= x.BF[i]
					}
					sl += floats.Sum(t.Rt.Data[j]) * (1. / float64(x.GammaNCats))
				}
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

// CalcLogLikeWorkBack this is intended for a worker that will be executing this per site
func CalcLogLikeWorkBack(t *Tree, nb *Node, x *DiscreteModel, jobs <-chan int, results chan<- float64) {
	numstates := x.NumStates
	for j := range jobs {
		sl := 0.0
		going := true
		cur := nb
		for going {
			if len(cur.Chs) > 0 {
				CalcLogLikeNode(cur, x, j)
			}
			if cur == t.Rt {
				for i := 0; i < numstates; i++ {
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

// CalcLikeWorkMarked this is intended to calculate only on the marked nodes back to teh root
func CalcLikeWorkMarked(t *Tree, x *DiscreteModel, jobs <-chan int, results chan<- LikeResult) {
	numstates := x.NumStates
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				if n.Marked == true {
					CalcLikeNode(n, x, j)
				}
			}
			if t.Rt == n && n.Marked == true {
				for i := 0; i < numstates; i++ {
					t.Rt.Data[j][i] *= x.BF[i]
				}
				sl = floats.Sum(t.Rt.Data[j])
			} else {
				sl = floats.Sum(t.Rt.Data[j])
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

// CalcLikeWorkMarkedGamma this is intended to calculate only on the marked nodes back to teh root
func CalcLikeWorkMarkedGamma(t *Tree, x *DiscreteModel, jobs <-chan int, results chan<- LikeResult) {
	numstates := x.NumStates
	for j := range jobs {
		sl := 0.0
		for _, g := range x.GammaCats {
			for _, n := range t.Post {
				if len(n.Chs) > 0 {
					if n.Marked == true {
						CalcLikeNodeGamma(n, x, j, g)
					}
				}
				if t.Rt == n && n.Marked == true {
					for i := 0; i < numstates; i++ {
						t.Rt.Data[j][i] *= x.BF[i]
					}
					sl += floats.Sum(t.Rt.Data[j]) * (1. / float64(x.GammaNCats))
				} else {
					sl += floats.Sum(t.Rt.Data[j]) * (1. / float64(x.GammaNCats))
				}
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

// CalcLogLikeWorkMarked this is intended to calculate only on the marked nodes back to teh root
func CalcLogLikeWorkMarked(t *Tree, x *DiscreteModel, jobs <-chan int, results chan<- float64) {
	numstates := x.NumStates
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				if n.Marked == true {
					CalcLogLikeNode(n, x, j)
				}
			}
			if t.Rt == n && n.Marked == true {
				for i := 0; i < numstates; i++ {
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
func CalcLogLikeNode(nd *Node, model *DiscreteModel, site int) {
	numstates := model.NumStates
	for i := 0; i < numstates; i++ {
		nd.Data[site][i] = 0.
	}
	x1 := 0.0
	x2 := make([]float64, numstates)
	for _, c := range nd.Chs {
		if math.IsNaN(c.Len) {
			c.Len = 0.0
		}
		if len(c.Chs) == 0 {
			P := model.GetPMap(c.Len)
			for i := 0; i < numstates; i++ {
				x1 = 0.0
				for j := 0; j < numstates; j++ {
					x1 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] += math.Log(x1)
			}
		} else {
			PL := model.GetPMapLogged(c.Len)
			for i := 0; i < numstates; i++ {
				for j := 0; j < numstates; j++ {
					//x2[j] = math.Log(P.At(i, j)) + c.Data[site][j]
					x2[j] = PL.At(i, j) + c.Data[site][j]
				}
				nd.Data[site][i] += floats.LogSumExp(x2)
			}
		}
	}
}

// CalcLogLikeNodeGamma calculates likelihood for node
func CalcLogLikeNodeGamma(nd *Node, model *DiscreteModel, site int, gammav float64) {
	numstates := model.NumStates
	for i := 0; i < numstates; i++ {
		nd.Data[site][i] = 0.
	}
	x1 := 0.0
	x2 := make([]float64, numstates)
	for _, c := range nd.Chs {
		if math.IsNaN(c.Len) {
			c.Len = 0.0
		}
		if len(c.Chs) == 0 {
			P := model.GetPMap(c.Len * gammav)
			for i := 0; i < numstates; i++ {
				x1 = 0.0
				for j := 0; j < numstates; j++ {
					x1 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] += math.Log(x1)
			}
		} else {
			PL := model.GetPMapLogged(c.Len * gammav)
			for i := 0; i < numstates; i++ {
				for j := 0; j < numstates; j++ {
					//x2[j] = math.Log(P.At(i, j)) + c.Data[site][j]
					x2[j] = PL.At(i, j) + c.Data[site][j]
				}
				nd.Data[site][i] += floats.LogSumExp(x2)
			}
		}
	}
}

// CalcLikeNode calculate the likelihood of a node
func CalcLikeNode(nd *Node, model *DiscreteModel, site int) {
	numstates := model.NumStates
	for i := 0; i < numstates; i++ {
		nd.Data[site][i] = 1.
	}
	x1 := 0.0
	x2 := 0.0
	for _, c := range nd.Chs {
		P := model.GetPMap(c.Len)
		if len(c.Chs) == 0 {
			for i := 0; i < numstates; i++ {
				x1 = 0.0
				for j := 0; j < numstates; j++ {
					x1 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] *= x1
			}
		} else {
			for i := 0; i < numstates; i++ {
				x2 = 0.0
				for j := 0; j < numstates; j++ {
					x2 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] *= x2
			}
		}
	}
}

func CalcSupLikeNode(nd *Node, model *DiscreteModel, site int) {
	numstates := model.NumStates
	for i := 0; i < numstates; i++ {
		nd.BData[site][i].SetFloat64(1.0)
	}
	x1 := NewSupFlo(0.0, 0)
	x2 := NewSupFlo(0.0, 0)
	tempSup := NewSupFlo(0.0, 0)
	for _, c := range nd.Chs {
		P := model.GetPMap(c.Len)
		if len(c.Chs) == 0 {
			for i := 0; i < numstates; i++ {
				x1.SetFloat64(0.0)
				for j := 0; j < numstates; j++ {
					tempSup.SetMantExp(c.BData[site][j].GetMant()*P.At(i, j), c.BData[site][j].GetExp())
					x1.AddEq(tempSup)
				}
				nd.BData[site][i].MulEq(x1)
			}
		} else {
			for i := 0; i < numstates; i++ {
				x2.SetFloat64(0.0)
				for j := 0; j < numstates; j++ {
					tempSup.SetMantExp(c.BData[site][j].GetMant()*P.At(i, j), c.BData[site][j].GetExp())
					x2.AddEq(tempSup)
				}
				nd.BData[site][i].MulEq(x2)
			}
		}
	}
}

// CalcLikeNodeGamma calculate the likelihood of a node
func CalcLikeNodeGamma(nd *Node, model *DiscreteModel, site int, gammav float64) {
	numstates := model.NumStates
	for i := 0; i < numstates; i++ {
		nd.Data[site][i] = 1.
	}
	x1 := 0.0
	x2 := 0.0
	for _, c := range nd.Chs {
		P := model.GetPMap(c.Len * gammav) //the only gamma bit, arg
		if len(c.Chs) == 0 {
			for i := 0; i < numstates; i++ {
				x1 = 0.0
				for j := 0; j < numstates; j++ {
					x1 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] *= x1
			}
		} else {
			for i := 0; i < numstates; i++ {
				x2 = 0.0
				for j := 0; j < numstates; j++ {
					x2 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] *= x2
			}
		}
	}
}

/*
 * calculate the conditionals for ancestral calc or branch lengths

 toward tip
tpcond  X
        | | rvcond
        | ^
        | |
        v |
        | | rvtpcond
rtcond  x
 toward root
*/

// TPconditionals regular tip conditionals
func TPconditionals(x *DiscreteModel, node *Node, patternval []float64) {
	numstates := x.NumStates
	if len(node.Chs) > 0 {
		for s := range patternval {
			for j := 0; j < numstates; j++ {
				node.TpConds[s][j] = 1.
				for _, i := range node.Chs {
					node.TpConds[s][j] *= i.RtConds[s][j]
				}
			}
		}
	}
}

// RTconditionals tipconds calculated to the rt (including BL)
func RTconditionals(x *DiscreteModel, node *Node, patternval []float64) {
	numstates := x.NumStates
	p := x.GetPCalc(node.Len)
	for s := range patternval {
		for j := 0; j < numstates; j++ {
			templike := 0.0
			for k := 0; k < numstates; k++ {
				templike += p.At(j, k) * node.TpConds[s][k]
			}
			node.RtConds[s][j] = templike
		}
	}
}

// RVconditionals take par RvTpConds and put get bl
func RVconditionals(x *DiscreteModel, node *Node, patternval []float64) {
	numstates := x.NumStates
	p := x.GetPCalc(node.Par.Len)
	for s := range patternval {
		for j := 0; j < numstates; j++ {
			node.Par.RvTpConds[s][j] = 0.0
			for k := 0; k < numstates; k++ {
				node.Par.RvTpConds[s][j] += p.At(j, k) * node.Par.RvConds[s][k]
			}
		}
	}
}

// RVTPconditionals ...
func RVTPconditionals(x *DiscreteModel, node *Node, patternval []float64) {
	numstates := x.NumStates
	for s := range patternval {
		for j := 0; j < numstates; j++ {
			node.RvConds[s][j] = node.Par.RvTpConds[s][j]
		}
		for _, oc := range node.Par.Chs {
			if node == oc {
				continue
			}
			for j := 0; j < numstates; j++ {
				node.RvConds[s][j] *= oc.RtConds[s][j]
			}
		}
	}
}

// CalcLikeFrontBack ...
func CalcLikeFrontBack(x *DiscreteModel, tree *Tree, patternval []float64) {
	numstates := x.NumStates
	for _, n := range tree.Post {
		if len(n.Chs) != 0 {
			n.TpConds = make([][]float64, len(patternval))
		}
		n.RvTpConds = make([][]float64, len(patternval))
		n.RvConds = make([][]float64, len(patternval))
		n.RtConds = make([][]float64, len(patternval))
		for i := 0; i < len(patternval); i++ {
			if len(n.Chs) != 0 {
				n.TpConds[i] = make([]float64, numstates)
				for j := 0; j < numstates; j++ {
					n.TpConds[i][j] = 1.0
				}
			}
			n.RvTpConds[i] = make([]float64, numstates)
			n.RvConds[i] = make([]float64, numstates)
			n.RtConds[i] = make([]float64, numstates)
			for j := 0; j < numstates; j++ {
				n.RvTpConds[i][j] = 1.0
				n.RvConds[i][j] = 1.0
				n.RtConds[i][j] = 1.0
			}
		}
	}
	//loglike := 0.
	for _, c := range tree.Post {
		//calculate the tip conditionals
		TPconditionals(x, c, patternval)
		//take the tip cond to the rt
		RTconditionals(x, c, patternval) // calculate from tpcond to rtcond
		/*if c == tree.Rt { // turn on if you want likelihoods
			for s := range patternval {
				tempretlike := 0.
				for i := 0; i < 4; i++ {
					tempretlike += (c.TpConds[s][i] * x.BF[i])
				}
				//fmt.Println("site", s, "log(L):", math.Log(tempretlike), "like:", tempretlike, "pattern:", patternval[s])
				//loglike -= math.Log(math.Pow(tempretlike, patternval[s]))
			}
		}*/
	}
	//fmt.Println(loglike)
	// prepare the rvcond
	for _, c := range tree.Pre {
		if c != tree.Rt { //need to set the root at 1.0s
			RVconditionals(x, c, patternval)
			RVTPconditionals(x, c, patternval)
		}
	}
}

// CalcAncStates for each node based on the calculations above
func CalcAncStates(x *DiscreteModel, tree *Tree, patternval []float64) (retstates map[*Node][][]float64) {
	numstates := x.NumStates
	CalcLikeFrontBack(x, tree, patternval)
	retstates = make(map[*Node][][]float64)
	// initialize the data storage for return
	for _, c := range tree.Pre {
		if len(c.Chs) == 0 {
			continue
		}
		ndata := make([][]float64, len(patternval))
		retstates[c] = ndata
	}
	// start reconstruction
	for i := 0; i < len(patternval); i++ {
		//for _, c := range tree.Tips {
		//	fmt.Println(c.Nam, c.TpConds[i])
		//}
		//fmt.Println("p", i, patternval[i])
		for _, c := range tree.Pre {
			if len(c.Chs) == 0 {
				continue
			}
			//fmt.Println(c.Newick(true))
			retstates[c][i] = make([]float64, numstates)
			if c == tree.Rt {
				su := 0.
				for j, s := range c.RtConds[i] {
					su += (s * x.BF[j])
				}
				for j, s := range c.RtConds[i] {
					//fmt.Print((s*x.BF[j])/su, " ")
					retstates[c][i][j] = (s * x.BF[j]) / su
				}
				//fmt.Print("\n")
			} else {
				p := x.GetPCalc(c.Len)
				//need subtree 1
				s1probs := c.TpConds
				//need subtree 2
				s2probs := c.RvConds
				tv := make([]float64, numstates)
				for j := 0; j < numstates; j++ {
					for k := 0; k < numstates; k++ {
						tv[j] += (s1probs[i][j] * p.At(j, k) * s2probs[i][k])
					}
					tv[j] *= x.BF[j]
				}
				su := 0.
				for _, s := range tv {
					su += s
				}
				for j, s := range tv {
					//fmt.Print(s/su, " ")
					retstates[c][i][j] = s / su
				}
				//fmt.Print("\n")
			}
		}
	}
	return
}

type JointConfig struct {
	Score  float64
	Config map[*Node]int
}

func Max(slice []float64) float64 {
	if len(slice) == 0 {
		panic("empty slice") // or handle it as you prefer
	}

	maxVal := slice[0]
	for _, value := range slice[1:] {
		if value > maxVal {
			maxVal = value
		}
	}
	return maxVal
}

func CalcJointAncStates(x *DiscreteModel, tree *Tree, patternval []float64) (rootconfigs [][]JointConfig) {
	numstates := x.NumStates
	//CalcLikeFrontBack(x, tree, patternval)
	rootconfigs = make([][]JointConfig, 0, 0)
	// initialize the data storage for return
	for _, c := range tree.Pre {
		if len(c.Chs) == 0 {
			continue
		}
	}
	// start reconstruction
	for i := 0; i < len(patternval); i++ { //len(patternval); i++ {
		curstates := make(map[*Node][]int)
		jconfigs := make(map[*Node][]JointConfig)
		//for _, c := range tree.Tips {
		//	fmt.Println(c.Nam, c.TpConds[i])
		//}
		//fmt.Println("p", i, patternval[i])
		for _, c := range tree.Post {
			if len(c.Chs) == 0 {
				for j, k := range c.Data[i] {
					if k == 1 {
						curstates[c] = append(curstates[c], j)
						jc := JointConfig{Score: 1.0, Config: make(map[*Node]int)}
						jc.Config[c] = j
						jconfigs[c] = append(jconfigs[c], jc)
					}
				}
				//fmt.Println(c.Data[i], curstates[c])
				continue
			}
			//assume you are calaculating log like marginals
			maxv := Max(c.Data[i])
			//insert skipping things, use c.data
			for j := 0; j < numstates; j++ {
				if maxv-c.Data[i][j] < 2 {
					curstates[c] = append(curstates[c], j)
				}
			}
			//fmt.Println(c, c.Data[i], curstates[c])
			//end skipping things
			for _, k := range curstates[c] {
				//bifurcating
				jchsc0 := []JointConfig{}
				jchsc1 := []JointConfig{}
				for c1, cc := range c.Chs {
					p := x.GetPCalc(cc.Len)
					for _, l := range jconfigs[cc] {
						jc := JointConfig{Score: 1.0, Config: make(map[*Node]int)}
						jc.Score = l.Score * p.At(k, l.Config[cc])
						for x, y := range l.Config {
							jc.Config[x] = y
						}
						if c1 == 0 {
							jchsc0 = append(jchsc0, jc)
						} else {
							jchsc1 = append(jchsc1, jc)
						}
					}
				}
				for _, l := range jchsc0 {
					for _, m := range jchsc1 {
						//fmt.Println(" l", l)
						//fmt.Println(" m", m)
						tc := JointConfig{Score: 1.0, Config: make(map[*Node]int)}
						if c == tree.Rt {
							tc.Score = l.Score * m.Score * x.BF[k]
						} else {
							tc.Score = l.Score * m.Score
						}
						tc.Config[c] = k
						for n, mm := range l.Config {
							tc.Config[n] = mm
						}
						for n, mm := range m.Config {
							tc.Config[n] = mm
						}
						//fmt.Println(" tc", tc.Config)
						jconfigs[c] = append(jconfigs[c], tc)
					}
				}
			}
			if c == tree.Rt {
				//get them in the right order
				//only include < 2lnL
				sort.Slice(jconfigs[tree.Rt], func(x, y int) bool {
					return jconfigs[tree.Rt][x].Score > jconfigs[tree.Rt][y].Score // Descending order
				})
				//fmt.Println(jconfigs[tree.Rt])
				maxV := math.Log(jconfigs[tree.Rt][0].Score)
				lastV := 0
				for j := range jconfigs[tree.Rt] {
					jconfigs[tree.Rt][j].Score = math.Log(jconfigs[tree.Rt][j].Score)
					//fmt.Println(jconfigs[tree.Rt][j].Score)
					if maxV-jconfigs[tree.Rt][j].Score < 2 {
						lastV = j
					}
				}
				lastV += 1
				rootconfigs = append(rootconfigs, jconfigs[tree.Rt][0:lastV])
			}
		}
	}
	return
}

// CalcJointAncStatesGamma same as above but with gamma, slower, given gamma vals
func CalcJointAncStatesGamma(x *DiscreteModel, tree *Tree, patternval []float64) (rootconfigs [][]JointConfig) {
	numstates := x.NumStates
	//CalcLikeFrontBack(x, tree, patternval)
	rootconfigs = make([][]JointConfig, 0, 0)
	// initialize the data storage for return
	for _, c := range tree.Pre {
		if len(c.Chs) == 0 {
			continue
		}
	}
	// start reconstruction
	for i := 0; i < len(patternval); i++ { //len(patternval); i++ {
		curstates := make(map[*Node][]int)
		jconfigs := make(map[*Node][]JointConfig)
		//for _, c := range tree.Tips {
		//	fmt.Println(c.Nam, c.TpConds[i])
		//}
		//fmt.Println("p", i, patternval[i])
		//THIS IS THE GAMMA BIT, NEED TO FIGURE OUT WHERE TO PUT IT
		//tsl := make([]float64, x.GammaNCats)
		//for p, g := range x.GammaCats {

		for _, c := range tree.Post {
			if len(c.Chs) == 0 {
				for j, k := range c.Data[i] {
					if k == 1 {
						curstates[c] = append(curstates[c], j)
						jc := JointConfig{Score: 1.0, Config: make(map[*Node]int)}
						jc.Config[c] = j
						jconfigs[c] = append(jconfigs[c], jc)
					}
				}
				//fmt.Println(c.Data[i], curstates[c])
				continue
			}
			//assume you are calaculating log like marginals
			maxv := Max(c.Data[i])
			//insert skipping things, use c.data
			for j := 0; j < numstates; j++ {
				if maxv-c.Data[i][j] < 2 {
					curstates[c] = append(curstates[c], j)
				}
			}
			//fmt.Println(c, c.Data[i], curstates[c])
			//end skipping things
			for _, k := range curstates[c] {
				//bifurcating
				jchsc0 := []JointConfig{}
				jchsc1 := []JointConfig{}
				for c1, cc := range c.Chs {
					p := x.GetPCalc(cc.Len)
					for _, l := range jconfigs[cc] {
						jc := JointConfig{Score: 1.0, Config: make(map[*Node]int)}
						jc.Score = l.Score * p.At(k, l.Config[cc])
						for x, y := range l.Config {
							jc.Config[x] = y
						}
						if c1 == 0 {
							jchsc0 = append(jchsc0, jc)
						} else {
							jchsc1 = append(jchsc1, jc)
						}
					}
				}
				for _, l := range jchsc0 {
					for _, m := range jchsc1 {
						//fmt.Println(" l", l)
						//fmt.Println(" m", m)
						tc := JointConfig{Score: 1.0, Config: make(map[*Node]int)}
						if c == tree.Rt {
							tc.Score = l.Score * m.Score * x.BF[k]
						} else {
							tc.Score = l.Score * m.Score
						}
						tc.Config[c] = k
						for n, mm := range l.Config {
							tc.Config[n] = mm
						}
						for n, mm := range m.Config {
							tc.Config[n] = mm
						}
						//fmt.Println(" tc", tc.Config)
						jconfigs[c] = append(jconfigs[c], tc)
					}
				}
			}
			if c == tree.Rt {
				//fmt.Println(jconfigs[tree.Rt])
				rootconfigs = append(rootconfigs, jconfigs[tree.Rt])
			}
		}
	}
	return
}

// CalcStochMap this has been taken from the multistate code so make sure it works for others
func CalcStochMap(x *DiscreteModel, tree *Tree, patternval []float64, time bool, from int, to int) (retstates map[*Node][][]float64) {
	CalcLikeFrontBack(x, tree, patternval)
	retstates = make(map[*Node][][]float64)
	// initialize the data storage for return
	for _, c := range tree.Pre {
		if len(c.Chs) == 0 {
			//	continue
		}
		ndata := make([][]float64, len(patternval))
		retstates[c] = ndata
	}
	// start reconstruction
	for i := 0; i < len(patternval); i++ {
		for _, c := range tree.Pre {
			retstates[c][i] = make([]float64, x.GetNumStates())
			if c == tree.Rt {
				continue
			} else {
				numM, ratM := x.GetStochMapMatrices(c.Len, from, to) //x.GetPCalc(c.Len)
				//need subtree 1
				s1probs := c.TpConds
				//need subtree 2
				s2probs := c.RvConds
				tvN := make([]float64, x.GetNumStates())
				tvR := make([]float64, x.GetNumStates())
				for j := 0; j < x.GetNumStates(); j++ {
					for k := 0; k < x.GetNumStates(); k++ {
						tvN[j] += (s1probs[i][j] * numM.At(j, k) * s2probs[i][k])
						tvR[j] += (s1probs[i][j] * ratM.At(j, k) * s2probs[i][k])
					}
					tvN[j] *= x.GetBF()[j]
					tvR[j] *= x.GetBF()[j]
				}
				if time { //time
					retstates[c][i] = tvR
				} else { //number
					retstates[c][i] = tvN
				}
			}
		}
	}
	return
}

// Subclade calculators for log and likelihoods
// these are for just calculating for one clade

// calcLogLikeOneSiteSubClade calc for just a clade, starting at a node
func calcLogLikeOneSiteSubClade(t *Tree, inn *Node, excl bool, x *DiscreteModel, site int) float64 {
	numstates := x.NumStates
	sl := 0.0
	arr := []*Node{}
	if excl == true {
		arr = t.Rt.PostorderArrayExcl(inn)
	} else {
		arr = inn.PostorderArray()
	}
	for _, n := range arr {
		if len(n.Chs) > 0 {
			CalcLogLikeNode(n, x, site)
		}
		var tn *Node
		if excl == true { // calc at rt
			tn = t.Rt
		} else {
			tn = inn
		}
		if tn == n {
			if tn == t.Rt { //only happens at the root
				for i := 0; i < numstates; i++ {
					n.Data[site][i] += math.Log(x.BF[i])
				}
				sl = floats.LogSumExp(n.Data[site])
			} else {
				//needs to get the branch length incorporated
				p := x.GetPMapLogged(n.Len)
				rtconds := make([]float64, x.GetNumStates())
				x2 := make([]float64, x.GetNumStates())
				for m := 0; m < numstates; m++ {
					for k := 0; k < numstates; k++ {
						x2[k] = p.At(m, k) + n.Data[site][k]
					}
					rtconds[m] = floats.LogSumExp(x2)
				}
				sl = floats.LogSumExp(rtconds)
			}
		}
	}
	return sl
}

func calcLogLikeOneSiteSubCladeGamma(t *Tree, inn *Node, excl bool, x *DiscreteModel, site int) float64 {
	numstates := x.NumStates
	arr := []*Node{}
	if excl == true {
		arr = t.Rt.PostorderArrayExcl(inn)
	} else {
		arr = inn.PostorderArray()
	}
	tsl := make([]float64, x.GammaNCats)
	for p, g := range x.GammaCats {
		for _, n := range arr {
			if len(n.Chs) > 0 {
				CalcLogLikeNodeGamma(n, x, site, g)
			}
			var tn *Node
			if excl == true { // calc at rt
				tn = t.Rt
			} else {
				tn = inn
			}
			if tn == n {
				if tn == t.Rt { //only happens at the root
					for i := 0; i < numstates; i++ {
						n.Data[site][i] += math.Log(x.BF[i])
					}
					tsl[p] = (floats.LogSumExp(n.Data[site]) + (math.Log(1) - math.Log(float64(x.GammaNCats))))
				} else {
					pc := x.GetPMapLogged(n.Len * g)
					rtconds := make([]float64, x.GetNumStates())
					x2 := make([]float64, x.GetNumStates())
					for m := 0; m < numstates; m++ {
						for k := 0; k < numstates; k++ {
							x2[k] = pc.At(m, k) + n.Data[site][k]
						}
						rtconds[m] = floats.LogSumExp(x2)
					}
					tsl[p] = (floats.LogSumExp(rtconds) + (math.Log(1) - math.Log(float64(x.GammaNCats))))
				}
			}
		}
	}
	return floats.LogSumExp(tsl)
}

// calcLogLikeOneSiteSubClade calc for just a clade, starting at a node
func calcLikeOneSiteSubClade(t *Tree, inn *Node, excl bool, x *DiscreteModel, site int) float64 {
	numstates := x.NumStates
	sl := 0.0
	arr := []*Node{}
	if excl == true {
		arr = t.Rt.PostorderArrayExcl(inn)
	} else {
		arr = inn.PostorderArray()
	}
	for _, n := range arr {
		if len(n.Chs) > 0 {
			CalcLikeNode(n, x, site)
		}
		var tn *Node
		if excl == true { // calc at rt
			tn = t.Rt
		} else {
			tn = inn
		}
		if tn == n {
			if tn == t.Rt { //only happens at the root
				for i := 0; i < numstates; i++ {
					n.Data[site][i] *= x.BF[i]
				}
				sl = floats.Sum(n.Data[site])
			} else {
				//needs to get the branch length incorporated
				p := x.GetPCalc(n.Len)
				rtconds := make([]float64, x.GetNumStates())
				for j := 0; j < numstates; j++ {
					templike := 0.0
					for k := 0; k < numstates; k++ {
						templike += p.At(j, k) * n.Data[site][k]
					}
					rtconds[j] = templike
				}
				sl = floats.Sum(rtconds)
			}
		}
	}
	return sl
}

func calcLikeOneSiteSubCladeGamma(t *Tree, inn *Node, excl bool, x *DiscreteModel, site int) float64 {
	numstates := x.NumStates
	sl := 0.0
	arr := []*Node{}
	if excl == true {
		arr = t.Rt.PostorderArrayExcl(inn)
	} else {
		arr = inn.PostorderArray()
	}
	for _, g := range x.GammaCats {
		for _, n := range arr {
			if len(n.Chs) > 0 {
				CalcLikeNodeGamma(n, x, site, g)
			}
			var tn *Node
			if excl == true { // calc at rt
				tn = t.Rt
			} else {
				tn = inn
			}
			if tn == n {
				if tn == t.Rt { //only happens at the root
					for i := 0; i < numstates; i++ {
						n.Data[site][i] *= x.BF[i]
					}
					//sl = floats.Sum(n.Data[site])
					sl += floats.Sum(n.Data[site]) * (1. / float64(x.GammaNCats))
				} else {
					//needs to get the branch length incorporated
					p := x.GetPCalc(n.Len * g)
					rtconds := make([]float64, x.GetNumStates())
					for j := 0; j < numstates; j++ {
						templike := 0.0
						for k := 0; k < numstates; k++ {
							templike += p.At(j, k) * n.Data[site][k]
						}
						rtconds[j] = templike
					}
					sl += floats.Sum(rtconds) * (1. / float64(x.GammaNCats))
				}
			}
		}
	}
	return sl
}

// PCalcLogLikePatternsSubClade parallel log likeliohood calculation including patterns
func PCalcLogLikePatternsSubClade(t *Tree, n *Node, excl bool, x *DiscreteModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	x.EmptyPLDict()
	var lkfun1 func(t *Tree, inn *Node, excl bool, x *DiscreteModel, site int) float64
	var lkfun2 func(t *Tree, inn *Node, excl bool, x *DiscreteModel, jobs <-chan int, results chan<- LikeResult)
	if x.GammaNCats != 0 { // gamma
		lkfun1 = calcLogLikeOneSiteSubCladeGamma
		lkfun2 = calcLogLikeSubCladeWorkGamma
	} else {
		lkfun1 = calcLogLikeOneSiteSubClade
		lkfun2 = calcLogLikeSubCladeWork
	}
	fl += lkfun1(t, n, excl, x, 0) * patternval[0]
	for i := 0; i < wks; i++ {
		go lkfun2(t, n, excl, x, jobs, results)
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

// PCalcLikePatternsSubClade parallel log likeliohood calculation including patterns
func PCalcLikePatternsSubClade(t *Tree, n *Node, excl bool, x *DiscreteModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	var lkfun1 func(t *Tree, inn *Node, excl bool, x *DiscreteModel, site int) float64
	var lkfun2 func(t *Tree, inn *Node, excl bool, x *DiscreteModel, jobs <-chan int, results chan<- LikeResult)
	if x.GammaNCats != 0 { // gamma
		lkfun1 = calcLikeOneSiteSubCladeGamma
		lkfun2 = calcLikeSubCladeWorkGamma
	} else {
		lkfun1 = calcLikeOneSiteSubClade
		lkfun2 = calcLikeSubCladeWork
	}
	fl += math.Log(lkfun1(t, n, excl, x, 0)) * patternval[0]
	for i := 0; i < wks; i++ {
		go lkfun2(t, n, excl, x, jobs, results)
	}
	for i := 1; i < nsites; i++ {
		jobs <- i
	}
	close(jobs)
	rr := LikeResult{}
	for i := 1; i < nsites; i++ {
		rr = <-results
		fl += (math.Log(rr.value) * patternval[rr.site])
	}
	return
}

// calcLikeSubCladeWork this is intended for a worker that will be executing this per site
func calcLikeSubCladeWork(t *Tree, inn *Node, excl bool, x *DiscreteModel, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	numstates := x.NumStates
	arr := []*Node{}
	if excl == true {
		arr = t.Rt.PostorderArrayExcl(inn)
	} else {
		arr = inn.PostorderArray()
	}
	for j := range jobs {
		sl := 0.0
		for _, n := range arr {
			if len(n.Chs) > 0 {
				CalcLikeNode(n, x, j)
			}
			var tn *Node
			if excl == true { // calc at rt
				tn = t.Rt
			} else {
				tn = inn
			}
			if tn == n {
				if tn == t.Rt { //only happens at the root
					for i := 0; i < numstates; i++ {
						n.Data[j][i] *= x.BF[i]
					}
					sl = floats.Sum(n.Data[j])
				} else {
					p := x.GetPCalc(n.Len)
					rtconds := make([]float64, x.GetNumStates())
					for m := 0; m < numstates; m++ {
						templike := 0.0
						for k := 0; k < numstates; k++ {
							templike += p.At(m, k) * n.Data[j][k]
						}
						rtconds[m] = templike
					}
					sl = floats.Sum(rtconds)
				}
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

// calcLikeSubCladeWork this is intended for a worker that will be executing this per site
func calcLikeSubCladeWorkGamma(t *Tree, inn *Node, excl bool, x *DiscreteModel, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	numstates := x.NumStates
	arr := []*Node{}
	if excl == true {
		arr = t.Rt.PostorderArrayExcl(inn)
	} else {
		arr = inn.PostorderArray()
	}
	for j := range jobs {
		sl := 0.0
		for _, g := range x.GammaCats {
			for _, n := range arr {
				if len(n.Chs) > 0 {
					CalcLikeNodeGamma(n, x, j, g)
				}
				var tn *Node
				if excl == true { // calc at rt
					tn = t.Rt
				} else {
					tn = inn
				}
				if tn == n {
					if tn == t.Rt { //only happens at the root
						for i := 0; i < numstates; i++ {
							n.Data[j][i] *= x.BF[i]
						}
						//sl = floats.Sum(n.Data[j])
						sl += floats.Sum(n.Data[j]) * (1. / float64(x.GammaNCats))
					} else {
						p := x.GetPCalc(n.Len * g)
						rtconds := make([]float64, x.GetNumStates())
						for m := 0; m < numstates; m++ {
							templike := 0.0
							for k := 0; k < numstates; k++ {
								templike += p.At(m, k) * n.Data[j][k]
							}
							rtconds[m] = templike
						}
						sl += floats.Sum(rtconds) * (1. / float64(x.GammaNCats))
					}
				}
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

// CalcLogLikeSubCladeWork this is intended for a worker that will be executing this per site
func calcLogLikeSubCladeWork(t *Tree, inn *Node, excl bool, x *DiscreteModel, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	numstates := x.NumStates
	arr := []*Node{}
	if excl == true {
		arr = t.Rt.PostorderArrayExcl(inn)
	} else {
		arr = inn.PostorderArray()
	}
	for j := range jobs {
		sl := 0.0
		for _, n := range arr {
			if len(n.Chs) > 0 {
				CalcLogLikeNode(n, x, j)
			}
			var tn *Node
			if excl == true { // calc at rt
				tn = t.Rt
			} else {
				tn = inn
			}
			if tn == n {
				if tn == t.Rt { //only happens at the root
					for i := 0; i < numstates; i++ {
						n.Data[j][i] += math.Log(x.BF[i])
					}
					sl = floats.LogSumExp(n.Data[j])
				} else {
					p := x.GetPMapLogged(n.Len)
					rtconds := make([]float64, x.GetNumStates())
					x2 := make([]float64, x.GetNumStates())
					for m := 0; m < numstates; m++ {
						for k := 0; k < numstates; k++ {
							x2[k] = p.At(m, k) + n.Data[j][k]
						}
						rtconds[m] = floats.LogSumExp(x2)
					}
					sl = floats.LogSumExp(rtconds)
				}
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

func calcLogLikeSubCladeWorkGamma(t *Tree, inn *Node, excl bool, x *DiscreteModel, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	numstates := x.NumStates
	arr := []*Node{}
	if excl == true {
		arr = t.Rt.PostorderArrayExcl(inn)
	} else {
		arr = inn.PostorderArray()
	}
	for j := range jobs {
		tsl := make([]float64, x.GammaNCats)
		for p, g := range x.GammaCats {

			for _, n := range arr {
				if len(n.Chs) > 0 {
					CalcLogLikeNodeGamma(n, x, j, g)
				}
				var tn *Node
				if excl == true { // calc at rt
					tn = t.Rt
				} else {
					tn = inn
				}
				if tn == n {
					if tn == t.Rt { //only happens at the root
						for i := 0; i < numstates; i++ {
							n.Data[j][i] += math.Log(x.BF[i])
						}
						tsl[p] = (floats.LogSumExp(n.Data[j]) + (math.Log(1) - math.Log(float64(x.GammaNCats))))
					} else {
						pc := x.GetPMapLogged(n.Len * g)
						rtconds := make([]float64, x.GetNumStates())
						x2 := make([]float64, x.GetNumStates())
						for m := 0; m < numstates; m++ {
							for k := 0; k < numstates; k++ {
								x2[k] = pc.At(m, k) + n.Data[j][k]
							}
							rtconds[m] += floats.LogSumExp(x2)
						}
						tsl[p] = (floats.LogSumExp(rtconds) + (math.Log(1) - math.Log(float64(x.GammaNCats))))
					}
				}
			}
		}
		results <- LikeResult{value: floats.LogSumExp(tsl), site: j}
	}
}
