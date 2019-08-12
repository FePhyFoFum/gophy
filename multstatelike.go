package gophy

import (
	"math"

	"gonum.org/v1/gonum/floats"
)

/*
 This is for calculating likelihoods for multistate
*/

/*=======================================================================
	These are the parallel calls for log like and likelihoods
  =======================================================================
*/

// PCalcLikePatternsMS this will calculate like in parallel
//     It will make a thread pool then calculate the LogLike for the first site
//     with CalcLikeOneSiteMS and the rest in the pool with CalcLikeWorkMS
//     This is with patterns.
func PCalcLikePatternsMS(t *Tree, x StateModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += math.Log(CalcLikeOneSiteMS(t, x, 0)) * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcLikeWorkMS(t, x, jobs, results)
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

// PCalcLogLikePatternsMS this will calculate loglike in parallel
//     It will make a thread pool then calculate the LogLike for the first site
//     with CalcLogLikeOneSiteMS and the rest in the pool with CalcLogLikeWorkMS
//     This is with patterns. Probably use this one
func PCalcLogLikePatternsMS(t *Tree, x StateModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += CalcLogLikeOneSiteMS(t, x, 0) * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcLogLikeWorkMS(t, x, jobs, results)
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
	//fmt.Println(fl)
	//TODO: check this and add as an option to the model
	//fl = fl - math.Log(CalcAscBiasVal(t, x))
	//fmt.Println(fl)
	return
}

// PCalcLogLikeMS this will calculate log like in parallel
//     It will make a thread pool then calculate the LogLike for the first site
//     with CalcLogLikeOneSiteMS and the rest in the pool with CalcLogLikeWorkMS
//     This is without patterns.
func PCalcLogLikeMS(t *Tree, x StateModel, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += CalcLogLikeOneSiteMS(t, x, 0)
	for i := 0; i < wks; i++ {
		go CalcLogLikeWorkMS(t, x, jobs, results)
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

// PCalcLikeMS this will calculate log like in parallel
//     It will make a thread pool then calculate the LogLike for the first site
//     with CalcLikeOneSiteMS and the rest in the pool with CalcLikeWorkMS
//     This is without patterns.
func PCalcLikeMS(t *Tree, x StateModel, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += math.Log(CalcLikeOneSiteMS(t, x, 0))
	for i := 0; i < wks; i++ {
		go CalcLikeWorkMS(t, x, jobs, results)
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

// PCalcLikePatternsMarkedMS this will calculate like in parallel
//     It will make a thread pool then calculate the LogLike for the first site
//     with CalcLikeOneSiteMarkedMS and the rest in the pool with CalcLikeWorkMarkedMS
//     This is with patterns. It will assume that nodes are going to be marked in between
//     calls to this function.
func PCalcLikePatternsMarkedMS(t *Tree, x StateModel, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	//x.EmptyPDict()
	fl += math.Log(CalcLikeOneSiteMarkedMS(t, x, 0)) * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcLikeWorkMarkedMS(t, x, jobs, results)
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

// PCalcLogLikeBackMS this will calculate loglike in parallel
//     It will make a thread pool then calculate the LogLike for the first site
//     with CalcLogLikeOneSiteBackMS and the rest in the pool with CalcLogLikeWorkBackMS
//     This is without patterns. It will start at a node n and go back to the root.
func PCalcLogLikeBackMS(t *Tree, n *Node, x StateModel, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += CalcLogLikeOneSiteBackMS(t, n, x, 0)
	for i := 0; i < wks; i++ {
		go CalcLogLikeWorkBackMS(t, n, x, jobs, results)
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

// PCalcLogLikeMarkedMS this will calculate loglike in parallel
//     It will make a thread pool then calculate the LogLike for the first site
//     with CalcLogLikeOneSiteMarkedMS and the rest in the pool with CalcLogLikeWorkMarkedMS
//     This is without patterns. It will assume that nodes are going to be marked in between
//     calls to this function.
func PCalcLogLikeMarkedMS(t *Tree, x StateModel, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	//x.EmptyPDict()
	fl += CalcLogLikeOneSiteMarkedMS(t, x, 0)
	for i := 0; i < wks; i++ {
		go CalcLogLikeWorkMarkedMS(t, x, jobs, results)
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

/*=======================================================================
	These are the non parallel calls to calculate likelihoods for one site
  =======================================================================
*/

// CalcLogLikeOneSiteMS is used to calculate the loglike of one site
// probably used to populate the PDict in the DNA Model so that we can reuse the calculations
func CalcLogLikeOneSiteMS(t *Tree, x StateModel, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			CalcLogLikeNodeMS(n, x, site)
		}
		if t.Rt == n {
			for i := 0; i < x.GetNumStates(); i++ {
				t.Rt.Data[site][i] += math.Log(x.GetBF()[i])
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLikeOneSiteMS is used to calculate the like of one site
// probably used to populate the PDict in the DNA Model so that we can reuse the calculations
func CalcLikeOneSiteMS(t *Tree, x StateModel, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			CalcLikeNodeMS(n, x, site)
		}
		if t.Rt == n {
			for i := 0; i < x.GetNumStates(); i++ {
				t.Rt.Data[site][i] *= x.GetBF()[i]
			}
			sl = floats.Sum(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLogLikeOneSiteBackMS is used to calculate the loglike of one site
// probably used to populate the PDict in the DNA Model so that we can reuse the calculations
// This starts at node nb and goes to the root.
func CalcLogLikeOneSiteBackMS(t *Tree, nb *Node, x StateModel, site int) float64 {
	sl := 0.0
	going := true
	cur := nb
	for going {
		if len(cur.Chs) > 0 {
			CalcLogLikeNodeMS(cur, x, site)
		}
		if cur == t.Rt {
			for i := 0; i < x.GetNumStates(); i++ {
				t.Rt.Data[site][i] += math.Log(x.GetBF()[i])
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
			going = false
			break
		}
		cur = cur.Par
	}
	return sl
}

// CalcLogLikeOneSiteMarkedMS calculates log likelihood for one site using CalcLogLikeNodeMS
//    This is outside of teh worker pool (probably doing this first to populate the P matrix)
func CalcLogLikeOneSiteMarkedMS(t *Tree, x StateModel, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			if n.Marked == true {
				CalcLogLikeNodeMS(n, x, site)
				if n != t.Rt {
					n.Par.Marked = true
				}
			}
		}
		if t.Rt == n && n.Marked == true {
			for i := 0; i < x.GetNumStates(); i++ {
				t.Rt.Data[site][i] += math.Log(x.GetBF()[i])
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
		} else {
			sl = floats.LogSumExp(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLikeOneSiteMarkedMS calculates likelihood for one site using CalcLikeNodeMS
//    This is outside of teh worker pool (probably doing this first to populate the P matrix)
func CalcLikeOneSiteMarkedMS(t *Tree, x StateModel, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			if n.Marked == true {
				CalcLikeNodeMS(n, x, site)
				if n != t.Rt {
					n.Par.Marked = true
				}
			}
		}
		if t.Rt == n && n.Marked == true {
			for i := 0; i < x.GetNumStates(); i++ {
				t.Rt.Data[site][i] *= x.GetBF()[i]
			}
			sl = floats.Sum(t.Rt.Data[site])
		} else {
			sl = floats.Sum(t.Rt.Data[site])
		}
	}
	return sl
}

/*=======================================================================
	These are the worker (parallel) calls to calculate likelihoods for one site
  =======================================================================
*/

// CalcLogLikeWorkMS calculates log likelihood for one site using CalcLogLikeNodeMS
//    Work refers to this being part of the worker pool
func CalcLogLikeWorkMS(t *Tree, x StateModel, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				CalcLogLikeNodeMS(n, x, j)
			}
			if t.Rt == n {
				for i := 0; i < x.GetNumStates(); i++ {
					t.Rt.Data[j][i] += math.Log(x.GetBF()[i])
				}
				sl = floats.LogSumExp(t.Rt.Data[j])
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

// CalcLikeWorkMS calculates likelihood for one site using CalcLikeNodeMS
//    Work refers to this being part of the worker pool
func CalcLikeWorkMS(t *Tree, x StateModel, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				CalcLikeNodeMS(n, x, j)
			}
			if t.Rt == n {
				for i := 0; i < x.GetNumStates(); i++ {
					t.Rt.Data[j][i] *= x.GetBF()[i]
				}
				sl = floats.Sum(t.Rt.Data[j])
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

// CalcLogLikeWorkBackMS only calculates the log likelihood back to the root for multistate
//    Work refers to this being part of the worker pool
func CalcLogLikeWorkBackMS(t *Tree, nb *Node, x StateModel, jobs <-chan int, results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		going := true
		cur := nb
		for going {
			if len(cur.Chs) > 0 {
				CalcLogLikeNodeMS(cur, x, j)
			}
			if cur == t.Rt {
				for i := 0; i < x.GetNumStates(); i++ {
					t.Rt.Data[j][i] += math.Log(x.GetBF()[i])
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

// CalcLogLikeWorkMarkedMS this should only calculate log like of the marked nodes back to the root
//    Work refers to this being part of the worker pool
func CalcLogLikeWorkMarkedMS(t *Tree, x StateModel, jobs <-chan int, results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				if n.Marked == true {
					CalcLogLikeNodeMS(n, x, j)
				}
			}
			if t.Rt == n && n.Marked == true {
				for i := 0; i < x.GetNumStates(); i++ {
					t.Rt.Data[j][i] += math.Log(x.GetBF()[i])
				}
				sl = floats.LogSumExp(t.Rt.Data[j])
			} else {
				sl = floats.LogSumExp(t.Rt.Data[j])
			}
		}
		results <- sl
	}
}

// CalcLikeWorkMarkedMS this should only calculate like of the marked nodes back to the root
//    Work refers to this being part of the worker pool
func CalcLikeWorkMarkedMS(t *Tree, x StateModel, jobs <-chan int, results chan<- LikeResult) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				if n.Marked == true {
					CalcLikeNodeMS(n, x, j)
				}
			}
			if t.Rt == n && n.Marked == true {
				for i := 0; i < x.GetNumStates(); i++ {
					t.Rt.Data[j][i] *= x.GetBF()[i]
				}
				sl = floats.Sum(t.Rt.Data[j])
			} else {
				sl = floats.Sum(t.Rt.Data[j])
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

/*======================================================================
This will calculate the likelihood of of each possible type of invariant
site for use in ascertainment bias correction.
========================================================================
*/

// CalcAscBiasVal creates dummy invariant characters for each possible state and sums their likelihoods to yield a Prob(invariant)
// Used in ascertainment bias correction like: L_marginal = L_unconditioned / Prob(invariant)
func CalcAscBiasVal(t *Tree, x StateModel) float64 {
	sl := 0.0
	if t.Pre[1].Num == 0 { //number the nodes if this hasn't been done
		count := 0
		for _, n := range t.Pre {
			n.Num = count
			count++
		}
	}
	stateVecs := make(map[int][]float64) // key is node num, val is [states]
	for state := 0; state < x.GetNumStates(); state++ {
		for _, n := range t.Post {
			if len(n.Chs) == 0 {
				if _, ok := stateVecs[n.Num]; !ok {
					stateVecs[n.Num] = make([]float64, x.GetNumStates()) // zero valued initialized
				} else {
					for i := 0; i < x.GetNumStates(); i++ {
						stateVecs[n.Num][i] = 0.0
					}
				}
				stateVecs[n.Num][state] = 1.0
			} else if len(n.Chs) > 0 {
				if _, ok := stateVecs[n.Num]; !ok {
					stateVecs[n.Num] = make([]float64, x.GetNumStates())
					for i := 0; i < x.GetNumStates(); i++ {
						stateVecs[n.Num][i] = 1.0
					}
				} else {
					for i := 0; i < x.GetNumStates(); i++ {
						stateVecs[n.Num][i] = 1.0
					}
				}
				x1 := 0.0
				x2 := 0.0
				for _, c := range n.Chs {
					P := x.GetPMap(c.Len)
					if len(c.Chs) == 0 {
						for i := 0; i < x.GetNumStates(); i++ {
							x1 = 0.0
							for j := 0; j < x.GetNumStates(); j++ {
								x1 += P.At(i, j) * stateVecs[c.Num][j]
							}
							stateVecs[n.Num][i] *= x1
						}
					} else {
						for i := 0; i < x.GetNumStates(); i++ {
							x2 = 0.0
							for j := 0; j < x.GetNumStates(); j++ {
								x2 += P.At(i, j) * stateVecs[c.Num][j] //c.Data[site][j]
							}
							stateVecs[n.Num][i] *= x2 //floats.LogSumExp(x2)

						}
					}
				}
			}
			if t.Rt == n {
				for i := 0; i < x.GetNumStates(); i++ {
					stateVecs[t.Rt.Num][i] *= x.GetBF()[i]
				}
				sl += floats.Sum(stateVecs[t.Rt.Num])
			}

		}
	}
	return 1 - sl //TODO: check this
}

/*=======================================================================
	These are the parts that actually calculate likelihood for one node.
	These are all called from the functions above
  =======================================================================
*/

// CalcLogLikeNodeMS calculates log likelihood for node for multistate
func CalcLogLikeNodeMS(nd *Node, model StateModel, site int) {
	for i := 0; i < model.GetNumStates(); i++ {
		nd.Data[site][i] = 0.
	}
	x1 := 0.0
	//x2 := []float64{0.0, 0.0, 0.0, 0.0}
	var x2 []float64
	for i := 0; i < model.GetNumStates(); i++ {
		x2 = append(x2, 0.0)
	}
	for _, c := range nd.Chs {
		P := model.GetPMap(c.Len)
		if len(c.Chs) == 0 {
			for i := 0; i < model.GetNumStates(); i++ {
				x1 = 0.0
				for j := 0; j < model.GetNumStates(); j++ {
					x1 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] += math.Log(x1)
			}
		} else {
			for i := 0; i < model.GetNumStates(); i++ {
				for j := 0; j < model.GetNumStates(); j++ {
					x2[j] = math.Log(P.At(i, j)) + c.Data[site][j]
				}
				nd.Data[site][i] += floats.LogSumExp(x2)
			}
		}
	}
}

// CalcLikeNodeMS calculate the likelihood of a node for multistate
func CalcLikeNodeMS(nd *Node, model StateModel, site int) {
	for i := 0; i < model.GetNumStates(); i++ {
		nd.Data[site][i] = 1.
	}
	x1 := 0.0
	x2 := 0.0
	for _, c := range nd.Chs {
		P := model.GetPMap(c.Len)
		if len(c.Chs) == 0 {
			for i := 0; i < model.GetNumStates(); i++ {
				x1 = 0.0
				for j := 0; j < model.GetNumStates(); j++ {
					x1 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] *= x1
			}
		} else {
			for i := 0; i < model.GetNumStates(); i++ {
				x2 = 0.0
				for j := 0; j < model.GetNumStates(); j++ {
					x2 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] *= x2 //floats.LogSumExp(x2)
			}
		}
	}
}

/*=======================================================================
	These are for calculating ancestral states, stochastic mapping,
	and branch length optimization. To do this we calculating additional
	vectors for each edge in order to speed up the calculation (or
	to calclate the derivatives)
  =======================================================================

 toward tip
tpcond  X
        | | rvtpcond
        | ^
        | |
        v |
        | | rvcond
rtcond  x
 toward root

 calculate
 postorder
   TPconditionalsMS
   RTconditionalsMS
 preorder
   RVconditionalsMS
   RVTPconditionalsMS
*/

// TPconditionalsMS calculate TpConds vector for multistate
func TPconditionalsMS(x StateModel, node *Node, patternval []float64) {
	if len(node.Chs) > 0 {
		for s := range patternval {
			for j := 0; j < x.GetNumStates(); j++ {
				node.TpConds[s][j] = 1.
				for _, i := range node.Chs {
					node.TpConds[s][j] *= i.RtConds[s][j]
				}
			}
		}
	}
}

// RTconditionalsMS calculate RtConds vector for multistate
func RTconditionalsMS(x StateModel, node *Node, patternval []float64) {
	p := x.GetPCalc(node.Len)
	for s := range patternval {
		for j := 0; j < x.GetNumStates(); j++ {
			templike := 0.0
			for k := 0; k < x.GetNumStates(); k++ {
				templike += p.At(j, k) * node.TpConds[s][k]
			}
			node.RtConds[s][j] = templike
		}
	}
}

// RVconditionalsMS calculate RvTpConds vector for the Parent of node and
//     rely on the Par RvConds for multistate
func RVconditionalsMS(x StateModel, node *Node, patternval []float64) {
	p := x.GetPCalc(node.Par.Len)
	for s := range patternval {
		for j := 0; j < x.GetNumStates(); j++ {
			node.Par.RvTpConds[s][j] = 0.0
			for k := 0; k < x.GetNumStates(); k++ {
				node.Par.RvTpConds[s][j] += p.At(j, k) * node.Par.RvConds[s][k]
			}
		}
	}
}

// RVTPconditionalsMS calculate RvConds vector multistate
func RVTPconditionalsMS(x StateModel, node *Node, patternval []float64) {
	for s := range patternval {
		for j := 0; j < x.GetNumStates(); j++ {
			node.RvConds[s][j] = node.Par.RvTpConds[s][j]
		}
		for _, oc := range node.Par.Chs {
			if node == oc {
				continue
			}
			for j := 0; j < x.GetNumStates(); j++ {
				node.RvConds[s][j] *= oc.RtConds[s][j]
			}
		}
	}
}

// CalcLikeFrontBackMS ...
func CalcLikeFrontBackMS(x StateModel, tree *Tree, patternval []float64) {
	for _, n := range tree.Post {
		if len(n.Chs) != 0 {
			n.TpConds = make([][]float64, len(patternval))
		}
		n.RvTpConds = make([][]float64, len(patternval))
		n.RvConds = make([][]float64, len(patternval))
		n.RtConds = make([][]float64, len(patternval))
		for i := 0; i < len(patternval); i++ {
			if len(n.Chs) != 0 {
				n.TpConds[i] = make([]float64, x.GetNumStates())
				for j := 0; j < x.GetNumStates(); j++ {
					n.TpConds[i][j] = 1.0
				}
			}
			n.RvTpConds[i] = make([]float64, x.GetNumStates())
			n.RvConds[i] = make([]float64, x.GetNumStates())
			n.RtConds[i] = make([]float64, x.GetNumStates())
			for j := 0; j < x.GetNumStates(); j++ {
				n.RvTpConds[i][j] = 1.0
				n.RvConds[i][j] = 1.0
				n.RtConds[i][j] = 1.0
			}
		}
	}
	//loglike := 0.
	for _, c := range tree.Post {
		//calculate the tip conditionals
		TPconditionalsMS(x, c, patternval)
		//take the tip cond to the rt
		RTconditionalsMS(x, c, patternval) // calculate from tpcond to rtcond
		/*if c == tree.Rt { // turn on if you want likelihoods
			for s := range patternval {
				tempretlike := 0.
				for i := 0; i < 4; i++ {
					tempretlike += (c.TpConds[s][i] * x.GetBF()[i])
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
			RVconditionalsMS(x, c, patternval)
			RVTPconditionalsMS(x, c, patternval)
		}
	}
}

// CalcAncStatesMS for each node based on the calculations above
func CalcAncStatesMS(x StateModel, tree *Tree, patternval []float64) (retstates map[*Node][][]float64) {
	CalcLikeFrontBackMS(x, tree, patternval)
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
		for _, c := range tree.Pre {
			if len(c.Chs) == 0 {
				continue
			}
			//fmt.Println(c.Newick(true))
			retstates[c][i] = make([]float64, x.GetNumStates())
			if c == tree.Rt {
				su := 0.
				for j, s := range c.RtConds[i] {
					su += (s * x.GetBF()[j])
				}
				for j, s := range c.RtConds[i] {
					retstates[c][i][j] = (s * x.GetBF()[j]) / su
				}
				//fmt.Print("\n")
			} else {
				p := x.GetPCalc(c.Len)
				//need subtree 1
				s1probs := c.TpConds
				//need subtree 2
				s2probs := c.RvConds
				tv := make([]float64, x.GetNumStates())
				for j := 0; j < x.GetNumStates(); j++ {
					for k := 0; k < x.GetNumStates(); k++ {
						tv[j] += (s1probs[i][j] * p.At(j, k) * s2probs[i][k])
					}
					tv[j] *= x.GetBF()[j]
				}
				su := 0.
				for _, s := range tv {
					su += s
				}
				for j, s := range tv {
					retstates[c][i][j] = s / su
				}
			}
		}
	}
	return
}

// CalcStochMapMS for each node based on the calculations above
func CalcStochMapMS(x StateModel, tree *Tree, patternval []float64, time bool, from int, to int) (retstates map[*Node][][]float64) {
	CalcLikeFrontBackMS(x, tree, patternval)
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
