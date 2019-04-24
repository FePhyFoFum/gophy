package gophy

import (
	"math"

	"gonum.org/v1/gonum/floats"
)

/*
 This is for calculating likelihoods for multistate with multiple models
 (each edge has a model assigned)
*/

/*=======================================================================
	These are the parallel calls for log like and likelihoods
  =======================================================================
*/

// PCalcLikePatternsMSMUL this will calculate like in parallel
//     It will make a thread pool then calculate the LogLike for the first site
//     with CalcLikeOneSiteMS and the rest in the pool with CalcLikeWorkMS
//     This is with patterns.
func PCalcLikePatternsMSMUL(t *Tree, models []StateModel, nodemodels map[*Node]int, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	for _, x := range models {
		x.EmptyPDict()
	}
	fl += math.Log(CalcLikeOneSiteMSMUL(t, models, nodemodels, 0)) * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcLikeWorkMSMUL(t, models, nodemodels, jobs, results)
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

// PCalcLogLikePatternsMSMUL this will calculate loglike in parallel
//     It will make a thread pool then calculate the LogLike for the first site
//     with CalcLogLikeOneSiteMS and the rest in the pool with CalcLogLikeWorkMS
//     This is with patterns. Probably use this one
func PCalcLogLikePatternsMSMUL(t *Tree, models []StateModel, nodemodels map[*Node]int, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	for _, x := range models {
		x.EmptyPDict()
	}
	fl += CalcLogLikeOneSiteMSMUL(t, models, nodemodels, 0) * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcLogLikeWorkMSMUL(t, models, nodemodels, jobs, results)
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

// PCalcLogLikeMSMUL this will calculate log like in parallel
//     It will make a thread pool then calculate the LogLike for the first site
//     with CalcLogLikeOneSiteMS and the rest in the pool with CalcLogLikeWorkMS
//     This is without patterns.
func PCalcLogLikeMSMUL(t *Tree, models []StateModel, nodemodels map[*Node]int, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	for _, x := range models {
		x.EmptyPDict()
	}
	fl += CalcLogLikeOneSiteMSMUL(t, models, nodemodels, 0)
	for i := 0; i < wks; i++ {
		go CalcLogLikeWorkMSMUL(t, models, nodemodels, jobs, results)
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

// PCalcLikeMSMUL this will calculate log like in parallel
//     It will make a thread pool then calculate the LogLike for the first site
//     with CalcLikeOneSiteMS and the rest in the pool with CalcLikeWorkMS
//     This is without patterns.
func PCalcLikeMSMUL(t *Tree, models []StateModel, nodemodels map[*Node]int, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	for _, x := range models {
		x.EmptyPDict()
	}
	fl += math.Log(CalcLikeOneSiteMSMUL(t, models, nodemodels, 0))
	for i := 0; i < wks; i++ {
		go CalcLikeWorkMSMUL(t, models, nodemodels, jobs, results)
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

/*=======================================================================
	These are the non parallel calls to calculate likelihoods for one site
  =======================================================================
*/

// CalcLogLikeOneSiteMSMUL is used to calculate the loglike of one site
// probably used to populate the PDict in the DNA Model so that we can reuse the calculations
func CalcLogLikeOneSiteMSMUL(t *Tree, models []StateModel, nodemodels map[*Node]int, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			CalcLogLikeNodeMSMUL(n, models, nodemodels, site)
		}
		if t.Rt == n {
			for i := 0; i < 4; i++ {
				t.Rt.Data[site][i] += math.Log(models[nodemodels[n]].GetBF()[i])
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLikeOneSiteMSMUL is used to calculate the like of one site
// probably used to populate the PDict in the DNA Model so that we can reuse the calculations
func CalcLikeOneSiteMSMUL(t *Tree, models []StateModel, nodemodels map[*Node]int, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		x := models[nodemodels[n]]
		if len(n.Chs) > 0 {
			CalcLikeNodeMSMUL(n, models, nodemodels, site)
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

/*=======================================================================
	These are the worker (parallel) calls to calculate likelihoods for one site
  =======================================================================
*/

// CalcLogLikeWorkMSMUL calculates log likelihood for one site using CalcLogLikeNodeMS
//    Work refers to this being part of the worker pool
func CalcLogLikeWorkMSMUL(t *Tree, models []StateModel, nodemodels map[*Node]int, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			x := models[nodemodels[n]]
			if len(n.Chs) > 0 {
				CalcLogLikeNodeMSMUL(n, models, nodemodels, j)
			}
			if t.Rt == n {
				for i := 0; i < 4; i++ {
					t.Rt.Data[j][i] += math.Log(x.GetBF()[i])
				}
				sl = floats.LogSumExp(t.Rt.Data[j])
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

// CalcLikeWorkMSMUL calculates likelihood for one site using CalcLikeNodeMS
//    Work refers to this being part of the worker pool
func CalcLikeWorkMSMUL(t *Tree, models []StateModel, nodemodels map[*Node]int, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			x := models[nodemodels[n]]
			if len(n.Chs) > 0 {
				CalcLikeNodeMSMUL(n, models, nodemodels, j)
			}
			if t.Rt == n {
				for i := 0; i < x.GetNumStates(); i++ {
					t.Rt.Data[j][i] *= models[nodemodels[n]].GetBF()[i]
				}
				sl = floats.Sum(t.Rt.Data[j])
			}
		}
		results <- LikeResult{value: sl, site: j}
	}
}

/*=======================================================================
	These are the parts that actually calculate likelihood for one node.
	These are all called from the functions above
  =======================================================================
*/

// CalcLogLikeNodeMSMUL calculates log likelihood for node for multistate
func CalcLogLikeNodeMSMUL(nd *Node, models []StateModel, nodemodels map[*Node]int, site int) {
	for i := 0; i < 4; i++ {
		nd.Data[site][i] = 0.
	}
	x1 := 0.0
	x2 := []float64{0.0, 0.0, 0.0, 0.0}
	for _, c := range nd.Chs {
		P := models[nodemodels[c]].GetPMap(c.Len)
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

// CalcLikeNodeMSMUL calculate the likelihood of a node for multistate
func CalcLikeNodeMSMUL(nd *Node, models []StateModel, nodemodels map[*Node]int, site int) {
	for i := 0; i < models[nodemodels[nd]].GetNumStates(); i++ {
		nd.Data[site][i] = 1.
	}
	x1 := 0.0
	x2 := 0.0
	for _, c := range nd.Chs {
		x := models[nodemodels[c]]
		P := x.GetPMap(c.Len)
		if len(c.Chs) == 0 {
			for i := 0; i < x.GetNumStates(); i++ {
				x1 = 0.0
				for j := 0; j < x.GetNumStates(); j++ {
					x1 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] *= x1
			}
		} else {
			for i := 0; i < x.GetNumStates(); i++ {
				x2 = 0.0
				for j := 0; j < x.GetNumStates(); j++ {
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

// TPconditionalsMSMUL calculate TpConds vector for multistate
func TPconditionalsMSMUL(models []StateModel, nodemodels map[*Node]int, node *Node, patternval []float64) {
	if len(node.Chs) > 0 {
		for s := range patternval {
			for j := 0; j < models[nodemodels[node]].GetNumStates(); j++ {
				node.TpConds[s][j] = 1.
				for _, i := range node.Chs {
					node.TpConds[s][j] *= i.RtConds[s][j]
				}
			}
		}
	}
}

// RTconditionalsMSMUL calculate RtConds vector for multistate
func RTconditionalsMSMUL(models []StateModel, nodemodels map[*Node]int, node *Node, patternval []float64) {
	x := models[nodemodels[node]]
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

// RVconditionalsMSMUL calculate RvTpConds vector for the Parent of node and
//     rely on the Par RvConds for multistate
func RVconditionalsMSMUL(models []StateModel, nodemodels map[*Node]int, node *Node, patternval []float64) {
	x := models[nodemodels[node]]
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

// RVTPconditionalsMSMUL calculate RvConds vector multistate
func RVTPconditionalsMSMUL(models []StateModel, nodemodels map[*Node]int, node *Node, patternval []float64) {
	x := models[nodemodels[node]]
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

// CalcLikeFrontBackMSMUL ...
func CalcLikeFrontBackMSMUL(models []StateModel, nodemodels map[*Node]int, tree *Tree, patternval []float64) {
	for _, n := range tree.Post {
		x := models[nodemodels[n]]
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
		TPconditionalsMSMUL(models, nodemodels, c, patternval)
		//take the tip cond to the rt
		RTconditionalsMSMUL(models, nodemodels, c, patternval) // calculate from tpcond to rtcond
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
			RVconditionalsMSMUL(models, nodemodels, c, patternval)
			RVTPconditionalsMSMUL(models, nodemodels, c, patternval)
		}
	}
}

// CalcAncStatesMSMUL for each node based on the calculations above
func CalcAncStatesMSMUL(models []StateModel, nodemodels map[*Node]int, tree *Tree, patternval []float64) (retstates map[*Node][][]float64) {
	CalcLikeFrontBackMSMUL(models, nodemodels, tree, patternval)
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
			x := models[nodemodels[c]]
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

// CalcStochMapMSMUL for each node based on the calculations above
func CalcStochMapMSMUL(models []StateModel, nodemodels map[*Node]int, tree *Tree, patternval []float64, time bool, from int, to int) (retstates map[*Node][][]float64) {
	CalcLikeFrontBackMSMUL(models, nodemodels, tree, patternval)
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
			x := models[nodemodels[c]]
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
