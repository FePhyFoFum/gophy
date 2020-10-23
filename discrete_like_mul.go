package gophy

import (
	"math"

	"gonum.org/v1/gonum/floats"
)

/*
 This is for calculating likelihoods for nucleotides with multiple models
*/

//PCalcLikePatternsMul parallel caclulation of likelihood with patterns
func PCalcLikePatternsMul(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int,
	patternval []float64, wks int) (fl float64) {
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
	fl += math.Log(calcLikeOneSiteMul(t, models, nodemodels, 0)) * patternval[0]
	for i := 0; i < wks; i++ {
		go calcLikeWorkMul(t, models, nodemodels, jobs, results)
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

//PCalcLikePatternsMulSubClade parallel log likeliohood calculation including patterns
func PCalcLikePatternsMulSubClade(t *Tree, n *Node, excl bool, models []*DiscreteModel,
	nodemodels map[*Node]int, patternval []float64, wks int) (fl float64) {
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
	fl += calcLikeOneSiteMulSubClade(t, n, excl, models, nodemodels, 0) * patternval[0]
	for i := 0; i < wks; i++ {
		go calcLikeMulSubCladeWork(t, n, excl, models, nodemodels, jobs, results)
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
func calcLikeMulSubCladeWork(t *Tree, inn *Node, excl bool, models []*DiscreteModel,
	nodemodels map[*Node]int, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	arr := []*Node{}
	if excl == true {
		arr = t.Rt.PostorderArrayExcl(inn)
	} else {
		arr = inn.PostorderArray()
	}
	for j := range jobs {
		sl := 0.0
		for _, n := range arr {
			x := models[nodemodels[n]]
			numstates := x.NumStates
			if len(n.Chs) > 0 {
				calcLikeNodeMul(n, x, j)
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

// calcLogLikeOneSiteSubClade calc for just a clade, starting at a node
func calcLikeOneSiteMulSubClade(t *Tree, inn *Node, excl bool, models []*DiscreteModel,
	nodemodels map[*Node]int, site int) float64 {
	sl := 0.0
	arr := []*Node{}
	if excl == true {
		arr = t.Rt.PostorderArrayExcl(inn)
	} else {
		arr = inn.PostorderArray()
	}
	for _, n := range arr {
		x := models[nodemodels[n]]
		numstates := x.NumStates
		if len(n.Chs) > 0 {
			calcLikeNodeMul(n, x, site)
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

//PCalcLogLikePatternsMul ...
func PCalcLogLikePatternsMul(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	for _, x := range models {
		x.EmptyPDict()
		x.EmptyPLDict()
	}
	fl += CalcLogLikeOneSiteMul(t, models, nodemodels, 0) * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcLogLikeWorkMul(t, models, nodemodels, jobs, results)
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

//CalcLikeOneSiteMul just one site
func calcLikeOneSiteMul(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		x := models[nodemodels[n]]
		numstates := x.NumStates
		if len(n.Chs) > 0 {
			calcLikeNodeMul(n, x, site)
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

//CalcLogLikeOneSiteMul just one site
func CalcLogLikeOneSiteMul(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		x := models[nodemodels[n]]
		numstates := x.NumStates
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

//calcLikeWorkMul this is the worker
func calcLikeWorkMul(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			x := models[nodemodels[n]]
			numstates := x.NumStates
			if len(n.Chs) > 0 {
				calcLikeNodeMul(n, x, j)
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

//CalcLogLikeWorkMul this is the worker
func CalcLogLikeWorkMul(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			x := models[nodemodels[n]]
			numstates := x.NumStates
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

//calcLikeNodeMul calculate the likelihood of a node
func calcLikeNodeMul(nd *Node, model *DiscreteModel, site int) {
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

//PCalcLikePatternsMarkedMul parallel likelihood caclulation with patterns and just update the values
func PCalcLikePatternsMarkedMul(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int, patternval []float64, wks int) (fl float64) {
	fl = 0.0
	nsites := len(patternval)
	jobs := make(chan int, nsites)
	//results := make(chan float64, nsites)
	results := make(chan LikeResult, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	//x.EmptyPDict()
	fl += math.Log(CalcLikeOneSiteMarkedMul(t, models, nodemodels, 0)) * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcLikeWorkMarkedMul(t, models, nodemodels, jobs, results)
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

// CalcLikeOneSiteMarkedMul this uses the marked machinery to recalculate
func CalcLikeOneSiteMarkedMul(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		x := models[nodemodels[n]]
		numstates := x.NumStates
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

// CalcLikeWorkMarkedMul this is intended to calculate only on the marked nodes back to teh root
func CalcLikeWorkMarkedMul(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int, jobs <-chan int, results chan<- LikeResult) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			x := models[nodemodels[n]]
			numstates := x.NumStates
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

/*
 * calculate the conditionals for ancestral calc or branch lengths
 */

//RTMultconditionals ...
func RTMultconditionals(models []*DiscreteModel, nodemodels map[*Node]int, node *Node, patternval []float64) {
	p := models[nodemodels[node]].GetPCalc(node.Len)
	numstates := models[nodemodels[node]].NumStates
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

//RVMultconditionals ...
func RVMultconditionals(models []*DiscreteModel, nodemodels map[*Node]int, node *Node, patternval []float64) {
	p := models[nodemodels[node]].GetPCalc(node.Par.Len)
	numstates := models[nodemodels[node]].NumStates
	for s := range patternval {
		for j := 0; j < numstates; j++ {
			node.Par.RvTpConds[s][j] = 0.0
			for k := 0; k < numstates; k++ {
				node.Par.RvTpConds[s][j] += p.At(j, k) * node.Par.RvConds[s][k]
			}
		}
	}
}

//CalcLikeFrontBackMult ...
func CalcLikeFrontBackMult(models []*DiscreteModel, nodemodels map[*Node]int, tree *Tree, patternval []float64) {
	for _, n := range tree.Post {
		if len(n.Chs) != 0 {
			n.TpConds = make([][]float64, len(patternval))
		}
		n.RvTpConds = make([][]float64, len(patternval))
		n.RvConds = make([][]float64, len(patternval))
		n.RtConds = make([][]float64, len(patternval))
		for i := 0; i < len(patternval); i++ {
			if len(n.Chs) != 0 {
				n.TpConds[i] = []float64{1.0, 1.0, 1.0, 1.0}
			}
			n.RvTpConds[i] = []float64{1.0, 1.0, 1.0, 1.0}
			n.RvConds[i] = []float64{1.0, 1.0, 1.0, 1.0}
			n.RtConds[i] = []float64{1.0, 1.0, 1.0, 1.0}
		}
	}
	//loglike := 0.
	for _, c := range tree.Post {
		//calculate the tip conditionals
		TPconditionals(models[nodemodels[c]], c, patternval)
		//take the tip cond to the rt
		RTMultconditionals(models, nodemodels, c, patternval) // calculate from tpcond to rtcond
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
			RVMultconditionals(models, nodemodels, c, patternval)
			RVTPconditionals(models[nodemodels[c]], c, patternval)
		}
	}
}
