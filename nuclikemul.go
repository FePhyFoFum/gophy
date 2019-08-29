package gophy

import (
	"math"

	"gonum.org/v1/gonum/floats"
)

/*
 This is for calculating likelihoods for nucleotides with multiple models
*/

//PCalcLikePatternsMul parallel caclulation of likelihood with patterns
func PCalcLikePatternsMul(t *Tree, models []*DNAModel, nodemodels map[*Node]int, patternval []float64, wks int) (fl float64) {
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
	fl += math.Log(CalcLikeOneSiteMul(t, models, nodemodels, 0)) * patternval[0]
	for i := 0; i < wks; i++ {
		go CalcLikeWorkMul(t, models, nodemodels, jobs, results)
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

//PCalcLogLikePatternsMul...
func PCalcLogLikePatternsMul(t *Tree, models []*DNAModel, nodemodels map[*Node]int, patternval []float64, wks int) (fl float64) {
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
func CalcLikeOneSiteMul(t *Tree, models []*DNAModel, nodemodels map[*Node]int, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		x := models[nodemodels[n]]
		if len(n.Chs) > 0 {
			CalcLikeNodeMul(n, x, site)
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

//CalcLogLikeOneSiteMul just one site
func CalcLogLikeOneSiteMul(t *Tree, models []*DNAModel, nodemodels map[*Node]int, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		x := models[nodemodels[n]]
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

//CalcLikeWorkMul this is the worker
func CalcLikeWorkMul(t *Tree, models []*DNAModel, nodemodels map[*Node]int, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			x := models[nodemodels[n]]
			if len(n.Chs) > 0 {
				CalcLikeNodeMul(n, x, j)
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

//CalcLogLikeWorkMul this is the worker
func CalcLogLikeWorkMul(t *Tree, models []*DNAModel, nodemodels map[*Node]int, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			x := models[nodemodels[n]]
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

//CalcLikeNodeMul calculate the likelihood of a node
func CalcLikeNodeMul(nd *Node, model *DNAModel, site int) {
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

//PCalcLikePatternsMarkedMul parallel likelihood caclulation with patterns and just update the values
func PCalcLikePatternsMarkedMul(t *Tree, models []*DNAModel, nodemodels map[*Node]int, patternval []float64, wks int) (fl float64) {
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

// CalcLikeOneSiteMarked this uses the marked machinery to recalculate
func CalcLikeOneSiteMarkedMul(t *Tree, models []*DNAModel, nodemodels map[*Node]int, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		x := models[nodemodels[n]]
		if len(n.Chs) > 0 {
			if n.Marked == true {
				CalcLikeNode(n, x, site)
				if n != t.Rt {
					n.Par.Marked = true
				}
			}
		}
		if t.Rt == n && n.Marked == true {
			for i := 0; i < 4; i++ {
				t.Rt.Data[site][i] *= x.BF[i]
			}
			sl = floats.Sum(t.Rt.Data[site])
		} else {
			sl = floats.Sum(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLikeWorkMarked this is intended to calculate only on the marked nodes back to teh root
func CalcLikeWorkMarkedMul(t *Tree, models []*DNAModel, nodemodels map[*Node]int, jobs <-chan int, results chan<- LikeResult) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			x := models[nodemodels[n]]
			if len(n.Chs) > 0 {
				if n.Marked == true {
					CalcLikeNode(n, x, j)
				}
			}
			if t.Rt == n && n.Marked == true {
				for i := 0; i < 4; i++ {
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
func RTMultconditionals(models []*DNAModel, nodemodels map[*Node]int, node *Node, patternval []float64) {
	p := models[nodemodels[node]].GetPCalc(node.Len)
	for s := range patternval {
		for j := 0; j < 4; j++ {
			templike := 0.0
			for k := 0; k < 4; k++ {
				templike += p.At(j, k) * node.TpConds[s][k]
			}
			node.RtConds[s][j] = templike
		}
	}
}

//RVMultconditionals ...
func RVMultconditionals(models []*DNAModel, nodemodels map[*Node]int, node *Node, patternval []float64) {
	p := models[nodemodels[node]].GetPCalc(node.Par.Len)
	for s := range patternval {
		for j := 0; j < 4; j++ {
			node.Par.RvTpConds[s][j] = 0.0
			for k := 0; k < 4; k++ {
				node.Par.RvTpConds[s][j] += p.At(j, k) * node.Par.RvConds[s][k]
			}
		}
	}
}

//CalcLikeFrontBackMult ...
func CalcLikeFrontBackMult(models []*DNAModel, nodemodels map[*Node]int, tree *Tree, patternval []float64) {
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
		TPconditionals(c, patternval)
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
			RVTPconditionals(c, patternval)
		}
	}
}
