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
	for _,x := range models{
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

//CalcLikeOneSiteMul just one site
func CalcLikeOneSiteMul(t *Tree, models []*DNAModel, nodemodels map[*Node]int,site int) float64 {
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

//CalcLikeNode calculate the likelihood of a node
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