package gophy

import (
	"math"

	"gonum.org/v1/gonum/mathext"

	"gonum.org/v1/gonum/stat/distuv"

	"gonum.org/v1/gonum/floats"
)

//LikeResult likelihood value and site
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

//PCalcLike parallel calculate likelihood
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

//PCalcLikePatterns parallel caclulation of likelihood with patterns
func PCalcLikePatterns(t *Tree, x *DNAModel, patternval []float64, wks int) (fl float64) {
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

//TODO: this whole bit needs to be checcked

//PCalcLikePatternsMarked parallel likelihood caclulation with patterns and just update the values
func PCalcLikePatternsMarked(t *Tree, x *DNAModel, patternval []float64, wks int) (fl float64) {
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

//PCalcLogLikePatterns parallel log likeliohood calculation including patterns
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

//PCalcLogLikeMarked parallel calculation of loglike with just updating
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

//CalcLikeOneSite just one site
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

// CalcLikeOneSiteMarked this uses the marked machinery to recalculate
func CalcLikeOneSiteMarked(t *Tree, x *DNAModel, site int) float64 {
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

//CalcLikeWork this is the worker
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

// CalcLikeWorkMarked this is intended to calculate only on the marked nodes back to teh root
func CalcLikeWorkMarked(t *Tree, x *DNAModel, jobs <-chan int, results chan<- LikeResult) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
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

//CalcLikeNode calculate the likelihood of a node
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

func pointGamma(p float64, a float64, b float64) float64 {
	c := distuv.ChiSquared{K: 2 * (a)}
	return ((c.Quantile(p)) / (2 * b))
}

func GetGammaCats(alpha float64, cats int, median bool) []float64 {
	K := float64(cats)
	a, b := alpha, alpha
	factor := a / b * K
	rK := make([]float64, cats)
	if median == true {
		gap05 := 1.0 / (2.0 * K)
		for i := 0; i < cats; i++ {
			rK[i] = pointGamma((float64(i)*2.0+1.)*gap05, a, b)
		}
		t := 0.
		for i := 0; i < cats; i++ {
			t += rK[i]
		}
		for i := 0; i < cats; i++ {
			rK[i] *= factor / t
		}
		return rK
	} else {
		freqK := make([]float64, cats)
		for i := 0; i < cats-1; i++ {
			freqK[i] = pointGamma((float64(i)+1.0)/K, a, b)
		}
		for i := 0; i < cats-1; i++ {
			freqK[i] = mathext.GammaIncReg(a+1, freqK[i]*b)
		}
		rK[0] = freqK[0] * factor
		rK[cats-1] = (1 - freqK[cats-2]) * factor
		for i := 1; i < cats-1; i++ {
			rK[i] = (freqK[i] - freqK[i-1]) * factor
		}
		return rK
	}
}

/*
 * calculate the conditionals for ancestral calc or branch lengths

 toward tip
tpcond  X
        | | rvtpcond
        | ^
        | |
        v |
        | | rvcond
rtcond  x
 toward root
*/
func TPconditionals(node *Node, patternval []float64) {
	if len(node.Chs) > 0 {
		for s := range patternval {
			for j := 0; j < 4; j++ {
				node.TpConds[s][j] = 1.
				for _, i := range node.Chs {
					node.TpConds[s][j] *= i.RtConds[s][j]
				}
			}
		}
	}
}

func RTconditionals(x *DNAModel, node *Node, patternval []float64) {
	p := x.GetPCalc(node.Len)
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

func RVconditionals(x *DNAModel, node *Node, patternval []float64) {
	p := x.GetPCalc(node.Par.Len)
	for s := range patternval {
		for j := 0; j < 4; j++ {
			node.Par.RvTpConds[s][j] = 0.0
			for k := 0; k < 4; k++ {
				node.Par.RvTpConds[s][j] += p.At(j, k) * node.Par.RvConds[s][k]
			}
		}
	}
}

func RVTPconditionals(node *Node, patternval []float64) {
	for s := range patternval {
		for j := 0; j < 4; j++ {
			node.RvConds[s][j] = node.Par.RvTpConds[s][j]
		}
		for _, oc := range node.Par.Chs {
			if node == oc {
				continue
			}
			for j := 0; j < 4; j++ {
				node.RvConds[s][j] *= oc.RtConds[s][j]
			}
		}
	}
}

// CalcLikeFrontBack ...
func CalcLikeFrontBack(x *DNAModel, tree *Tree, patternval []float64) {
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
			RVTPconditionals(c, patternval)
		}
	}
}

// CalcAncStates for each node based on the calculations above
func CalcAncStates(x *DNAModel, tree *Tree, patternval []float64) (retstates map[*Node][][]float64) {
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
			retstates[c][i] = []float64{0.0, 0.0, 0.0, 0.0}
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
				tv := []float64{0.0, 0.0, 0.0, 0.0}
				for j := 0; j < 4; j++ {
					for k := 0; k < 4; k++ {
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
