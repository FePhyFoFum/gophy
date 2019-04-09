package gophy

import (
	"math"

	"gonum.org/v1/gonum/floats"
)

/*
 This is for calculating likelihoods for nucleotides
*/

// PCalcLogLikeMS this will calculate log like in parallel
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

//PCalcLikeMS parallel calculate likelihood
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

//PCalcLikePatterns parallel caclulation of likelihood with patterns
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

//TODO: this whole bit needs to be checcked

//PCalcLikePatternsMarked parallel likelihood caclulation with patterns and just update the values
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

//PCalcLogLikePatterns parallel log likeliohood calculation including patterns
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
	return
}

// PCalcLogLikeBack a bit of a shortcut. Could do better, but walks back from the n node to the root
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

//PCalcLogLikeMarked parallel calculation of loglike with just updating
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

// CalcLogLikeOneSite just calculate the likelihood of one site
// probably used to populate the PDict in the DNA Model so that we can reuse the calculations
func CalcLogLikeOneSiteMS(t *Tree, x StateModel, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			CalcLogLikeNodeMS(n, x, site)
		}
		if t.Rt == n {
			for i := 0; i < 4; i++ {
				t.Rt.Data[site][i] += math.Log(x.GetBF()[i])
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
		}
	}
	return sl
}

//CalcLikeOneSiteMS just one site
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

// CalcLogLikeOneSiteBackMS like the one above but from nb to the root only
func CalcLogLikeOneSiteBackMS(t *Tree, nb *Node, x StateModel, site int) float64 {
	sl := 0.0
	going := true
	cur := nb
	for going {
		if len(cur.Chs) > 0 {
			CalcLogLikeNodeMS(cur, x, site)
		}
		if cur == t.Rt {
			for i := 0; i < 4; i++ {
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

// CalcLogLikeOneSiteMarkedMS this uses the marked machinery to recalculate
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
			for i := 0; i < 4; i++ {
				t.Rt.Data[site][i] += math.Log(x.GetBF()[i])
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
		} else {
			sl = floats.LogSumExp(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLikeOneSiteMarkedMS this uses the marked machinery to recalculate
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
			for i := 0; i < 4; i++ {
				t.Rt.Data[site][i] *= x.GetBF()[i]
			}
			sl = floats.Sum(t.Rt.Data[site])
		} else {
			sl = floats.Sum(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLogLikeWorkMS this is intended for a worker that will be executing this per site
func CalcLogLikeWorkMS(t *Tree, x StateModel, jobs <-chan int, results chan<- LikeResult) { //results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				CalcLogLikeNodeMS(n, x, j)
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

//CalcLikeWorkMS this is the worker
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

// CalcLogLikeWorkBack this is intended for a worker that will be executing this per site
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
				for i := 0; i < 4; i++ {
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

// CalcLikeWorkMarked this is intended to calculate only on the marked nodes back to teh root
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
				for i := 0; i < 4; i++ {
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

// CalcLogLikeWorkMarked this is intended to calculate only on the marked nodes back to teh root
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
				for i := 0; i < 4; i++ {
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

// CalcLogLikeNode calculates likelihood for node
func CalcLogLikeNodeMS(nd *Node, model StateModel, site int) {
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

// CalcLikeFrontBack ...
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
		//for _, c := range tree.Tips {
		//	fmt.Println(c.Nam, c.TpConds[i])
		//}
		//fmt.Println("p", i, patternval[i])
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
					//fmt.Print((s*x.GetBF()[j])/su, " ")
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
					//fmt.Print(s/su, " ")
					retstates[c][i][j] = s / su
				}
				//fmt.Print("\n")
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
		//for _, c := range tree.Tips {
		//	fmt.Println(c.Nam, c.TpConds[i])
		//}
		//fmt.Println("p", i, patternval[i])
		for _, c := range tree.Pre {
			if len(c.Chs) == 0 {
				//	continue
			}
			//fmt.Println(c.Newick(true))
			retstates[c][i] = make([]float64, x.GetNumStates())
			if c == tree.Rt {
				/*
					su := 0.
					for j, s := range c.RtConds[i] {
						su += (s * x.GetBF()[j])
					}
					for j, s := range c.RtConds[i] {
						//fmt.Print((s*x.GetBF()[j])/su, " ")
						retstates[c][i][j] = (s * x.GetBF()[j]) / su
					}
				*/
				//fmt.Print("\n")
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
				//fmt.Println(tvN)
				//fmt.Println(tvR)
				if time { //time
					retstates[c][i] = tvR
				} else { //number
					retstates[c][i] = tvN
				}
				//su := 0.
				//for _, s := range tv {
				//	su += s
				//}
				//for j, s := range tv {
				//	//fmt.Print(s/su, " ")
				//	retstates[c][i][j] = s / su
				//}
				//fmt.Print("\n")
			}
		}
	}
	return
}
