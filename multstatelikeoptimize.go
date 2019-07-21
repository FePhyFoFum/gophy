package gophy

import (
	"fmt"
	"math"
	"os"

	"gonum.org/v1/gonum/optimize"
)

// OptimizeMULT1R ...
func OptimizeMS1R(t *Tree, x StateModel, patternvals []float64, wks int) {
	count := 0
	//start := time.Now()
	fcn := func(mds []float64) float64 {
		if mds[0] < 0 {
			return 1000000000000
		}
		x.SetupQJC1Rate(mds[0])
		lnl := PCalcLikePatternsMS(t, x, patternvals, wks)
		if count%100 == 0 {
			//curt := time.Now()
			//fmt.Println(count, lnl, curt.Sub(start))
			//start = time.Now()
		}
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	//settings.MajorIterations = 100
	settings.Concurrent = 0
	//settings.FuncEvaluations = 100
	//settings.FunctionThreshold = 0.1
	//settings.GradientThreshold = 0.00001
	settings.Recorder = nil
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := []float64{0.0441} //1. / float64(x.GetNumStates())}
	res, err := optimize.Minimize(p, p0, nil, nil)
	if err != nil {
		fmt.Println(err)
	}
	x.SetupQJC1Rate(res.X[0])
}

// OptimizeMKMS optimize GTR
//    symmetrical and scale the last rate to 1
func OptimizeMKMS(t *Tree, x StateModel, startv float64, patternvals []float64, sym bool, wks int) {
	count := 0
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		x.SetupQMk(mds, sym)
		lnl := PCalcLikePatternsMS(t, x, patternvals, wks)
		count++
		return -lnl
	}
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	if x.GetNumStates() < 3 && sym == true {
		fmt.Println("NEED TO DO 2 STATES")
		os.Exit(0)
	}
	//p0 := make([]float64, (((x.GetNumStates()*x.GetNumStates())-x.GetNumStates())/2)-1) // scaled
	p0 := make([]float64, (((x.GetNumStates() * x.GetNumStates()) - x.GetNumStates()) / 2)) //unscaled
	if sym == false {
		p0 = make([]float64, ((x.GetNumStates() * x.GetNumStates()) - x.GetNumStates())) // unsym
	}
	for i := range p0 {
		p0[i] = startv
	}
	res, err := optimize.Minimize(p, p0, nil, nil)
	if err != nil {
		fmt.Println(err)
	}
	fmt.Println(res.F, res.X)
	x.SetupQMk(res.X, sym)
}

//AdjustBLNRMS This is a single edge NR
func AdjustBLNRMS(node *Node, x StateModel, patternvals []float64, t *Tree, wks int, threshold float64) {
	xmin := 10e-8
	xmax := 2.0
	guess := node.Len
	if guess < xmin || guess > xmax {
		guess = 0.1
	}
	startL := PCalcLikePatternsMS(t, x, patternvals, wks)
	startLen := node.Len
	//need subtree 1
	s1probs := node.TpConds
	//need subtree 2
	s2probs := node.RvConds
	for z := 0; z < 10; z++ {
		t := node.Len
		p := x.GetPCalc(t)
		x.DecomposeQ()
		d1p := x.ExpValueFirstD(t)
		d2p := x.ExpValueSecondD(t)
		d1 := 0.
		d2 := 0.
		like := 1.
		for s := range patternvals {
			templike := 0.
			tempd1 := 0.
			tempd2 := 0.
			for j := 0; j < x.GetNumStates(); j++ {
				for k := 0; k < x.GetNumStates(); k++ {
					templike += (s1probs[s][j] * p.At(j, k) * s2probs[s][k] * x.GetBF()[j])
					tempd1 += (s1probs[s][j] * d1p.At(j, k) * s2probs[s][k] * x.GetBF()[j])
					tempd2 += (s1probs[s][j] * d2p.At(j, k) * s2probs[s][k] * x.GetBF()[j])
				}
			}
			d1 += ((tempd1 / templike) * patternvals[s])
			d2 += (((tempd2 / templike) - (math.Pow(tempd1, 2) / math.Pow(templike, 2))) * patternvals[s])
			like += math.Log(templike)
		}
		if len(node.Chs) == 2 && (node.Chs[0].Nam == "taxon_1" || node.Chs[1].Nam == "taxon_1") {
			//fmt.Println(like, t, t-(d1/d2), d1)
		}
		if (t - (d1 / d2)) < 0 {
			node.Len = 10e-12
			break
		} else {
			node.Len = (t - (d1 / d2))
		}
		if math.Abs(d1) < threshold {
			break
		}
	}
	endL := PCalcLikePatternsMS(t, x, patternvals, wks)
	//make sure that we actually made the likelihood better
	if startL > endL {
		node.Len = startLen
	}
	if len(node.Chs) == 2 && (node.Chs[0].Nam == "taxon_1" || node.Chs[1].Nam == "taxon_1") {
		fmt.Println("-")
	}
}

// OptimizeBLNRMS Newton-Raphson for each branch. Does 4 passes
func OptimizeBLNRMS(t *Tree, x StateModel, patternvals []float64, wks int) {
	for _, c := range t.Pre {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBackMS(x, t, patternvals)
		AdjustBLNRMS(c, x, patternvals, t, wks, 10e-12)
	}
	for _, c := range t.Post {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBackMS(x, t, patternvals)
		AdjustBLNRMS(c, x, patternvals, t, wks, 10e-12)
	}

	for _, c := range t.Pre {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBackMS(x, t, patternvals)
		AdjustBLNRMS(c, x, patternvals, t, wks, 10e-12)
	}
	for _, c := range t.Post {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBackMS(x, t, patternvals)
		AdjustBLNRMS(c, x, patternvals, t, wks, 10e-12)
	}
}
