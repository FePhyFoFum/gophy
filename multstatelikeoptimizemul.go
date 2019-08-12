package gophy

import (
	"fmt"
	"math"
	"os"

	"gonum.org/v1/gonum/optimize"
)

// OptimizeMS1RMul multistate one rate JC Multiple models reconstruction
func OptimizeMS1RMul(t *Tree, models []StateModel, nodemodels map[*Node]int, patternvals []float64, wks int) {
	count := 0
	//start := time.Now()
	fcn := func(mds []float64) float64 {
		if mds[0] < 0 {
			return 1000000000000
		}
		for i := range models {
			models[i].SetupQJC1Rate(mds[i])
		}
		lnl := PCalcLikePatternsMSMUL(t, models, nodemodels, patternvals, wks)
		count++
		return -lnl
	}
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := make([]float64, len(models)) //}
	for i := range models {
		p0[i] = 1. / float64(models[i].GetNumStates())
	}
	res, err := optimize.Minimize(p, p0, nil, nil)
	if err != nil {
		fmt.Println(err)
	}
	for i := range models {
		models[i].SetupQJC1Rate(res.X[i])
	}
}

// OptimizeMKMSMul optimize "unscaled" GTR
//    symmetrical and scale the last rate to 1
func OptimizeMKMSMul(t *Tree, models []StateModel, nodemodels map[*Node]int, startv float64, patternvals []float64, sym bool, wks int) {
	count := 0
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		for i, j := range models {
			np := ((j.GetNumStates() * j.GetNumStates()) - j.GetNumStates()) / 2
			if sym == false {
				np = ((j.GetNumStates() * j.GetNumStates()) - j.GetNumStates())
			}
			cn := i * np
			j.SetupQMk(mds[cn:cn+np], sym)
		}
		lnl := PCalcLikePatternsMSMUL(t, models, nodemodels, patternvals, wks)
		count++
		return -lnl
	}
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	if models[0].GetNumStates() < 3 {
		fmt.Println("NEED TO DO 2 STATES")
		os.Exit(0)
	}
	p0 := make([]float64, 0)
	for _, j := range models {
		np := ((j.GetNumStates() * j.GetNumStates()) - j.GetNumStates()) / 2
		if sym == false {
			np = ((j.GetNumStates() * j.GetNumStates()) - j.GetNumStates())
		}
		for j := 0; j < np; j++ {
			p0 = append(p0, startv)
		}
	}
	res, err := optimize.Minimize(p, p0, nil, nil)
	if err != nil {
		fmt.Println(err)
	}
	fmt.Println(res.F, res.X)
	for i, j := range models {
		np := ((j.GetNumStates() * j.GetNumStates()) - j.GetNumStates()) / 2
		if sym == false {
			np = ((j.GetNumStates() * j.GetNumStates()) - j.GetNumStates())
		}
		cn := i * np
		j.SetupQMk(res.X[cn:cn+np], sym)
	}
}

//AdjustBLNRMSMul This is a single edge NR
func AdjustBLNRMSMul(node *Node, models []StateModel, nodemodels map[*Node]int, patternvals []float64, t *Tree, wks int, threshold float64) {
	xmin := 10e-8
	xmax := 2.0
	guess := node.Len
	if guess < xmin || guess > xmax {
		guess = 0.1
	}
	startL := PCalcLikePatternsMSMUL(t, models, nodemodels, patternvals, wks)
	x := models[nodemodels[node]]
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
			for j := 0; j < 4; j++ { //TODO: CHANGE THIS FROM 4!!
				for k := 0; k < 4; k++ {
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
	endL := PCalcLikePatternsMSMUL(t, models, nodemodels, patternvals, wks)
	//make sure that we actually made the likelihood better
	if startL > endL {
		node.Len = startLen
	}
	if len(node.Chs) == 2 && (node.Chs[0].Nam == "taxon_1" || node.Chs[1].Nam == "taxon_1") {
		fmt.Println("-")
	}
}

// OptimizeBLNRMSMul Newton-Raphson for each branch. Does 4 passes
func OptimizeBLNRMSMul(t *Tree, models []StateModel, nodemodels map[*Node]int, patternvals []float64, wks int) {
	for _, c := range t.Pre {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBackMSMUL(models, nodemodels, t, patternvals)
		AdjustBLNRMSMul(c, models, nodemodels, patternvals, t, wks, 10e-12)
	}
	for _, c := range t.Post {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBackMSMUL(models, nodemodels, t, patternvals)
		AdjustBLNRMSMul(c, models, nodemodels, patternvals, t, wks, 10e-12)
	}

	for _, c := range t.Pre {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBackMSMUL(models, nodemodels, t, patternvals)
		AdjustBLNRMSMul(c, models, nodemodels, patternvals, t, wks, 10e-12)
	}
	for _, c := range t.Post {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBackMSMUL(models, nodemodels, t, patternvals)
		AdjustBLNRMSMul(c, models, nodemodels, patternvals, t, wks, 10e-12)
	}
}

// OptimizeGTRMSMul optimize GTR
func OptimizeGTRMSMul(t *Tree, models []StateModel, nodemodels map[*Node]int, patternvals []float64, wks int) {
	count := 0
	//start := time.Now()
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		for i, j := range models {
			cn := i * 5
			//fmt.Println(mds[cn:cn+5])
			j.SetScaledRateMatrix(mds[cn:cn+5], true)
			j.SetupQGTR()
		}
		lnl := PCalcLikePatternsMSMUL(t, models, nodemodels, patternvals, wks)
		if count%100 == 0 {
			//curt := time.Now()
			//fmt.Println(count, lnl, curt.Sub(start))
			//start = time.Now()
		}
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := make([]float64, 0)
	for range models {
		for j := 0; j < 5; j++ {
			p0 = append(p0, 1.0)
		}
	}
	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		//fmt.Println(err)
	}
	fmt.Println("   ", res.F)
	for i, j := range models {
		cn := i * 5
		fmt.Println(res.X[cn : cn+5])
		j.SetScaledRateMatrix(res.X[cn:cn+5], true)
		j.SetupQGTR()
	}
}

// OptimizeGTRBPMul optimize GTR and base composition for the different parts
func OptimizeGTRBPMSMul(t *Tree, models []StateModel, nodemodels map[*Node]int, patternvals []float64, wks int) {
	count := 0
	//start := time.Now()
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		for i, j := range models {
			cn := i * 8
			//fmt.Println(mds[cn:cn+5])
			bf := mds[cn+5 : cn+5+3]
			bf = append(bf, 1-SumFloatVec(bf))
			for _, k := range bf {
				if k > 1 || k < 0 {
					return 1000000000000
				}
			}
			j.SetScaledRateMatrix(mds[cn:cn+5], true)
			j.SetBaseFreqs(bf)
			j.SetupQGTR()
		}
		lnl := PCalcLikePatternsMSMUL(t, models, nodemodels, patternvals, wks)
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
	//settings.Concurrent = 0
	//settings.FuncEvaluations = 1000
	//FC := optimize.FunctionConverge{}
	//FC.Relative = 0.001
	//settings.Converger = &FC
	//settings.Recorder = nil
	/*grad := func(grad, x []float64) []float64 {
		return fd.Gradient(grad, fcn, x, nil)
	}*/
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := make([]float64, 0)
	for range models {
		for j := 0; j < 5; j++ {
			p0 = append(p0, 1.0)
		}
		for j := 0; j < 3; j++ {
			p0 = append(p0, 0.25)
		}
	}
	fmt.Println(p0)
	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		//fmt.Println(err)
	}
	fmt.Println("   ", res.F)
	for i, j := range models {
		cn := i * 8
		bf := res.X[cn+5 : cn+5+3]
		bf = append(bf, 1-SumFloatVec(bf))
		fmt.Println(res.X[cn : cn+8])
		j.SetScaledRateMatrix(res.X[cn:cn+5], true)
		j.SetBaseFreqs(bf)
		j.SetupQGTR()
	}
}

// OptimizeGTRMSCompSharedRM optimize GTR base composition but share rate matrix for the different parts
func OptimizeGTRMSCompSharedRM(t *Tree, models []StateModel, nodemodels map[*Node]int, patternvals []float64, wks int) {
	count := 0
	//start := time.Now()
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		for i, j := range models {
			j.SetRateMatrix(mds[0:5])
			cn := (i * 3) + 5
			bf := mds[cn : cn+3]
			bf = append(bf, 1-SumFloatVec(bf))
			for _, k := range bf {
				if k > 1 || k < 0 {
					return 1000000000000
				}
			}
			j.SetBaseFreqs(bf)
			j.SetupQGTR()
		}
		lnl := PCalcLikePatternsMSMUL(t, models, nodemodels, patternvals, wks)
		if count%100 == 0 {
			//curt := time.Now()
			//fmt.Println(count, lnl, curt.Sub(start))
			//start = time.Now()
		}
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	/*grad := func(grad, x []float64) []float64 {
		return fd.Gradient(grad, fcn, x, nil)
	}*/
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := make([]float64, 0)
	for j := 0; j < 5; j++ {
		p0 = append(p0, 1.0)
	}
	for range models {
		for j := 0; j < 3; j++ {
			p0 = append(p0, 0.25)
		}
	}
	fmt.Println(p0)
	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		//fmt.Println(err)
	}
	fmt.Println("   ", res.F)
	for i, j := range models {
		j.SetRateMatrix(res.X[0:5])
		cn := (i * 3) + 5
		bf := res.X[cn : cn+3]
		bf = append(bf, 1-SumFloatVec(bf))
		j.SetBaseFreqs(bf)
		j.SetupQGTR()
	}
}
