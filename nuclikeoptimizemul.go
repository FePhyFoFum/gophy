package gophy

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/optimize"
)

// OptimizeBLSMul optimize all branch lengths
func OptimizeBLSMul(t *Tree, models []*DNAModel, nodemodels map[*Node]int, patternvals []float64, wks int) float64 {
	count := 0
	//start := time.Now()
	fcn := func(bl []float64) float64 {
		for _, i := range bl {
			if i < 0 {
				return 1000000000000
			}
		}
		for x, n := range t.Post {
			if n.Len != bl[x] {
				n.Len = bl[x]
				n.Marked = true
				if len(n.Chs) == 0 {
					n.Par.Marked = true
				}
			}
		}

		lnl := PCalcLikePatternsMarkedMul(t, models, nodemodels, patternvals, wks)
		for _, j := range t.Post {
			j.Marked = false
		}
		if count%100 == 0 {
			//curt := time.Now()
			//fmt.Println(count, lnl, curt.Sub(start))
			//start = time.Now()
		}
		count++
		return -lnl
	}

	settings := optimize.Settings{}
	//settings.MajorIterations = 1000
	//settings.Concurrent = 0
	//settings.FuncEvaluations = 10000
	//FC := optimize.FunctionConverge{}
	//FC.Absolute = 0.1
	//settings.Converger = &FC
	//settings.Recorder = nil
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	for _, n := range t.Post {
		p0 = append(p0, n.Len)
	}
	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		fmt.Println(err)
	}
	//fmt.Println("  ", res.F)
	for x, n := range t.Post {
		n.Len = res.X[x]
	}
	return res.F
}

//AdjustBLNRMult This is a single edge NR
func AdjustBLNRMult(node *Node, models []*DNAModel, nodemodels map[*Node]int, patternvals []float64, t *Tree, wks int, threshold float64) {
	xmin := 10e-8
	xmax := 2.0
	guess := node.Len
	if guess < xmin || guess > xmax {
		guess = 0.1
	}
	startL := PCalcLikePatternsMul(t, models, nodemodels, patternvals, wks)
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
			for j := 0; j < 4; j++ {
				for k := 0; k < 4; k++ {
					templike += (s1probs[s][j] * p.At(j, k) * s2probs[s][k] * x.BF[j])
					tempd1 += (s1probs[s][j] * d1p.At(j, k) * s2probs[s][k] * x.BF[j])
					tempd2 += (s1probs[s][j] * d2p.At(j, k) * s2probs[s][k] * x.BF[j])
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
	endL := PCalcLikePatterns(t, x, patternvals, wks)
	//make sure that we actually made the likelihood better
	if startL > endL {
		node.Len = startLen
	}
	if len(node.Chs) == 2 && (node.Chs[0].Nam == "taxon_1" || node.Chs[1].Nam == "taxon_1") {
		fmt.Println("-")
	}
}

// OptimizeBLNRMult Newton-Raphson for each branch. Does 4 passes
func OptimizeBLNRMult(t *Tree, models []*DNAModel, nodemodels map[*Node]int, patternvals []float64, wks int) {
	for _, c := range t.Pre {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBackMult(models, nodemodels, t, patternvals)
		AdjustBLNRMult(c, models, nodemodels, patternvals, t, wks, 10e-6)
	}
	for _, c := range t.Post {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBackMult(models, nodemodels, t, patternvals)
		AdjustBLNRMult(c, models, nodemodels, patternvals, t, wks, 10e-12)
	}
	for _, c := range t.Pre {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBackMult(models, nodemodels, t, patternvals)
		AdjustBLNRMult(c, models, nodemodels, patternvals, t, wks, 10e-12)
	}
	for _, c := range t.Post {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBackMult(models, nodemodels, t, patternvals)
		AdjustBLNRMult(c, models, nodemodels, patternvals, t, wks, 10e-20)
	}
}

// OptimizeGTRMul optimize GTR
func OptimizeGTRMul(t *Tree, models []*DNAModel, nodemodels map[*Node]int, patternvals []float64, wks int) {
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
			j.SetRateMatrix(mds[cn : cn+5])
			j.SetupQGTR()
		}
		lnl := PCalcLikePatternsMul(t, models, nodemodels, patternvals, wks)
		if count%100 == 0 {
			//curt := time.Now()
			//fmt.Println(count, lnl, curt.Sub(start))
			//start = time.Now()
		}
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	//settings.MajorIterations = 10000
	//settings.Concurrent = 0
	//settings.FuncEvaluations = 1000
	//FC := optimize.FunctionConverge{}
	//FC.Relative = 0.1
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
	}
	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		//fmt.Println(err)
	}
	fmt.Println("   ", res.F)
	for i, j := range models {
		cn := i * 5
		fmt.Println(res.X[cn : cn+5])
		j.SetRateMatrix(res.X[cn : cn+5])
		j.SetupQGTR()
	}
}

// OptimizeGTRBPMul optimize GTR and base composition for the different parts
func OptimizeGTRBPMul(t *Tree, models []*DNAModel, nodemodels map[*Node]int, patternvals []float64, wks int) {
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
			j.SetRateMatrix(mds[cn : cn+5])
			j.SetBaseFreqs(bf)
			j.SetupQGTR()
		}
		lnl := PCalcLikePatternsMul(t, models, nodemodels, patternvals, wks)
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
		j.SetRateMatrix(res.X[cn : cn+5])
		j.SetBaseFreqs(bf)
		j.SetupQGTR()
	}
}

// OptimizeGTRCompSharedRM optimize GTR base composition but share rate matrix for the different parts
func OptimizeGTRCompSharedRM(t *Tree, models []*DNAModel, nodemodels map[*Node]int, patternvals []float64, wks int) {
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
		lnl := PCalcLikePatternsMul(t, models, nodemodels, patternvals, wks)
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
