package gophy

import (
	"fmt"
	"math"
	"time"

	"github.com/go-nlopt/nlopt"
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/optimize"
)

// OptimizeBLSMul optimize all branch lengths
func OptimizeBLSMul(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int, patternvals []float64, wks int) float64 {
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
func AdjustBLNRMult(node *Node, models []*DiscreteModel, nodemodels map[*Node]int, patternvals []float64, t *Tree, wks int, threshold float64) {
	xmin := 10e-8
	xmax := 2.0
	guess := node.Len
	if guess < xmin || guess > xmax {
		guess = 0.1
	}
	startL := PCalcLogLikePatternsMul(t, models, nodemodels, patternvals, wks)
	x := models[nodemodels[node]]
	numstates := x.NumStates
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
			for j := 0; j < numstates; j++ {
				for k := 0; k < numstates; k++ {
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
	endL := PCalcLogLikePatternsMul(t, models, nodemodels, patternvals, wks)
	//make sure that we actually made the likelihood better
	if startL > endL {
		node.Len = startLen
	}
	if len(node.Chs) == 2 && (node.Chs[0].Nam == "taxon_1" || node.Chs[1].Nam == "taxon_1") {
		fmt.Println("-")
	}
}

// OptimizeBLNRMult Newton-Raphson for each branch. Does 4 passes
func OptimizeBLNRMult(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int, patternvals []float64, wks int) {
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

func OptimizeGammaMult(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int,
	patternvals []float64, log bool, wks int) {
	var lkfun func(*Tree, []*DiscreteModel, map[*Node]int, []float64, int) float64
	if log {
		lkfun = PCalcLogLikePatternsGammaMul
	} else {
		lkfun = PCalcLikePatternsGammaMul
	}
	count := 0
	fcn := func(mds []float64) float64 {
		if mds[0] < 0 {
			return 1000000000
		}
		for i := range mds {
			models[i].GammaAlpha = mds[i]
			models[i].GammaCats = GetGammaCats(models[i].GammaAlpha, models[i].GammaNCats, false)
		}
		lnl := lkfun(t, models, nodemodels, patternvals, wks)
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10e-5
	FC.Iterations = 100
	settings.Converger = &FC
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := make([]float64, len(models))
	for i := range p0 {
		p0[i] = 1.0
	}
	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		fmt.Println(err)
	}
	for i := range res.X {
		models[i].GammaAlpha = res.X[i]
		models[i].GammaCats = GetGammaCats(models[i].GammaAlpha, models[i].GammaNCats, false)
	}
	fmt.Println(res)
}

// OptimizeGammaBLS optimize all branch lengths
func OptimizeGammaBLSMult(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int, patternvals []float64, wks int) {
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
				/*n.Marked = true
				if len(n.Chs) == 0 {
					n.Par.Marked = true
				}*/
			}
		}
		lnl := PCalcLikePatternsGammaMul(t, models, nodemodels, patternvals, wks)
		count++
		return -lnl
	}
	/*grad := func(grad, x []float64) {
		fd.Gradient(grad, fcn, x, nil)
	}*/
	settings := optimize.Settings{} //DefaultSettings()
	//settings.MajorIterations = 10
	settings.Concurrent = 0
	//settings.FuncEvaluations = 10
	//settings.GradientThreshold = 0.1
	settings.Recorder = nil
	//FC := optimize.FunctionConverge{}
	//FC.Absolute = 10
	//FC.Relative = 10
	//FC.Iterations = 10
	//settings.FunctionConverge = &FC
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	for _, n := range t.Post {
		p0 = append(p0, n.Len)
	}
	res, err := optimize.Minimize(p, p0, nil, nil)
	if err != nil {
		fmt.Println(err)
	}
	fmt.Println(res.F)
	for x, n := range t.Post {
		n.Len = res.X[x]
	}
}

func OptimizeGammaBLSNLMult(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int,
	patternvals []float64, log bool, wks int) {
	//start := time.Now()
	var lkfun func(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int,
		patternval []float64, wks int) (fl float64)
	if log {
		lkfun = PCalcLogLikePatternsGammaMul
	} else {
		lkfun = PCalcLikePatternsGammaMul
	}
	fcn := func(bl, gradient []float64) float64 {
		for _, i := range bl {
			if i < 0 {
				return 1000000000000
			}
		}
		for x, n := range t.Post {
			if n.Len != bl[x] {
				n.Len = bl[x]
			}
		}
		lnl := lkfun(t, models, nodemodels, patternvals, wks)
		return -lnl
	}
	p0 := make([]float64, len(t.Post))
	lbounds := make([]float64, len(p0))
	hbounds := make([]float64, len(p0))
	for i, n := range t.Post {
		p0[i] = n.Len
		lbounds[i] = 0.
		hbounds[i] = 5.
	}

	opt, err := nlopt.NewNLopt(nlopt.LN_AUGLAG, uint(len(p0)))
	opt.SetLowerBounds(lbounds)
	opt.SetUpperBounds(hbounds)
	opt.SetMaxEval(500)
	opt.SetFtolAbs(10e-5)
	opt.SetMinObjective(fcn)
	xopt, minf, err := opt.Optimize(p0)
	if err != nil {
		fmt.Println(err)
	}
	fmt.Println(minf)
	for x, n := range t.Post {
		n.Len = xopt[x]
	}
}

func OptimizeGammaAndBLMult(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int,
	patternvals []float64, log bool, wks int) {
	OptimizeGammaMult(t, models, nodemodels, patternvals, log, wks)
	OptimizeGammaBLSNLMult(t, models, nodemodels, patternvals, log, wks)
	OptimizeGammaMult(t, models, nodemodels, patternvals, log, wks)
	OptimizeGammaBLSNLMult(t, models, nodemodels, patternvals, log, wks)
	//OptimizeGammaMult(t, models, nodemodels, patternvals, log, wks)
	//OptimizeGammaBLSNLMult(t, models, nodemodels, patternvals, log, wks)
}

func OptimizeSharedGammaAndBLMult(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int,
	patternvals []float64, log bool, wks int) {
	OptimizeGammaBLSNLMult(t, models, nodemodels, patternvals, log, wks)
	OptimizeGammaBLSNLMult(t, models, nodemodels, patternvals, log, wks)
	//OptimizeGammaBLSNLMult(t, models, nodemodels, patternvals, log,wks)
}

// OptimizeGTRDNAMul optimize GTR
func OptimizeGTRDNAMul(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int, patternvals []float64, wks int) {
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

// OptimizeGTRBPDNAMul optimize GTR and base composition for the different parts
func OptimizeGTRBPDNAMul(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int, usemodelvals bool,
	patternvals []float64, log bool, wks int) {
	var lkfun func(*Tree, []*DiscreteModel, map[*Node]int, []float64, int) float64
	if log {
		lkfun = PCalcLogLikePatternsMul
	} else {
		lkfun = PCalcLikePatternsMul
	}
	count := 0
	//start := time.Now()
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		cur := 0
		for _, j := range models {
			j.SetRateMatrix(mds[cur : cur+5])
			cur += 5
			bf := []float64{mds[cur], mds[cur+1], mds[cur+2]}
			cur += 3
			bf = append(bf, 1-SumFloatVec(bf))
			for _, k := range bf {
				if k > 1 || k < 0 {
					return 1000000000000
				}
			}

			j.SetBaseFreqs(bf)
			j.SetupQGTR()
		}
		lnl := lkfun(t, models, nodemodels, patternvals, wks)
		if count%100 == 0 {
			//curt := time.Now()
			//fmt.Println(count, lnl, curt.Sub(start))
			//start = time.Now()
		}
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10e-8
	FC.Relative = 10e-8
	FC.Iterations = 100
	settings.Converger = &FC
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := make([]float64, 0)
	if usemodelvals == false {
		for range models {
			for j := 0; j < 5; j++ {
				p0 = append(p0, 1.0)
			}
			for j := 0; j < 3; j++ {
				p0 = append(p0, 0.25)
			}
		}
	} else {
		for i := range models {
			p0 = append(p0, models[i].R.At(0, 1))
			p0 = append(p0, models[i].R.At(0, 2))
			p0 = append(p0, models[i].R.At(0, 3))
			p0 = append(p0, models[i].R.At(1, 2))
			p0 = append(p0, models[i].R.At(1, 3))
			for j := 0; j < 3; j++ {
				p0 = append(p0, models[i].BF[j])
			}
		}
	}

	fmt.Println(p0)
	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		//fmt.Println(err)
	}
	fmt.Println("   ", res.F)
	cur := 0
	for _, j := range models {
		j.SetRateMatrix(res.X[cur : cur+5])
		cur += 5
		bf := []float64{res.X[cur], res.X[cur+1], res.X[cur+2]}
		cur += 3
		bf = append(bf, 1-SumFloatVec(bf))
		j.SetBaseFreqs(bf)
		j.SetupQGTR()
	}
}

// OptimizeGTRDNACompSharedRM optimize GTR base composition but share rate matrix for the different parts
func OptimizeGTRDNACompSharedRM(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int,
	usemodelvals bool, patternvals []float64, log bool, wks int) {
	var lkfun func(*Tree, []*DiscreteModel, map[*Node]int, []float64, int) float64
	if log {
		lkfun = PCalcLogLikePatternsMul
	} else {
		lkfun = PCalcLikePatternsMul
	}
	count := 0
	//start := time.Now()
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		cur := 5
		//fmt.Println(" ++", mds)
		for _, j := range models {
			j.SetRateMatrix(mds[0:5])
			bf := []float64{mds[cur], mds[cur+1], mds[cur+2]}
			cur += 3
			bf = append(bf, 1-floats.Sum(bf))
			for _, k := range bf {
				if k > 1 || k < 0 {
					return 1000000000000
				}
			}
			j.SetBaseFreqs(bf)
			j.SetupQGTR()
		}
		//lnl := PCalcLikePatternsMul(t, models, nodemodels, patternvals, wks)
		lnl := lkfun(t, models, nodemodels, patternvals, wks)

		if count%100 == 0 {
			//curt := time.Now()
			//fmt.Println(count, lnl, curt.Sub(start))
			//start = time.Now()
		}
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10e-5 //10e-10
	FC.Relative = 10e-6
	FC.Iterations = 200
	settings.Converger = &FC
	/*grad := func(grad, x []float64) []float64 {
		return fd.Gradient(grad, fcn, x, nil)
	}*/
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := make([]float64, 0)
	if usemodelvals == false {
		for j := 0; j < 5; j++ {
			p0 = append(p0, 1.0)
		}
		for range models {
			for j := 0; j < 3; j++ {
				p0 = append(p0, 0.25)
			}
		}
	} else {
		p0 = append(p0, models[0].R.At(0, 1))
		p0 = append(p0, models[0].R.At(0, 2))
		p0 = append(p0, models[0].R.At(0, 3))
		p0 = append(p0, models[0].R.At(1, 2))
		p0 = append(p0, models[0].R.At(1, 3))
		cur := 5
		for i := range models {
			for j := 0; j < 3; j++ {
				p0 = append(p0, models[i].BF[j])
			}
			cur += 3
		}
	}

	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		//fmt.Println(err)
	}
	//fmt.Println("   ", res.F)
	cur := 5
	for _, j := range models {
		j.SetRateMatrix(res.X[0:5])
		bf := []float64{res.X[cur], res.X[cur+1], res.X[cur+2]}
		cur += 3
		bf = append(bf, 1-floats.Sum(bf))
		j.SetBaseFreqs(bf)
		j.SetupQGTR()
	}
}

func OptimizeGTRDNACompSharedRMNL(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int,
	usemodelvals bool, patternvals []float64, log bool, wks int) {
	var lkfun func(*Tree, []*DiscreteModel, map[*Node]int, []float64, int) float64
	if log {
		lkfun = PCalcLogLikePatternsMul
	} else {
		lkfun = PCalcLikePatternsMul
	}
	count := 0
	//start := time.Now()
	fcn := func(mds, gradient []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		cur := 5
		//fmt.Println(" ++", mds)
		for _, j := range models {
			j.SetRateMatrix(mds[0:5])
			bf := []float64{mds[cur], mds[cur+1], mds[cur+2]}
			cur += 3
			bf = append(bf, 1-floats.Sum(bf))
			for _, k := range bf {
				if k > 1 || k < 0 {
					return 1000000000000
				}
			}
			j.SetBaseFreqs(bf)
			j.SetupQGTR()
		}
		lnl := lkfun(t, models, nodemodels, patternvals, wks)

		if count%100 == 0 {
			//curt := time.Now()
			//fmt.Println(count, lnl)
			//start = time.Now()
		}
		count++
		return -lnl
	}
	p0 := make([]float64, 0)
	if usemodelvals == false {
		for j := 0; j < 5; j++ {
			p0 = append(p0, 1.0)
		}
		for range models {
			for j := 0; j < 3; j++ {
				p0 = append(p0, 0.25)
			}
		}
	} else {
		p0 = append(p0, models[0].R.At(0, 1))
		p0 = append(p0, models[0].R.At(0, 2))
		p0 = append(p0, models[0].R.At(0, 3))
		p0 = append(p0, models[0].R.At(1, 2))
		p0 = append(p0, models[0].R.At(1, 3))
		cur := 5
		for i := range models {
			for j := 0; j < 3; j++ {
				p0 = append(p0, models[i].BF[j])
			}
			cur += 3
		}
	}
	opt, err := nlopt.NewNLopt(nlopt.LN_AUGLAG, uint(len(p0)))
	opt.SetMaxEval(1000)
	opt.SetFtolAbs(10e-5)
	opt.SetMinObjective(fcn)
	res, _, err := opt.Optimize(p0)
	if err != nil {
		//fmt.Println(err)
	}
	//fmt.Println("   ", res.F)
	cur := 5
	for _, j := range models {
		j.SetRateMatrix(res[0:5])
		bf := []float64{res[cur], res[cur+1], res[cur+2]}
		cur += 3
		bf = append(bf, 1-floats.Sum(bf))
		j.SetBaseFreqs(bf)
		j.SetupQGTR()
	}
}

// OptimizeAACompSharedRM optimize GTR base composition but share rate matrix for the different parts
func OptimizeAACompSharedRM(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int,
	patternvals []float64, log bool, wks int) {
	var lkfun func(*Tree, []*DiscreteModel, map[*Node]int, []float64, int) float64
	if log {
		lkfun = PCalcLogLikePatternsMul
	} else {
		lkfun = PCalcLikePatternsMul
	}
	count := 0
	start := time.Now()
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		cur := 0
		//fmt.Println(" ++", mds)
		for _, j := range models {
			bf := make([]float64, 20)
			for k := 0; k < 19; k++ {
				bf[k] = mds[cur]
				cur++
			}
			bf[19] = 1 - floats.Sum(bf[0:19])
			j.SetBaseFreqs(bf)
			j.SetupQGTR()
		}
		//lnl := PCalcLikePatternsMul(t, models, nodemodels, patternvals, wks)
		lnl := lkfun(t, models, nodemodels, patternvals, wks)

		if count%100 == 0 {
			curt := time.Now()
			fmt.Println(count, lnl, curt.Sub(start))
			start = time.Now()
		}
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10e-8 //10e-10
	FC.Relative = 10e-10
	FC.Iterations = 100
	settings.Converger = &FC
	/*grad := func(grad, x []float64) []float64 {
		return fd.Gradient(grad, fcn, x, nil)
	}*/
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := make([]float64, 0)

	cur := 0
	for i := range models {
		for j := 0; j < 19; j++ { // 20 -1
			p0 = append(p0, models[i].BF[j])
			cur++
		}
	}

	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		//fmt.Println(err)
	}
	//fmt.Println("   ", res.F)
	cur = 0
	for _, j := range models {
		bf := make([]float64, 19)
		for j := 0; j < 19; j++ {
			bf[j] = res.X[j]
		}
		bf = append(bf, 1-floats.Sum(bf))
		j.SetBaseFreqs(bf)
		j.SetupQGTR()
	}
}

func OptimizeAACompSharedRMNLOPT(t *Tree, models []*DiscreteModel, nodemodels map[*Node]int,
	patternvals []float64, log bool, wks int) {
	var lkfun func(*Tree, []*DiscreteModel, map[*Node]int, []float64, int) float64
	if log {
		lkfun = PCalcLogLikePatternsMul
	} else {
		lkfun = PCalcLikePatternsMul
	}
	count := 0
	start := time.Now()
	fcn := func(mds, gradient []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		cur := 0
		fmt.Println(" ++", mds)
		for _, j := range models {
			bf := make([]float64, 20)
			for k := 0; k < 19; k++ {
				bf[k] = mds[cur]
				cur++
			}
			fs := floats.Sum(bf[0:19])
			if fs > 1 {
				return 1000000000000
			}
			bf[19] = 1 - fs
			j.SetBaseFreqs(bf)
			j.SetupQGTR()
		}
		//lnl := PCalcLikePatternsMul(t, models, nodemodels, patternvals, wks)
		lnl := lkfun(t, models, nodemodels, patternvals, wks)

		if count%100 == 0 {
			curt := time.Now()
			fmt.Println(count, lnl, curt.Sub(start))
			start = time.Now()
		}
		count++
		return -lnl
	}
	p0 := make([]float64, 0)
	cur := 0
	for i := range models {
		for j := 0; j < 19; j++ { // 20 -1
			p0 = append(p0, models[i].BF[j])
			cur++
		}
	}
	opt, err := nlopt.NewNLopt(nlopt.LN_AUGLAG, uint(len(p0)))
	opt.SetMaxEval(2000)
	opt.SetLowerBounds1(10e-8)
	opt.SetUpperBounds1(0.9999)
	opt.SetFtolAbs(10e-5)
	//opt.AddInequalityConstraint(confcn, 10e-8)
	opt.SetMinObjective(fcn)
	res, minf, err := opt.Optimize(p0)
	if err != nil {
		//fmt.Println(err)
	}
	fmt.Println("   ", minf)
	cur = 0
	for _, j := range models {
		bf := make([]float64, 19)
		for j := 0; j < 19; j++ {
			bf[j] = res[j]
		}
		bf = append(bf, 1-floats.Sum(bf))
		j.SetBaseFreqs(bf)
		j.SetupQGTR()
	}
}

// OptimizeGTRDNACompSharedRMSingleModel optimize GTR base composition but share rate matrix for the different parts
//   you send an int that will identify which model can be adjusted and only that one will be
func OptimizeGTRDNACompSharedRMSingleModel(t *Tree, models []*DiscreteModel,
	nodemodels map[*Node]int, usemodelvals bool, whichmodel int,
	patternvals []float64, log bool, wks int) {
	var lkfun func(*Tree, []*DiscreteModel, map[*Node]int, []float64, int) float64
	if log {
		lkfun = PCalcLogLikePatternsMul
	} else {
		lkfun = PCalcLikePatternsMul
	}
	count := 0
	//start := time.Now()
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		cur := 5
		//fmt.Println(" ++", mds)
		j := models[whichmodel]
		j.SetRateMatrix(mds[0:5])
		bf := []float64{mds[cur], mds[cur+1], mds[cur+2]}
		cur += 3
		bf = append(bf, 1-floats.Sum(bf))
		for _, k := range bf {
			if k > 1 || k < 0 {
				return 1000000000000
			}
		}
		j.SetBaseFreqs(bf)
		j.SetupQGTR()
		//lnl := PCalcLikePatternsMul(t, models, nodemodels, patternvals, wks)
		lnl := lkfun(t, models, nodemodels, patternvals, wks)
		if count%100 == 0 {
			//curt := time.Now()
			//fmt.Println(count, lnl, curt.Sub(start))
			//start = time.Now()
		}
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10e-8 //10e-10
	FC.Relative = 10e-10
	FC.Iterations = 100
	settings.Converger = &FC
	/*grad := func(grad, x []float64) []float64 {
		return fd.Gradient(grad, fcn, x, nil)
	}*/
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := make([]float64, 0)
	if usemodelvals == false {
		for j := 0; j < 5; j++ {
			p0 = append(p0, 1.0)
		}
		for j := 0; j < 3; j++ {
			p0 = append(p0, 0.25)
		}
	} else {
		p0 = append(p0, models[0].R.At(0, 1))
		p0 = append(p0, models[0].R.At(0, 2))
		p0 = append(p0, models[0].R.At(0, 3))
		p0 = append(p0, models[0].R.At(1, 2))
		p0 = append(p0, models[0].R.At(1, 3))
		cur := 5
		for j := 0; j < 3; j++ {
			p0 = append(p0, models[whichmodel].BF[j])
		}
		cur += 3
	}

	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		//fmt.Println(err)
	}
	//fmt.Println("   ", res.F)
	cur := 5
	j := models[whichmodel]
	j.SetRateMatrix(res.X[0:5])
	bf := []float64{res.X[cur], res.X[cur+1], res.X[cur+2]}
	cur += 3
	bf = append(bf, 1-floats.Sum(bf))
	j.SetBaseFreqs(bf)
	j.SetupQGTR()
}

// OptimizeGTRDNACompSharedRMSubClade optimize GTR base composition but share rate matrix for the different parts
func OptimizeGTRDNACompSharedRMSubClade(t *Tree, n *Node, excl bool, models []*DiscreteModel, nodemodels map[*Node]int,
	usemodelvals bool, patternvals []float64, wks int) {
	count := 0
	//start := time.Now()
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		cur := 5
		for _, j := range models {
			j.SetRateMatrix(mds[0:5])
			bf := []float64{mds[cur], mds[cur+1], mds[cur+2]}
			cur += 3
			bf = append(bf, 1-SumFloatVec(bf))
			for _, k := range bf {
				if k > 1 || k < 0 {
					return 1000000000000
				}
			}
			j.SetBaseFreqs(bf)
			j.SetupQGTR()
		}
		lnl := PCalcLikePatternsMulSubClade(t, n, excl, models, nodemodels, patternvals, wks)
		if count%100 == 0 {
			//curt := time.Now()
			//fmt.Println(count, lnl, curt.Sub(start))
			//start = time.Now()
		}
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10e-10
	FC.Iterations = 100
	settings.Converger = &FC
	/*grad := func(grad, x []float64) []float64 {
		return fd.Gradient(grad, fcn, x, nil)
	}*/
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := make([]float64, 0)
	if usemodelvals == false {
		for j := 0; j < 5; j++ {
			p0 = append(p0, 1.0)
		}
		for range models {
			for j := 0; j < 3; j++ {
				p0 = append(p0, 0.25)
			}
		}
	} else {
		p0 = append(p0, models[0].R.At(0, 1))
		p0 = append(p0, models[0].R.At(0, 2))
		p0 = append(p0, models[0].R.At(0, 3))
		p0 = append(p0, models[0].R.At(1, 2))
		p0 = append(p0, models[0].R.At(1, 3))
		for i := range models {
			for j := 0; j < 3; j++ {
				p0 = append(p0, models[i].BF[j])
			}
		}
	}
	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		//fmt.Println(err)
	}
	//fmt.Println("   ", res.F)
	cur := 5
	for _, j := range models {
		j.SetRateMatrix(res.X[0:5])
		bf := []float64{res.X[cur], res.X[cur+1], res.X[cur+2]}
		cur += 3
		bf = append(bf, 1-floats.Sum(bf))
		j.SetBaseFreqs(bf)
		j.SetupQGTR()
	}
}

//OptimizeAACompSharedRMSingleModel optimize AA base comp but share rate matrix for different parts
func OptimizeAACompSharedRMSingleModel(t *Tree, models []*DiscreteModel,
	nodemodels map[*Node]int, usemodelvals bool, whichmodel int,
	patternvals []float64, log bool, wks int) {
	var lkfun func(*Tree, []*DiscreteModel, map[*Node]int, []float64, int) float64
	if log {
		lkfun = PCalcLogLikePatternsMul
	} else {
		lkfun = PCalcLikePatternsMul
	}
	numstates := models[whichmodel].NumStates
	count := 0
	//start := time.Now()
	fcn := func(mds []float64) float64 {
		mds1 := make([]float64, numstates)
		tsum := 0.
		for j, i := range mds {
			if i < 0 {
				return 1000000000000
			}
			mds1[j] = i
			tsum += i
		}
		if tsum > 1.0 {
			return 1000000000000
		}
		j := models[whichmodel]
		mds1[numstates-1] = 1. - tsum
		j.SetBaseFreqs(mds1)
		j.SetupQGTR()
		//lnl := PCalcLikePatternsMul(t, models, nodemodels, patternvals, wks)
		lnl := lkfun(t, models, nodemodels, patternvals, wks)
		if count%100 == 0 {
			//curt := time.Now()
			//fmt.Println(count, lnl, curt.Sub(start))
			//start = time.Now()
		}
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10e-8 //10e-10
	FC.Relative = 10e-10
	FC.Iterations = 100
	settings.Converger = &FC

	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := make([]float64, numstates-1)
	for i := range p0 {
		p0[i] = 1. / float64(numstates)
	}
	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		fmt.Println(err)
	}
	mds1 := make([]float64, numstates)
	tsum := 0.
	for j, i := range res.X {
		mds1[j] = i
		tsum += i
	}
	mds1[numstates-1] = 1. - tsum
	j := models[whichmodel]
	j.SetBaseFreqs(mds1)
	j.SetupQGTR()
}
