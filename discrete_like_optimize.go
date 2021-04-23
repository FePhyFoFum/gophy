package gophy

import (
	"fmt"
	"math"
	"os"
	"time"

	"github.com/go-nlopt/nlopt"
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/optimize"
)

// OptimizeBL This uses the standard gonum optimizers. Not great.
func OptimizeBL(nd *Node, t *Tree, x *DiscreteModel, patternvals []float64, wks int) {
	count := 0
	start := time.Now()
	fcn := func(bl []float64) float64 {
		for _, i := range bl {
			if i < 0 {
				return 1000000000000
			}
		}
		if nd.Len != bl[0] {
			nd.Len = bl[0]
			nd.Marked = true
			if len(nd.Chs) == 0 {
				nd.Par.Marked = true
			}
		}

		//lnl := PCalcLogLike(t, x, nsites, wks)
		//lnl := PCalcLogLikeMarked(t, x, nsites, wks)
		lnl := PCalcLikePatternsMarked(t, x, patternvals, wks)
		for _, j := range t.Post {
			j.Marked = false
		}
		if count%100 == 0 {
			curt := time.Now()
			fmt.Println(count, lnl, curt.Sub(start))
			start = time.Now()
		}
		count++
		return -lnl
	}
	/*grad := func(grad, x []float64) {
		fd.Gradient(grad, fcn, x, nil)
	}*/
	settings := optimize.Settings{} //DefaultSettings()
	settings.MajorIterations = 10
	settings.Concurrent = 0
	settings.FuncEvaluations = 100
	//settings.FunctionThreshold = 0.1
	settings.GradientThreshold = 0.1
	settings.Recorder = nil
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10
	FC.Relative = 10
	FC.Iterations = 10
	//settings.FunctionConverge = &FC
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	p0 = append(p0, nd.Len)
	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		fmt.Println(err)
	}
	nd.Len = res.X[0]
}

// OptimizeGammaBLS optimize all branch lengths
func OptimizeGammaBLS(t *Tree, x *DiscreteModel, patternvals []float64, wks int) {
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

		//lnl := PCalcLogLike(t, x, nsites, wks)
		//lnl := PCalcLogLikeMarked(t, x, nsites, wks)
		lnl := PCalcLikePatternsGamma(t, x, patternvals, wks)
		/*for _, j := range t.Post {
			j.Marked = false
		}*/
		//if count%100 == 0 {
		//	curt := time.Now()
		//	fmt.Println(count, lnl, curt.Sub(start))
		//	start = time.Now()
		//}
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

func OptimizeGammaBLSNL(t *Tree, x *DiscreteModel, patternvals []float64, wks int) {
	//start := time.Now()
	fcn := func(bl, gradient []float64) float64 {
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

		//lnl := PCalcLogLike(t, x, nsites, wks)
		//lnl := PCalcLogLikeMarked(t, x, nsites, wks)
		lnl := PCalcLikePatternsGamma(t, x, patternvals, wks)
		/*for _, j := range t.Post {
			j.Marked = false
		}*/
		//if count%100 == 0 {
		//	curt := time.Now()
		//	fmt.Println(count, lnl, curt.Sub(start))
		//	start = time.Now()
		//}
		return -lnl
	}
	p0 := make([]float64, len(t.Post))
	lbounds := make([]float64, len(p0))
	hbounds := make([]float64, len(p0))
	for i, n := range t.Post {
		p0[i] = n.Len
		lbounds[i] = 0.
		hbounds[i] = 10.
	}

	opt, err := nlopt.NewNLopt(nlopt.LN_AUGLAG, uint(len(p0)))
	opt.SetLowerBounds(lbounds)
	opt.SetUpperBounds(hbounds)
	opt.SetMaxEval(1000)
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

//AdjustBLNR This is a single edge NR
func AdjustBLNR(node *Node, x *DiscreteModel, patternvals []float64, t *Tree, wks int, threshold float64) {
	numstates := x.NumStates
	xmin := 10e-8
	xmax := 2.0
	guess := node.Len
	if guess < xmin || guess > xmax {
		guess = 0.1
	}
	startL := PCalcLogLikePatterns(t, x, patternvals, wks) // TODO: make sure loglike gets the same, was like. should be the same
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
	endL := PCalcLogLikePatterns(t, x, patternvals, wks) // TODO: make sure loglike gets the same, was like. should be the same
	//make sure that we actually made the likelihood better
	if startL > endL || math.IsNaN(endL) {
		node.Len = startLen
	}
}

// OptimizeBLNR Newton-Raphson for each branch. Does 4 passes
func OptimizeBLNR(t *Tree, x *DiscreteModel, patternvals []float64, wks int) {
	for _, c := range t.Pre {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBack(x, t, patternvals)
		AdjustBLNR(c, x, patternvals, t, wks, 10e-12)
	}
	for _, c := range t.Post {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBack(x, t, patternvals)
		AdjustBLNR(c, x, patternvals, t, wks, 10e-12)
	}

	for _, c := range t.Pre {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBack(x, t, patternvals)
		AdjustBLNR(c, x, patternvals, t, wks, 10e-12)
	}
	for _, c := range t.Post {
		if c == t.Rt {
			continue
		}
		CalcLikeFrontBack(x, t, patternvals)
		AdjustBLNR(c, x, patternvals, t, wks, 10e-12)
	}
}

func OptimizeBLNRGN(t *Tree, x *DiscreteModel, patternvals []float64, wks int) {
	return
}

// OptimizeBLS optimize all branch lengths
func OptimizeBLS(t *Tree, x *DiscreteModel, patternvals []float64, wks int) {
	count := 0
	start := time.Now()
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

		lnl := PCalcLogLikePatterns(t, x, patternvals, wks)
		//lnl := PCalcLogLikeMarked(t, x, nsites, wks)
		//lnl := PCalcLikePatternsMarked(t, x, patternvals, wks)
		for _, j := range t.Post {
			j.Marked = false
		}
		if count%100 == 0 {
			curt := time.Now()
			fmt.Println(count, lnl, curt.Sub(start))
			start = time.Now()
		}
		count++
		return -lnl
	}
	/*grad := func(grad, x []float64) {
		fd.Gradient(grad, fcn, x, nil)
	}*/
	settings := optimize.Settings{} //DefaultSettings()
	settings.MajorIterations = 10
	settings.Concurrent = 0
	settings.FuncEvaluations = 10
	//settings.FunctionThreshold = 0.1
	settings.GradientThreshold = 0.1
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

// OptimizeGTRDNA optimize GTR
func OptimizeGTRDNA(t *Tree, x *DiscreteModel, patternvals []float64,
	sup bool, wks int) []float64 {
	var lkfun func(*Tree, *DiscreteModel, []float64, int) float64
	if sup {
		if x.GammaNCats != 0 {
			lkfun = PCalcLogLikePatternsGamma
		} else {
			lkfun = PCalcLogLikePatterns
		}
	} else {
		if x.GammaNCats != 0 {
			lkfun = PCalcLikePatternsGamma
		} else {
			lkfun = PCalcLikePatterns
		}
	}
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		x.SetRateMatrix(mds)
		x.SetupQGTR()
		//lnl := PCalcLikePatterns(t, x, patternvals, wks)
		lnl := lkfun(t, x, patternvals, wks)
		return -lnl
	}
	settings := optimize.Settings{}
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10e-4

	FC.Iterations = 75
	settings.Converger = &FC
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := []float64{1.0, 1.0, 1.0, 1.0, 1.0}
	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		fmt.Println(err)
	}
	//fmt.Println(res.F)
	x.SetRateMatrix(res.X)
	return res.X
}

//OptimizeBF optimizing the basefreq model but for a clade
func OptimizeBF(t *Tree, x *DiscreteModel, patternvals []float64, log bool, wks int) {
	numstates := x.NumStates
	var lkfun func(*Tree, *DiscreteModel, []float64, int) float64
	if log {
		lkfun = PCalcLogLikePatterns
	} else {
		lkfun = PCalcLikePatterns
	}
	count := 0
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
		mds1[numstates-1] = 1. - tsum
		x.SetBaseFreqs(mds1)
		x.SetupQGTR()
		//lnl := PCalcLikePatterns(t, x, patternvals, wks)
		lnl := lkfun(t, x, patternvals, wks)
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10e-3
	FC.Iterations = 75
	settings.Converger = &FC
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	//p0 := []float64{0.25, 0.25, 0.25} // 4-1
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
	x.SetBaseFreqs(mds1)
	x.SetupQGTR()
}

//OptimizeGTRDNASubClade optimizing the GTR model but for a subclade
func OptimizeGTRDNASubClade(t *Tree, n *Node, excl bool, x *DiscreteModel, patternvals []float64, wks int) {
	count := 0
	start := time.Now()
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		x.SetRateMatrix(mds)
		x.SetupQGTR()
		lnl := PCalcLogLikePatternsSubClade(t, n, excl, x, patternvals, wks)
		if count%100 == 0 {
			curt := time.Now()
			fmt.Println(count, lnl, curt.Sub(start))
			start = time.Now()
		}
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	settings.MajorIterations = 10
	settings.Concurrent = 0
	settings.FuncEvaluations = 10
	//settings.FunctionThreshold = 0.1
	settings.GradientThreshold = 0.1
	settings.Recorder = nil
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := []float64{1.0, 1.0, 1.0, 1.0, 1.0}
	res, err := optimize.Minimize(p, p0, nil, nil)
	if err != nil {
		fmt.Println(err)
	}
	fmt.Println(res.F)
	x.SetRateMatrix(res.X)
}

//OptimizeBFSubClade optimizing the basefreq model but for a subclade
func OptimizeBFSubClade(t *Tree, n *Node, excl bool, x *DiscreteModel, patternvals []float64, log bool, wks int) {
	numstates := x.NumStates
	var lkfun func(*Tree, *Node, bool, *DiscreteModel, []float64, int) float64
	if log {
		lkfun = PCalcLogLikePatternsSubClade
	} else {
		lkfun = PCalcLikePatternsSubClade
	}
	count := 0
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
		mds1[numstates-1] = 1. - tsum
		x.SetBaseFreqs(mds1)
		x.SetupQGTR()
		//lnl := PCalcLikePatternsSubClade(t, n, excl, x, patternvals, wks)
		lnl := lkfun(t, n, excl, x, patternvals, wks)
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10e-8
	FC.Relative = 10e-10
	FC.Iterations = 100
	settings.Converger = &FC
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	//p0 := []float64{0.25, 0.25, 0.25} // 4-1
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
	x.SetBaseFreqs(mds1)
	x.SetupQGTR()
}

//OptimizeBFDNARMSubClade optimizing the basefreq model but for a subclade
func OptimizeBFDNARMSubClade(t *Tree, n *Node, excl bool, x *DiscreteModel, patternvals []float64, wks int) {
	count := 0
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		x.SetRateMatrix(mds[0:5])
		bf := []float64{mds[5], mds[6], mds[7]}
		bf = append(bf, 1.-floats.Sum(bf))
		x.SetBaseFreqs(bf)
		x.SetupQGTR()
		lnl := PCalcLikePatternsSubClade(t, n, excl, x, patternvals, wks)
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10e-8
	FC.Relative = 10e-10
	FC.Iterations = 100
	settings.Converger = &FC
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := []float64{1.0, 1.0, 1.0, 1.0, 1.0, 0.25, 0.25, 0.25} // 4-1
	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		fmt.Println(err)
	}
	bf := []float64{res.X[5], res.X[6], res.X[7]}
	bf = append(bf, 1.0-floats.Sum(bf))
	x.SetRateMatrix(res.X[0:5])
	x.SetBaseFreqs(bf)
	x.SetupQGTR()
}

//OptimizeGamma ...
func OptimizeGamma(t *Tree, x *DiscreteModel, patternvals []float64, log bool, wks int) {
	var lkfun func(*Tree, *DiscreteModel, []float64, int) float64
	if log {
		lkfun = PCalcLogLikePatternsGamma
	} else {
		lkfun = PCalcLikePatternsGamma
	}
	count := 0
	fcn := func(mds []float64) float64 {
		if mds[0] < 0 {
			return 1000000000
		}
		x.GammaAlpha = mds[0]
		x.GammaCats = GetGammaCats(x.GammaAlpha, x.GammaNCats, false)
		lnl := lkfun(t, x, patternvals, wks)
		count++
		return -lnl
	}
	settings := optimize.Settings{}
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10e-4
	FC.Iterations = 100
	settings.Converger = &FC
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := []float64{1.0}
	res, err := optimize.Minimize(p, p0, &settings, nil)
	if err != nil {
		fmt.Println(err)
	}
	x.GammaAlpha = res.X[0]
	x.GammaCats = GetGammaCats(x.GammaAlpha, x.GammaNCats, false)
	fmt.Println(res)
}

//OptimizeGammaAndBL ...
func OptimizeGammaAndBL(t *Tree, x *DiscreteModel, patternvals []float64, log bool, wks int) {
	OptimizeGamma(t, x, patternvals, log, wks)
	OptimizeGammaBLSNL(t, x, patternvals, wks)
	OptimizeGamma(t, x, patternvals, log, wks)
	OptimizeGammaBLSNL(t, x, patternvals, wks)
	//OptimizeGamma(t, x, patternvals, log, wks)
	//OptimizeGammaBLSNL(t, x, patternvals, wks)
}

// OptimizeMS1R ...
// hold over from other things, probably change
// this is for specific multistate models
func OptimizeMS1R(t *Tree, x *DiscreteModel, patternvals []float64, wks int) {
	count := 0
	//start := time.Now()
	fcn := func(mds []float64) float64 {
		if mds[0] < 0 {
			return 1000000000000
		}
		x.SetupQJC1Rate(mds[0])
		lnl := PCalcLogLikePatterns(t, x, patternvals, wks)
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
func OptimizeMKMS(t *Tree, x *DiscreteModel, startv float64, patternvals []float64, sym bool, wks int) {
	count := 0
	fcn := func(mds []float64) float64 {
		for _, i := range mds {
			if i < 0 {
				return 1000000000000
			}
		}
		x.SetupQMk(mds, sym)
		lnl := PCalcLogLikePatterns(t, x, patternvals, wks)
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
	//fmt.Println(res.F, res.X)
	x.SetupQMk(res.X, sym)
}
