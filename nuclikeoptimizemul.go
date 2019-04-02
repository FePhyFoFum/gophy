package gophy

import (
	"fmt"

	"gonum.org/v1/gonum/diff/fd"

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

// OptimizeGTR optimize GTR
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
	grad := func(grad, x []float64) []float64 {
		return fd.Gradient(grad, fcn, x, nil)
	}
	p := optimize.Problem{Func: fcn, Grad: grad, Hess: nil}
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
	grad := func(grad, x []float64) []float64 {
		return fd.Gradient(grad, fcn, x, nil)
	}
	p := optimize.Problem{Func: fcn, Grad: grad, Hess: nil}
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
