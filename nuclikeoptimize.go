package gophy

import (
	"fmt"
	"time"

	"gonum.org/v1/gonum/optimize"
)

// OptimizeBL ...
func OptimizeBL(nd *Node, t *Tree, x *DNAModel, nsites int, wks int) {
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
		lnl := PCalcLogLikeMarked(t, x, nsites, wks)
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
	settings := optimize.DefaultSettingsLocal()
	//settings.UseInitialData = false
	settings.FunctionThreshold = 0.01
	//settings.GradientThreshold = 0.01
	settings.Concurrent = 0
	//settings.Recorder = nil
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10
	FC.Relative = 10
	FC.Iterations = 10
	settings.FunctionConverge = &FC
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	p0 = append(p0, nd.Len)
	res, err := optimize.Minimize(p, p0, settings, nil)
	if err != nil {
		//fmt.Println(err)
	}
	nd.Len = res.X[0]
}
