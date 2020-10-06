package clustering

import (
	"fmt"

	"github.com/FePhyFoFum/gophy"
	"gonum.org/v1/gonum/optimize"
)

//OptimizePreservationLam will optimize the poisson rate parameter in the preservation model
func OptimizePreservationLam(tree *gophy.Tree) (float64, float64) {
	//lam := 1.0 //2.4
	preNodes := tree.Pre
	fcn := func(p []float64) float64 {
		lam := p[0]
		large := 100000000000.0
		if lam <= 0.0 {
			return large
		}
		stratLL := ADPoissonTreeLoglike(preNodes, lam)
		lnl := stratLL
		return -lnl
	}
	settings := optimize.Settings{} //DefaultSettings()
	settings.MajorIterations = 10
	settings.Concurrent = 0
	settings.FuncEvaluations = 100
	settings.GradientThreshold = 0.1
	settings.Recorder = nil
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	p0 := []float64{1.0}
	meth := &optimize.NelderMead{}
	res, err := optimize.Minimize(p, p0, nil, meth)
	if err != nil {
		fmt.Println(err)
	}
	return -res.F, res.X[0]

}

/*
func OptimizeMorphStratHeights(tree *Node, lam float64) (float64, float64, []float64) {
	//lam := 1.0 //2.4
	fcn := func(heights []float64) float64 {
		large := 100000000000.0
		for _, i := range heights {
			if i <= 0.0 {
				return large
			}
		}
		preNodes := tree.PreorderArray()
		bad := AssignInternalNodeHeights(preNodes, heights)
		if bad {
			return large
		}
		morphLL := RootedLogLikeParallel(tree, true, 4)
		stratLL := ADPoissonTreeLoglike(preNodes, lam)
		lnl := morphLL + stratLL
		return -lnl
	}
	settings := optimize.Settings{} //DefaultSettings()
	settings.MajorIterations = 10
	settings.Concurrent = 0
	settings.FuncEvaluations = 100
	settings.GradientThreshold = 0.1
	settings.Recorder = nil
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	preNodes := tree.PreorderArray()
	for _, n := range preNodes {
		if n.ISTIP == false { // && n.ANC == false {
			p0 = append(p0, n.Height)
		}
	}
	meth := &optimize.NelderMead{}
	res, err := optimize.Minimize(p, p0, nil, meth)
	if err != nil {
		fmt.Println(err)
	}
	AssignInternalNodeHeights(preNodes, res.X)
	var retparams []float64
	retparams = append(retparams, lam)
	for _, bl := range res.X {
		retparams = append(retparams, bl)
	}
	return -res.F, float64(len(res.X)), retparams //res.X

}

func OptimizeLamMorphStratHeights(tree *Node) (float64, float64, []float64) {
	fcn := func(params []float64) float64 {
		lam := params[0]
		heights := params[1:]
		large := 100000000000.0
		for _, i := range heights {
			if i <= 0.0 {
				return large
			}
		}
		preNodes := tree.PreorderArray()
		bad := AssignInternalNodeHeights(preNodes, heights)
		if bad {
			return large
		}
		morphLL := RootedLogLikeParallel(tree, true, 4)
		stratLL := ADPoissonTreeLoglike(preNodes, lam)
		lnl := morphLL + stratLL
		return -lnl
	}
	settings := optimize.Settings{} //DefaultSettings()
	settings.MajorIterations = 10
	settings.Concurrent = 0
	settings.FuncEvaluations = 100
	settings.GradientThreshold = 0.1
	settings.Recorder = nil
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	p0 = append(p0, 1.0)
	preNodes := tree.PreorderArray()
	for _, n := range preNodes {
		if n.ISTIP == false { // && n.ANC == false {
			p0 = append(p0, n.Height)
		}
	}
	meth := &optimize.NelderMead{}
	res, err := optimize.Minimize(p, p0, nil, meth)
	if err != nil {
		fmt.Println(err)
	}
	AssignInternalNodeHeights(preNodes, res.X[1:])
	return -res.F, float64(len(res.X)), res.X

}

func OptimizeGlobalRateHeights(tree *Node, lam float64) (float64, []float64) {
	fcn := func(params []float64) float64 {
		rate := params[0]
		heights := params[1:]
		large := 100000000000.0
		for _, i := range heights {
			if i <= 0.0 {
				return large
			}
		}
		preNodes := tree.PreorderArray()
		AssignGlobalRate(preNodes, rate)
		bad := AssignInternalNodeHeights(preNodes, heights)
		if bad {
			return large
		}
		morphLL := RootedLogLikeParallel(tree, true, 4)
		stratLL := ADPoissonTreeLoglike(preNodes, lam)
		lnl := morphLL + stratLL
		return -lnl
	}
	settings := optimize.Settings{} //DefaultSettings()
	settings.MajorIterations = 10
	settings.Concurrent = 0
	settings.FuncEvaluations = 100
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
	p0 = append(p0, 0.5)
	preNodes := tree.PreorderArray()
	for _, n := range preNodes {
		if len(n.Chs) > 0 { //&& n.ANC == false {
			p0 = append(p0, n.Height)
		}
	}
	meth := &optimize.NelderMead{}
	res, err := optimize.Minimize(p, p0, nil, meth)
	if err != nil {
		fmt.Println(err)
	}
	AssignGlobalRate(preNodes, res.X[0])
	AssignInternalNodeHeights(preNodes, res.X[1:])
	return -res.F, res.X
}
*/
