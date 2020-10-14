package gophy

//
//TODO: add
//      CV
//      free nodes/ fossil nodes

import (
	"fmt"
	"math"
	"math/rand"
	"os"
	"runtime"

	"github.com/go-nlopt/nlopt"
	"gonum.org/v1/gonum/optimize"
)

// PLObj things for PL
type PLObj struct {
	Rates                []float64
	Durations            []float64
	Dates                []float64
	Mins                 []float64
	PenMins              []float64
	PenMaxs              []float64
	Maxs                 []float64
	CharDurations        []float64
	LogFactCharDurations []float64
	ParentsNdsInts       []int
	ChildrenVec          [][]int
	FreeNodes            []int
	FreeNodesM           map[int]bool
	NumNodes             int
	LogPen               bool
	PenaltyBoundary      float64
	Smoothing            float64
	RateGroups           map[int]int //[nodenum]rategroup
	CVNode               int
	Tree                 Tree
	Logs                 map[float64]float64
}

// OptimizeRD optimization
func (p *PLObj) OptimizeRD(params []float64, likeFunc func([]float64, bool) float64,
	LF bool, PL bool) *optimize.Result {
	count := 0
	//start := time.Now()
	fcn := func(pa []float64) float64 {
		for _, i := range pa {
			if i < 0 {
				return 1000000000000
			}
		}
		pl := likeFunc(pa, true)
		if pl == -1 {
			return 100000000000
		}
		if (count+1)%1000 == 0 {
			//os.Exit(0)
			//curt := time.Now()
			//fmt.Fprintln(os.Stderr, count, pl, curt.Sub(start))
			//start = time.Now()
		}
		count++
		return pl
	}
	var grad func(grad, x []float64)
	if LF == true {
		grad = func(grad, x []float64) {
			//fmt.Println("before", grad)
			p.calcLFGradient(x, &grad)
			//fmt.Println("after", grad)
			//os.Exit(0)

		}
	} else if PL == true {
		grad = func(grad, x []float64) {
			//fmt.Println("before", grad)
			p.calcPLGradient(x, &grad)
			//fmt.Println("after", grad)
			//os.Exit(0)

		}
	} else {
		grad = nil
	}
	/*grad := func(grad, x []float64) {
		fd.Gradient(grad, fcn, x, nil)
	}*/
	settings := optimize.Settings{}
	FC := optimize.FunctionConverge{}
	FC.Absolute = 10e-3
	//FC.Relative = 1
	FC.Iterations = 1000
	settings.Converger = &FC
	ps := optimize.Problem{Func: fcn, Grad: grad, Hess: nil}
	var p0 []float64
	p0 = params
	if LF || PL {
		res, err := optimize.Minimize(ps, p0, &settings, &optimize.LBFGS{})
		if err != nil {
			//fmt.Println(err)
		}
		return res
	}
	res, err := optimize.Minimize(ps, p0, &settings, nil)
	if err != nil {
		//fmt.Println(err)
	}
	return res
	//fmt.Println(res.F)
	//fmt.Println(res)
}

//TODO: getback to this
func (p *PLObj) OptimizeRDBOUNDED(params []float64) {
	var p0 []float64
	p0 = params

	opt, err := nlopt.NewNLopt(nlopt.LD_LBFGS, uint(len(params)))
	if err != nil {
		panic(err)
	}
	defer opt.Destroy()

	//get bounds from mins and maxs
	lbounds := make([]float64, len(params))
	hbounds := make([]float64, len(params))
	lbounds[0] = 0
	hbounds[0] = 100000000
	c := 1
	for _, i := range p.FreeNodes {
		lbounds[c] = p.Mins[i]
		hbounds[c] = p.Maxs[i]
		c++
	}
	opt.SetLowerBounds(lbounds)
	opt.SetUpperBounds(hbounds)
	var evals int
	fcn := func(pa, gradient []float64) float64 {
		evals++
		fmt.Println(pa)
		for _, i := range pa {
			if i < 0 {
				return 1000000000000
			}
		}
		pl := p.CalcLF(pa, true)
		fmt.Println(pl)
		if pl == -1 {
			return 100000000000
		}
		p.calcLFGradient(pa, &gradient)
		return pl
	}
	opt.SetMinObjective(fcn)
	xopt, minf, err := opt.Optimize(p0)
	fmt.Println(xopt, minf)
	/*
		optimizer := new(lbfgsb.Lbfgsb)
		optimizer.Init(len(params))
		optimizer.SetBoundsSparse(bounds)
		optimizer.SetFTolerance(1e-10)
		optimizer.SetGTolerance(1e-10)
		x := new(LikeFunction)
		x.p = p
		obj := x
		minimum, exitStatus := optimizer.Minimize(obj, p0)
		fmt.Println(minimum.F, minimum.X, minimum.G, exitStatus)
	*/
}

type LikeFunction struct {
	p *PLObj
}

func (sf LikeFunction) EvaluateFunction(pa []float64) float64 {
	for _, i := range pa {
		if i < 0 {
			return 1000000000000
		}
	}

	pl := sf.p.CalcLF(pa, true)
	fmt.Println(pa, pl)
	if pl == -1 {
		return 100000000000
	}
	return pl
}

func (sf LikeFunction) EvaluateGradient(point []float64) (gradient []float64) {
	gradient = make([]float64, len(point))
	sf.p.calcLFGradient(point, &gradient)
	return

}

// CalcPL ...
func (p *PLObj) CalcPL(params []float64, free bool) float64 {
	for _, i := range params {
		if i < 0 {
			return -1.
		}
	}

	pcount := 0
	fcount := 0
	for i := range p.Rates {
		if i == 0 { // skip the root
			continue
		}
		p.Rates[i] = params[pcount]
		pcount++
		fcount++
	}
	if free { // only set the free nodes
		for _, i := range p.FreeNodes {
			p.Dates[i] = params[pcount]
			pcount++
			fcount++
		}
	} else {
		for i := range p.Dates {
			p.Dates[i] = params[pcount]
			pcount++
			fcount++
		}
	}

	ret := p.SetDurations()
	if ret == false {
		return -1.
	}
	/*
		jobs := make(chan int, p.NumNodes-1)
		results := make(chan float64, p.NumNodes-1)
		for w := 0; w < 5; w++ {
			go p.PCalcRateLogLike(jobs, results)
		}
		njobs := 0
		for i := 1; i < p.NumNodes; i++ {
			jobs <- i
			njobs++
		}
		close(jobs)
		ll := 0.0
		for i := 0; i < njobs; i++ {
			ll += <-results
		}
	*/
	ll := p.CalcRateLogLike()
	tp := p.CalcPenalty()
	rp := p.CalcRoughnessPenalty()

	pl := (ll + tp) + p.Smoothing*rp
	return pl
}

//CalcLF ...
func (p *PLObj) CalcLF(params []float64, free bool) float64 {
	for _, i := range params {
		if i < 0 {
			return -1.
		}
	}

	pcount := 0
	fcount := 0
	for i := range p.Rates {
		p.Rates[i] = params[pcount]
		fcount++
	}
	pcount++
	if free { // only set the free nodes
		for _, i := range p.FreeNodes {
			p.Dates[i] = params[pcount]
			pcount++
			fcount++
		}
	} else {
		for i := range p.Dates {
			p.Dates[i] = params[pcount]
			pcount++
			fcount++
		}
	}
	ret := p.SetDurations()
	if ret == false {
		return -1.
	}
	/*
		jobs := make(chan int, p.NumNodes-1)
		results := make(chan float64, p.NumNodes-1)
		for w := 0; w < 3; w++ {
			go p.PCalcRateLogLike(jobs, results)
		}
		njobs := 0
		for i := 1; i < p.NumNodes; i++ {
			jobs <- i
			njobs++
		}
		close(jobs)
		ll := 0.0
		for i := 0; i < njobs; i++ {
			ll += <-results
		}
	*/
	ll := p.CalcRateLogLike()
	//ll := p.CalcRateLogLikeP()
	return ll
}

//CalcMultLF ...
func (p *PLObj) CalcMultLF(params []float64, free bool) float64 {
	for _, i := range params {
		if i < 0 {
			return -1.
		}
	}

	pcount := 0
	fcount := 0
	for i := range p.Rates {
		p.Rates[i] = params[p.RateGroups[i]]
		fcount++
		if pcount < p.RateGroups[i] {
			pcount = p.RateGroups[i]
		}
	}
	pcount++
	if free { // only set the free nodes
		for _, i := range p.FreeNodes {
			p.Dates[i] = params[pcount]
			pcount++
			fcount++
		}
	} else {
		for i := range p.Dates {
			p.Dates[i] = params[pcount]
			pcount++
			fcount++
		}
	}
	ret := p.SetDurations()
	if ret == false {
		return -1.
	}
	ll := p.CalcRateLogLike()
	//fmt.Println(ll)
	return ll
}

// CalcMultPL ...
func (p *PLObj) CalcMultPL(params []float64, free bool) float64 {
	for _, i := range params {
		if i < 0 {
			return -1.
		}
	}

	pcount := 0
	fcount := 0
	for i := range p.Rates {
		if i == 0 { // skip the root
			continue
		}
		p.Rates[i] = params[pcount]
		pcount++
		fcount++
	}
	if free { // only set the free nodes
		for _, i := range p.FreeNodes {
			p.Dates[i] = params[pcount]
			pcount++
			fcount++
		}
	} else {
		for i := range p.Dates {
			p.Dates[i] = params[pcount]
			pcount++
			fcount++
		}
	}
	ret := p.SetDurations()
	if ret == false {
		return -1.
	}

	ll := p.CalcRateLogLike()
	tp := p.CalcPenalty()
	rp := p.CalcRoughnessPenaltyMult()
	pl := (ll + tp) + p.Smoothing*rp
	return pl
}

func minLength(l float64, numsites float64) float64 {
	return math.Max(l, 1./numsites)
}

// SetValues ...
func (p *PLObj) SetValues(t Tree, numsites float64, minmap map[*Node]float64,
	maxmap map[*Node]float64, verbose bool) {
	p.Tree = t
	p.CVNode = -1
	p.LogPen = false
	p.NumNodes = 0
	p.FreeNodes = make([]int, 0)
	p.FreeNodesM = make(map[int]bool)
	//assign node numbers
	for i, n := range t.Pre {
		n.Num = i
		p.NumNodes++
		if n != t.Rt && len(n.Chs) > 0 {
			p.FreeNodes = append(p.FreeNodes, n.Num)
			p.FreeNodesM[n.Num] = true
		}
	}
	p.CharDurations = make([]float64, p.NumNodes)
	p.Durations = make([]float64, p.NumNodes)
	p.LogFactCharDurations = make([]float64, p.NumNodes)
	p.ParentsNdsInts = make([]int, p.NumNodes)
	p.ChildrenVec = make([][]int, p.NumNodes)
	p.Dates = make([]float64, p.NumNodes)
	p.Rates = make([]float64, p.NumNodes)
	for i, n := range t.Pre {
		if n.Par != nil {
			p.ParentsNdsInts[i] = n.Par.Num
		}
		p.ChildrenVec[i] = make([]int, len(n.Chs))
		for x, m := range n.Chs {
			p.ChildrenVec[i][x] = m.Num
		}
		p.CharDurations[i] = math.Round(minLength(n.Len, numsites) * numsites)
		p.LogFactCharDurations[i] = LogFact(p.CharDurations[i])
	}
	// set min maxs
	p.Maxs = make([]float64, p.NumNodes)
	p.Mins = make([]float64, p.NumNodes)
	p.PenMins = make([]float64, p.NumNodes)
	p.PenMaxs = make([]float64, p.NumNodes)
	for _, n := range t.Post {
		//if len(n.Chs) > 0 { //internal -- taken out for the tip dates
		if _, ok := minmap[n]; ok {
			p.Mins[n.Num] = minmap[n]
		} else {
			ymin := 0.0
			for _, j := range n.Chs {
				if len(j.Chs) > 0 {
					if p.Mins[j.Num] > ymin {
						ymin = p.Mins[j.Num]
					}
				}
			}
			p.Mins[n.Num] = ymin
		}
		if _, ok := maxmap[n]; ok {
			p.Maxs[n.Num] = maxmap[n]
		} else {
			ymax := 10000.0
			par := n
			for {
				par = par.Par
				tmax := ymax
				if _, ok := maxmap[par]; ok {
					tmax = maxmap[par]
				}
				if tmax < ymax {
					ymax = tmax
				}
				if par.Par == nil {
					break
				}
			}
			p.Maxs[n.Num] = ymax
		}
	}
	//fmt.Println(p.Mins)
	//fmt.Println(p.Maxs)
	//os.Exit(0)
	//end set min max

	//setup dates
	rtDt := minmap[t.Rt]
	p.Dates[t.Rt.Num] = rtDt
	setMaxNodesToPresent(t)
	for _, n := range t.Rt.Chs {
		p.feasibleTimes(n, rtDt)
		//p.feasibleTimesSplit(n, rtDt) //just splitting them, can go back to the other
	}
	durcheck := p.SetDurations()
	if verbose {
		fmt.Println(p.PrintNewickDurations(t) + ";")
	}
	if durcheck == false {
		fmt.Println(p.Durations)
		os.Exit(0)
	}
}

//RunLF langley fitch run with starting float, probably numsites/20.
func (p *PLObj) RunLF(startRate float64, verbose bool) *optimize.Result {
	params := make([]float64, len(p.FreeNodes)+1) // lf
	c := 0
	params[c] = startRate //lf
	c++                   //lf
	for _, i := range p.FreeNodes {
		params[c] = p.Dates[i]
		c++
	}
	//fmt.Fprintln(os.Stderr, params)
	val := p.CalcLF(params, true)
	if verbose {
		fmt.Fprintln(os.Stderr, "start LF:", val)
	}
	x := p.OptimizeRD(params, p.CalcLF, true, false)
	if verbose {
		fmt.Fprintln(os.Stderr, "end LF:", x.F)
	}
	//p.OptimizeRDBOUNDED(params)
	return x
}

//RunMLF multiple rate langley fitch with starting float, probably x.X[0] from LF
// and mrcagroup
func (p *PLObj) RunMLF(startRate float64, mrcagroups []*Node, t Tree,
	verbose bool) *optimize.Result {
	p.RateGroups = make(map[int]int)
	params := make([]float64, len(p.FreeNodes)+len(mrcagroups)+1)
	params[0] = startRate
	for i := range mrcagroups {
		params[i+1] = startRate
	}
	p.SetRateGroups(t, mrcagroups)
	c := len(mrcagroups) + 1
	for _, i := range p.FreeNodes {
		params[c] = p.Dates[i]
		c++
	}
	//fmt.Fprintln(os.Stderr, params)
	val := p.CalcMultLF(params, true)
	if verbose {
		fmt.Fprintln(os.Stderr, "start MLF:", val)
	}
	x := p.OptimizeRD(params, p.CalcMultLF, false, false)
	if verbose {
		fmt.Fprintln(os.Stderr, "end MLF:", x.F)
	}
	return x
}

//RunPL penalized likelihood run with a starting float, probably x.X[0] from LF
func (p *PLObj) RunPL(startRate float64, verbose bool) *optimize.Result {
	params := make([]float64, p.NumNodes-1+len(p.FreeNodes)) // pl
	c := 0
	for i := range p.Rates { //every edge has a rate
		if i == 0 { //root
			continue
		}
		params[c] = startRate
		c++
	}
	for _, i := range p.FreeNodes {
		params[c] = p.Dates[i]
		c++
	}
	val := p.CalcPL(params, true)
	if verbose {
		fmt.Fprintln(os.Stderr, "start PL:", val)
	}
	x := p.OptimizeRD(params, p.CalcPL, false, true)
	if verbose {
		fmt.Fprintln(os.Stderr, "end PL:", x.F)
	}
	return x
}

//RunMPL penalized likelihood run with a starting float from x.X[0] from LF
// and mrcagroup -- assuming that p.RateGroups has already been setup
func (p *PLObj) RunMPL(mrcagroups []*Node, t Tree, verbose bool) *optimize.Result {
	params := make([]float64, p.NumNodes-1+len(p.FreeNodes))
	c := 0
	//fmt.Fprintln(os.Stderr, p.Rates)
	for i, j := range p.Rates { //every edge has a rate
		if i == 0 { //root
			continue
		}
		params[c] = j
		c++
	}
	for _, i := range p.FreeNodes {
		params[c] = p.Dates[i]
		c++
	}
	val := p.CalcMultPL(params, true)
	if verbose {
		fmt.Fprintln(os.Stderr, "start MPL:", val)
	}
	x := p.OptimizeRD(params, p.CalcMultPL, false, true)
	if verbose {
		fmt.Fprintln(os.Stderr, "end MPL:", x.F)
	}
	return x
}

//RunCV ...
// this is just doing cross validation LOOCV
func (p *PLObj) RunCV(params []float64, verbose bool, likeFunc func([]float64, bool) float64) {
	chisq := 0.0
	for _, i := range p.Tree.Tips {
		p.CVNode = i.Num
		p.OptimizeRD(params, likeFunc, false, false)
		//fmt.Println(x.F)
		ratees := 0.
		par := i.Par.Num
		if par == 0 { //is the root
			ratees = p.Rates[i.GetSib().Num]
		} else {
			ratees = p.Rates[par]
		}
		d := p.Durations[p.CVNode]
		expe := ratees * d
		length := p.CharDurations[p.CVNode]
		sq := (length - expe) * (length - expe) // (expe-length)*(expe-length);
		if d < 1e-100 || expe < 1e-100 {
			chisq += 0
		} else {
			chisq += (sq / expe) //average chisq
		}
		if verbose {
			fmt.Fprintln(os.Stderr, i.Nam, " rate:", ratees, " expe:", expe, " dur:",
				d, " len: ", length, " sq: ", sq, " chisq: ", chisq)
		}
	}
	fmt.Fprintln(os.Stderr, "chisq", chisq)
}

//SetRateGroups send the mrcagroups that will identify the rate shifts
// these go preorder so have the deepest first if they are nested
func (p *PLObj) SetRateGroups(tree Tree, mrcagroups []*Node) {
	for i := 1; i < p.NumNodes; i++ {
		p.RateGroups[i] = 0
	}
	count := 1
	for _, nd := range mrcagroups {
		for _, i := range nd.PreorderArray() {
			p.RateGroups[i.Num] = count
		}
		count++
	}
}

func setMaxNodesToPresent(tree Tree) {
	for _, i := range tree.Tips {
		i.FData["maxnodes"] = 1
	}
	for _, i := range tree.Post {
		if len(i.Chs) > 0 {
			mx := 0.0
			for _, j := range i.Chs {
				mx = math.Max(mx, j.FData["maxnodes"])
			}
			i.FData["maxnodes"] = mx + 1
		}
	}
}

//split the difference
func (p *PLObj) feasibleTimesSplit(nd *Node, timeAnc float64) {
	if len(nd.Chs) == 0 {
		p.Dates[nd.Num] = 0
		if _, ok := nd.FData["tipdates"]; ok {
			p.Dates[nd.Num] = nd.FData["tipdates"]
		}
		return
	}
	depth := nd.FData["maxnodes"]
	p.Dates[nd.Num] = timeAnc - (timeAnc / depth)
	//fmt.Println(nd, timeAnc, depth, timeAnc/depth)
	for _, n := range nd.Chs {
		p.feasibleTimesSplit(n, p.Dates[nd.Num])
	}
}

func (p *PLObj) feasibleTimes(nd *Node, timeAnc float64) {
	if len(nd.Chs) == 0 {
		p.Dates[nd.Num] = 0
		if _, ok := nd.FData["tipdates"]; ok {
			p.Dates[nd.Num] = nd.FData["tipdates"]
		}
		return
	}
	//fmt.Println(nd, nd.Num)
	if p.Maxs[nd.Num] < timeAnc {
		timeAnc = p.Maxs[nd.Num]
	}
	p.Dates[nd.Num] = timeAnc - (timeAnc-p.Mins[nd.Num])*(0.02*rand.Float64()*0.96)/math.Log(2.+3.)
	//fmt.Println(p.Dates[nd.Num], timeAnc)
	for _, n := range nd.Chs {
		p.feasibleTimes(n, p.Dates[nd.Num])
	}
}

// CalcPenalty ...
func (p *PLObj) CalcPenalty() (rkt float64) {
	tpen := 0.0
	for i := 0; i < p.NumNodes; i++ {
		if p.PenMins[i] != -1 {
			tpen += 1. / (p.Dates[i] - p.PenMins[i])
			if math.IsInf(tpen, 1) {
				tpen = 0
			}
		}
		if p.PenMaxs[i] != -1 {
			tpen += 1. / (p.PenMaxs[i] - p.Dates[i])
			if math.IsInf(tpen, 1) {
				tpen = 0
			}
		}
	}
	rkt = p.PenaltyBoundary * tpen
	return
}

// CalcRateLogLike ...
func (p *PLObj) CalcRateLogLike() (ll float64) {
	ll = 0.0
	for i := 1; i < p.NumNodes; i++ {
		if i == p.CVNode { //skip the cv node
			continue
		}
		x := p.Rates[i] * p.Durations[i]
		c := p.CharDurations[i]
		l := 0.0
		if x > 0.0 {
			l = -(c*math.Log(x) - x - p.LogFactCharDurations[i])
		} else if x == 0 {
			if c > 0.0 {
				l = 1e+15
			} else if c == 0 {
				l = 0.0
			}
		}
		ll += l
	}
	return
}

// CalcRateLogLikeP ...
// NOT FASTER
func (p *PLObj) CalcRateLogLikeP() (ll float64) {
	//cpu :=
	//runtime.GOMAXPROCS(cpu)
	cpu := runtime.GOMAXPROCS(0)
	chunk := (p.NumNodes - 1) / cpu
	ll = 0.0
	lr := make(chan float64)
	for j := 0; j < cpu; j++ {
		//	for i := 1; i < p.NumNodes; i++ {
		go func(start int, r chan<- float64) {
			end := start + chunk
			if end >= p.NumNodes {
				end = p.NumNodes
			}
			l2 := 0.0
			for i := start; i < end; i++ {
				if i == 0 {
					continue
				}
				x := p.Rates[i] * p.Durations[i]
				c := p.CharDurations[i]
				l := 0.0
				if x > 0.0 {
					l = -(c*math.Log(x) - x - p.LogFactCharDurations[i])
				} else if x == 0 {
					if c > 0.0 {
						l = 1e+15
					} else if c == 0 {
						l = 0.0
					}
				}
				l2 += l
			}
			r <- l2
		}(j*chunk, lr)
	}
	for i := 0; i < cpu; i++ {
		ll += <-lr
	}
	return
}

func (p *PLObj) PCalcRateLogLike(jobs <-chan int, results chan<- float64) {
	for i := range jobs {
		//ll = 0.0
		//for i := 1; i < p.NumNodes; i++ {
		//if i == p.CVNode { //skip the cv node
		//	continue
		//}
		x := p.Rates[i] * p.Durations[i]
		c := p.CharDurations[i]
		l := 0.0
		if x > 0.0 {
			l = -(c*math.Log(x) - x - p.LogFactCharDurations[i])
		} else if x == 0 {
			if c > 0.0 {
				l = 1e+15
			} else if c == 0 {
				l = 0.0
			}
		}
		results <- l
	}
	return
}

// CalcRoughnessPenalty ...
func (p *PLObj) CalcRoughnessPenalty() (su float64) {
	su = 0
	tomy := 0.0
	s := 0.0
	ss := 0.0
	for i := 0; i < p.NumNodes; i++ {
		if i == p.CVNode { // skip cvnode
			continue
		}
		if i == 0 { //root
			for _, j := range p.ChildrenVec[i] {
				r := p.Rates[j]
				if p.LogPen {
					r = math.Log(r)
				}
				s += r
				ss += r * r
				tomy++
			}
		} else {
			if p.ParentsNdsInts[i] != 0 {
				sm := p.Rates[i] - p.Rates[p.ParentsNdsInts[i]]
				su += (sm * sm)
			}
		}
	}
	su += (ss - s*s/tomy) / tomy
	return
}

// CalcRoughnessPenaltyMult ...
func (p *PLObj) CalcRoughnessPenaltyMult() (su float64) {
	su = 0
	tomy := 0.0
	s := 0.0
	ss := 0.0
	for i := 0; i < p.NumNodes; i++ {
		if i == 0 { //root
			for _, j := range p.ChildrenVec[i] {
				r := p.Rates[j]
				if p.LogPen {
					r = math.Log(r)
				}
				s += r
				ss += r * r
				tomy++
			}
		} else {
			if p.ParentsNdsInts[i] != 0 && p.RateGroups[p.ParentsNdsInts[i]] == p.RateGroups[i] {
				sm := p.Rates[i] - p.Rates[p.ParentsNdsInts[i]]
				su += (sm * sm)
			} //need to figure out the rate penalty when for the edges that subtend the shift
		}
	}
	su += (ss - s*s/tomy) / tomy
	return
}

// SetDurations for nodes
func (p *PLObj) SetDurations() (success bool) {
	success = true
	for i := 0; i < p.NumNodes; i++ {
		if p.Dates[i] < p.Mins[i] {
			success = false
			return
		}
		if p.Dates[i] > p.Maxs[i] {
			success = false
			return
		}
		if i != 0 {
			p.Durations[i] = p.Dates[p.ParentsNdsInts[i]] - p.Dates[i]
			if p.Durations[i] < 0 {
				success = false
				return
			}
		}
		if p.Durations[i] == 0 || math.IsNaN(p.Durations[i]) {
			p.Durations[i] = math.SmallestNonzeroFloat64
		}
	}
	return
}

//PrintNewickDurations uses the NewickFloatBL
func (p *PLObj) PrintNewickDurations(t Tree) string {
	for _, n := range t.Pre {
		n.FData["duration"] = p.Durations[n.Num]
	}
	return t.Rt.NewickFloatBL("duration")
}

//without CV
func (p *PLObj) calcLFGradient(params []float64, g *[]float64) {
	//rate derivative
	count := 0
	rate := params[0]
	sumbl := 0.
	sumtimes := 0.
	for i := 0; i < p.NumNodes; i++ {
		sumbl += p.CharDurations[i]
		sumtimes += p.Durations[i]
	}
	(*g)[count] = -(sumbl/rate - sumtimes)

	if math.IsNaN((*g)[count]) {
		for i := 0; i < len(params); i++ {
			(*g)[i] = 1000000
		}
		return
	}
	count++

	//time derivative
	for _, i := range p.FreeNodes {
		//	for i := 0; i < p.NumNodes; i++ {
		//		if _, ok := p.FreeNodesM[i]; ok { //use to be free->(at) == 1
		(*g)[count] = 0
		if i != 0 { //not the root
			if p.CharDurations[i] == 0.0 {
				(*g)[count] = rate
			} else {
				(*g)[count] = -p.CharDurations[i]/p.Durations[i] + rate
			}
		}
		for j := 0; j < len(p.ChildrenVec[i]); j++ { //WAY FASTER THAN ABOVE
			curch := p.ChildrenVec[i][j]
			if p.CharDurations[curch] == 0.0 {
				(*g)[count] -= rate
			} else {
				(*g)[count] += p.CharDurations[curch]/p.Durations[curch] - rate
			}
		}

		(*g)[count] *= -1
		count++
		//		}
	}
}

func (p *PLObj) calcPLGradient(params []float64, g *[]float64) {
	icount := p.NumNodes - 1 //seems to be numnodes not -1, icount := p.NumNodes - 1
	if p.CVNode > 0 {
		icount--
	}

	//int subcv = calc_isum(p.cvnodes);
	//icount -= subcv;
	//time derivative
	for _, i := range p.FreeNodes {
		//	for i := 0; i < p.NumNodes; i++ {
		//		if _, ok := p.FreeNodesM[i]; ok { //if (free->at(curponi) == 1){// greater than 0 is free, negative is not free
		(*g)[icount] = 0.0
		if i != 0 { //not the root
			if p.CharDurations[i] == 0.0 {
				(*g)[icount] = p.Rates[i]
			} else {
				//g->at(icount) = -nd->char_duration/nd->duration+nd->rate;
				(*g)[icount] = -p.CharDurations[i]/p.Durations[i] + p.Rates[i]
			}
		}
		for j := 0; j < len(p.ChildrenVec[i]); j++ {
			curch := p.ChildrenVec[i][j]
			if p.CVNode != curch { //if cvnodes[curch] == 0 {
				if p.CharDurations[curch] == 0.0 {
					(*g)[icount] -= p.Rates[curch]
				} else {
					(*g)[icount] += p.CharDurations[curch]/p.Durations[curch] - p.Rates[curch]
				}
			}
		}
		(*g)[icount] *= -1
		icount++
		//		}
	}
	//fmt.Println(icount, params, g)
	//os.Exit(0)
	//rate derivative//this one is additive
	icount = 0
	for i := 1; i < p.NumNodes; i++ {
		if p.CVNode != i { //if(cvnodes[i] == 0){
			tomy := 0.
			meanr := 0.0
			tg := 0.0
			var lograte float64
			//if i != 0 { //not the root
			tg = p.CharDurations[i]/p.Rates[i] - p.Durations[i]
			if p.LogPen == true {
				lograte = math.Log(p.Rates[i])
			}
			if p.ParentsNdsInts[i] == 0 { //parent is the root
				curp := p.ParentsNdsInts[i]
				for j := 0; j < len(p.ChildrenVec[curp]); j++ {
					curch := p.ChildrenVec[curp][j]
					if p.CVNode != curch { //if(cvnodes[curch] == 0){
						tomy++ //++tomy;
						if p.LogPen {
							meanr += math.Log(p.Rates[curch])
						} else {
							meanr += p.Rates[curch]
						}
					}
				}
				meanr /= tomy
				if p.LogPen {
					tg += -(2 * p.Smoothing / p.Rates[i]) * (lograte - meanr) / tomy
				} else {
					tg += -(2 * p.Smoothing) * (p.Rates[i] - meanr) / tomy
				}
				//		    cout << "f "<<  tg << endl;
				for j := 0; j < len(p.ChildrenVec[i]); j++ {
					curch := p.ChildrenVec[i][j]
					if p.CVNode != curch { //if (cvnodes[curch] == 0){
						if p.LogPen {
							tg += 2 * p.Smoothing * (math.Log(p.Rates[curch]) - lograte) / p.Rates[i]
						} else {
							tg += 2 * p.Smoothing * (p.Rates[curch] - p.Rates[i])
						}
					}
				}
				//		    cout << tg << endl;
			} else {
				curp := p.ParentsNdsInts[i]
				if p.LogPen {
					tg += (-2 * p.Smoothing) * (lograte - math.Log(p.Rates[curp])) / p.Rates[i]
				} else {
					tg += (-2 * p.Smoothing) * (p.Rates[i] - p.Rates[curp])
				}
				for j := 0; j < len(p.ChildrenVec[i]); j++ {
					curch := p.ChildrenVec[i][j]
					if p.CVNode != curch { //if (cvnodes[curch] == 0){
						if p.LogPen {
							tg += 2 * p.Smoothing * (math.Log(p.Rates[curch]) - lograte) / p.Rates[i]
						} else {
							tg += 2 * p.Smoothing * (p.Rates[curch] - p.Rates[i])
						}
					}
				}
				//		    cout << "notroot " << -tg << endl;
			}
			(*g)[icount] = -tg //-1;
			icount++
			//}
		}
	}
}
