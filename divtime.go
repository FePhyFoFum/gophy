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
	NumNodes             int
	LogPen               bool
	PenaltyBoundary      float64
	Smoothing            float64
	RateGroups           map[int]int //[nodenum]rategroup
	CVNode               int
	Tree                 Tree
}

// OptimizeRD optimization
func (p *PLObj) OptimizeRD(params []float64, likeFunc func([]float64, bool) float64) *optimize.Result {
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
			//fmt.Println(count, pl, curt.Sub(start))
			//start = time.Now()
		}
		count++
		return pl
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
	ps := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	p0 = params
	res, err := optimize.Minimize(ps, p0, &settings, nil)
	if err != nil {
		fmt.Println(err)
	}
	return res
	//fmt.Println(res.F)
	//fmt.Println(res)
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

	/*jobs := make(chan int, p.NumNodes-1)
	results := make(chan float64, p.NumNodes-1)
	for w := 0; w < 2; w++ {
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
	}*/
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
	ll := p.CalcRateLogLike()
	//fmt.Println(ll)
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
	runtime.GOMAXPROCS(2)

	p.Tree = t
	p.NumNodes = 0
	p.FreeNodes = make([]int, 0)
	//assign node numbers
	for i, n := range t.Pre {
		n.Num = i
		p.NumNodes++
		if n != t.Rt && len(n.Chs) > 0 {
			p.FreeNodes = append(p.FreeNodes, n.Num)
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
	//end set min max

	//setup dates
	rtDt := minmap[t.Rt]
	p.Dates[t.Rt.Num] = rtDt
	setMaxNodesToPresent(t)
	for _, n := range t.Rt.Chs {
		//p.feasibleTimes(n, rtDt)
		p.feasibleTimesSplit(n, rtDt) //just splitting them, can go back to the other
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
	x := p.OptimizeRD(params, p.CalcLF)
	if verbose {
		fmt.Fprintln(os.Stderr, "end LF:", x.F)
	}
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
	x := p.OptimizeRD(params, p.CalcMultLF)
	if verbose {
		fmt.Fprintln(os.Stderr, "end MLF:", x.F)
	}
	return x
}

//RunPL penalized likelihood run with a starting float, probably x.X[0] from LF
func (p *PLObj) RunPL(startRate float64, verbose bool) *optimize.Result {
	params := make([]float64, p.NumNodes+len(p.FreeNodes)) // pl
	c := 0
	for i := range p.Rates { //every edge has a rate
		params[i] = startRate
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
	x := p.OptimizeRD(params, p.CalcPL)
	if verbose {
		fmt.Fprintln(os.Stderr, "end PL:", x.F)
	}
	return x
}

//RunMPL penalized likelihood run with a starting float from x.X[0] from LF
// and mrcagroup -- assuming that p.RateGroups has already been setup
func (p *PLObj) RunMPL(mrcagroups []*Node, t Tree, verbose bool) *optimize.Result {
	params := make([]float64, p.NumNodes+len(p.FreeNodes))
	c := 0
	//fmt.Fprintln(os.Stderr, p.Rates)
	for i, j := range p.Rates { //every edge has a rate
		params[i] = j
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
	x := p.OptimizeRD(params, p.CalcMultPL)
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
		p.OptimizeRD(params, likeFunc)
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
	ll = 0.0
	lr := make(chan float64)
	for i := 1; i < p.NumNodes; i++ {
		go func(i int, lr chan float64) {
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
			lr <- l
		}(i, lr)
	}
	for i := 1; i < p.NumNodes; i++ {
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

//
