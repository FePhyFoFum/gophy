package gophy

import (
	"fmt"
	"math"
	"math/rand"
	"time"

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
	NumNodes             int
	LogPen               bool
	PenaltyBoundary      float64
	Smoothing            float64
}

// OptimizePL PL optimization
func (p *PLObj) OptimizePL(params []float64) {
	count := 0
	start := time.Now()
	fcn := func(pa []float64) float64 {
		for _, i := range pa {
			if i < 0 {
				return 1000000000000
			}
		}
		pl := p.CalcPL(pa)
		if pl == -1 {
			return 100000000000
		}
		if count%10 == 0 {
			curt := time.Now()
			fmt.Println(count, pl, curt.Sub(start))
			start = time.Now()
		}
		count++
		return pl
	}
	/*grad := func(grad, x []float64) {
		fd.Gradient(grad, fcn, x, nil)
	}*/
	settings := optimize.Settings{}
	//settings.UseInitialData = false
	//settings.FunctionThreshold = 1e-12
	//settings.GradientThreshold = 0.01
	settings.Concurrent = 0
	//settings.Recorder = nil
	FC := optimize.FunctionConverge{}
	FC.Absolute = 1
	FC.Relative = 1
	FC.Iterations = 100
	//settings.FunctionConverge = &FC
	ps := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	p0 = params
	res, err := optimize.Minimize(ps, p0, &settings, nil)
	if err != nil {
		fmt.Println(err)
	}
	fmt.Println(res.F)
	fmt.Println(res)
}

// CalcPL ...
func (p *PLObj) CalcPL(params []float64) float64 {
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
	for i := range p.Dates {
		p.Dates[i] = params[pcount]
		pcount++
		fcount++
	}
	ret := p.SetDurations()
	if ret == false {
		return -1.
	}
	ll := p.CalcRateLike()
	tp := p.CalcPenalty()
	rp := p.CalcRoughnessPenalty()
	pl := (ll + tp) + p.Smoothing*rp
	return pl
}

// SetValues ...
func (p *PLObj) SetValues(t Tree, numsites float64, minmap map[*Node]float64,
	maxmap map[*Node]float64) {
	p.NumNodes = 0
	for i, n := range t.Pre {
		n.Num = i
		p.NumNodes++
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
		p.CharDurations[i] = math.Round(n.Len * numsites)
		p.LogFactCharDurations[i] = LogFact(p.CharDurations[i])
	}
	// set min maxs
	p.Maxs = make([]float64, p.NumNodes)
	p.Mins = make([]float64, p.NumNodes)
	p.PenMins = make([]float64, p.NumNodes)
	p.PenMaxs = make([]float64, p.NumNodes)
	for _, n := range t.Post {
		if len(n.Chs) > 0 { //internal
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
	}
	//end set min max
	fmt.Println(p.CharDurations)
	fmt.Println(p.LogFactCharDurations)
	fmt.Println(p.ParentsNdsInts)
	fmt.Println(p.ChildrenVec)
	fmt.Println(p.Maxs)
	fmt.Println(p.Mins)
	//setup dates
	rtDt := 10.0
	p.Dates[t.Rt.Num] = 10.0
	for _, n := range t.Rt.Chs {
		p.feasibleTimes(n, rtDt)
	}
	fmt.Println(p.Dates)
	fmt.Println(p.SetDurations())
	params := make([]float64, p.NumNodes*2)
	c := 0
	for i := range p.Rates {
		params[i] = 1.0
		c++
	}
	for i := range p.Dates {
		params[c] = p.Dates[i]
		c++
	}
	fmt.Println(p.CalcPL(params))
	p.OptimizePL(params)
}

func (p *PLObj) feasibleTimes(nd *Node, timeAnc float64) {
	if len(nd.Chs) == 0 {
		p.Dates[nd.Num] = 0
		return
	}
	fmt.Println(nd, nd.Num)
	if p.Maxs[nd.Num] < timeAnc {
		timeAnc = p.Maxs[nd.Num]
	}
	p.Dates[nd.Num] = timeAnc - (timeAnc-p.Mins[nd.Num])*(0.02*rand.Float64()*0.96)/math.Log(2.+3.)
	fmt.Println(p.Dates[nd.Num], timeAnc)
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

// CalcRateLike ...
func (p *PLObj) CalcRateLike() (ll float64) {
	ll = 0.0
	for i := 1; i < p.NumNodes; i++ {
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

// CalcRoughnessPenalty ...
func (p *PLObj) CalcRoughnessPenalty() (su float64) {
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
			if p.ParentsNdsInts[i] != 0 {
				sm := p.Rates[i] - p.Rates[p.ParentsNdsInts[i]]
				su += (sm * sm)
			}
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
