package main

import (
	"bufio"
	"flag"
	"fmt"
	"gophy"
	"log"
	"math"
	"os"
	"runtime/pprof"
	"time"

	"golang.org/x/exp/rand"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/optimize"
)

// OptimizeBL ...
func OptimizeBL(nd *gophy.Node, t *gophy.Tree, x *gophy.DNAModel, nsites int, wks int) {
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
	settings := optimize.DefaultSettings()
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
	res, err := optimize.Local(p, p0, settings, nil)
	if err != nil {
		//fmt.Println(err)
	}
	nd.Len = res.X[0]
}

func slidingWindow(x, w float64) float64 {
	return (rand.Float64() * ((x + w/2) - (x - w/2))) + (x - w/2.)
}

// MCMC simple MCMC for branch lengths
func MCMC(t *gophy.Tree, x *gophy.DNAModel, nsites int, wks int, outfilename string) {
	f, err := os.Create(outfilename)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	wr := bufio.NewWriter(f)
	w := 0.5
	printf := 100
	iter := 10000
	curlike := PCalcLogLike(t, x, nsites, wks)
	newlike := 0.0
	start := time.Now()
	//fullcalc := false
	for i := 0; i < iter; i++ {
		nd := rand.Intn(len(t.Post))
		for t.Post[nd] == t.Rt {
			nd = rand.Intn(len(t.Post))
		}
		odv := t.Post[nd].Len
		newv := slidingWindow(odv, w)
		t.Post[nd].Len = newv
		t.Post[nd].Marked = true
		if len(t.Post[nd].Chs) == 0 {
			t.Post[nd].Par.Marked = true
		}
		//if fullcalc {
		//newlike = PCalcLogLike(t, x, nsites, wks)
		//} else {
		//	newlike = PCalcLogLikeBack(t, t.Post[nd], x, nsites, wks)
		//}
		newlike = PCalcLogLikeMarked(t, x, nsites, wks)
		for _, n := range t.Post {
			n.Marked = false
		}
		r := math.Log(rand.Float64())
		if r < newlike-curlike {
			//if newlike-curlike > 0 {
			curlike = newlike
			//fullcalc = false
		} else {
			t.Post[nd].Len = odv
			t.Post[nd].Marked = true
			if len(t.Post[nd].Chs) == 0 {
				t.Post[nd].Par.Marked = true
			}
			//fullcalc = true
		}
		if i%printf == 0 {
			curt := time.Now()
			fmt.Fprintln(os.Stderr, i, curlike, "(last", printf, ":", curt.Sub(start), ")")
			fmt.Fprint(wr, t.Rt.Newick(true)+";\n")
			start = time.Now()
		}
	}
	end := time.Now()
	fmt.Println(end.Sub(start))
	err = wr.Flush()
	if err != nil {
		log.Fatal(err)
	}
}

// PCalcLogLike this will calculate log like in parallel
func PCalcLogLike(t *gophy.Tree, x *gophy.DNAModel, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += CalcLogLikeOneSite(t, x, 0)
	for i := 0; i < wks; i++ {
		go CalcLogLikeWork(t, x, jobs, results)
	}
	for i := 1; i < nsites; i++ {
		jobs <- i
	}
	close(jobs)
	for i := 1; i < nsites; i++ {
		fl += <-results
	}
	return
}

// PCalcLogLikeBack a bit of a shortcut. Could do better, but walks back from the n node to the root
func PCalcLogLikeBack(t *gophy.Tree, n *gophy.Node, x *gophy.DNAModel, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	x.EmptyPDict()
	fl += CalcLogLikeOneSiteBack(t, n, x, 0)
	for i := 0; i < wks; i++ {
		go CalcLogLikeWorkBack(t, n, x, jobs, results)
	}
	for i := 1; i < nsites; i++ {
		jobs <- i
	}
	close(jobs)
	for i := 1; i < nsites; i++ {
		fl += <-results
	}
	return
}

func PCalcLogLikeMarked(t *gophy.Tree, x *gophy.DNAModel, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	// populate the P matrix dictionary without problems of race conditions
	// just the first site
	//x.EmptyPDict()
	fl += CalcLogLikeOneSiteMarked(t, x, 0)
	for i := 0; i < wks; i++ {
		go CalcLogLikeWorkMarked(t, x, jobs, results)
	}
	for i := 1; i < nsites; i++ {
		jobs <- i
	}
	close(jobs)
	for i := 1; i < nsites; i++ {
		fl += <-results
	}
	return
}

// CalcLogLikeOneSite just calculate the likelihood of one site
// probably used to populate the PDict in the DNA Model so that we can reuse the calculations
func CalcLogLikeOneSite(t *gophy.Tree, x *gophy.DNAModel, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			CalcLogLikeNode(n, x, site)
		}
		if t.Rt == n {
			for i := 0; i < 4; i++ {
				t.Rt.Data[site][i] += math.Log(0.25)
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLogLikeOneSiteBack like the one above but from nb to the root only
func CalcLogLikeOneSiteBack(t *gophy.Tree, nb *gophy.Node, x *gophy.DNAModel, site int) float64 {
	sl := 0.0
	going := true
	cur := nb
	for going {
		if len(cur.Chs) > 0 {
			CalcLogLikeNode(cur, x, site)
		}
		if cur == t.Rt {
			for i := 0; i < 4; i++ {
				t.Rt.Data[site][i] += math.Log(0.25)
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
			going = false
			break
		}
		cur = cur.Par
	}
	return sl
}

func CalcLogLikeOneSiteMarked(t *gophy.Tree, x *gophy.DNAModel, site int) float64 {
	sl := 0.0
	for _, n := range t.Post {
		if len(n.Chs) > 0 {
			if n.Marked == true {
				CalcLogLikeNode(n, x, site)
				if n != t.Rt {
					n.Par.Marked = true
				}
			}
		}
		if t.Rt == n && n.Marked == true {
			for i := 0; i < 4; i++ {
				t.Rt.Data[site][i] += math.Log(0.25)
			}
			sl = floats.LogSumExp(t.Rt.Data[site])
		} else {
			sl = floats.LogSumExp(t.Rt.Data[site])
		}
	}
	return sl
}

// CalcLogLikeWork this is intended for a worker that will be executing this per site
func CalcLogLikeWork(t *gophy.Tree, x *gophy.DNAModel, jobs <-chan int, results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				CalcLogLikeNode(n, x, j)
			}
			if t.Rt == n {
				for i := 0; i < 4; i++ {
					t.Rt.Data[j][i] += math.Log(0.25)
				}
				sl = floats.LogSumExp(t.Rt.Data[j])
			}
		}
		results <- sl
	}
}

// CalcLogLikeWorkBack this is intended for a worker that will be executing this per site
func CalcLogLikeWorkBack(t *gophy.Tree, nb *gophy.Node, x *gophy.DNAModel, jobs <-chan int, results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		going := true
		cur := nb
		for going {
			if len(cur.Chs) > 0 {
				CalcLogLikeNode(cur, x, j)
			}
			if cur == t.Rt {
				for i := 0; i < 4; i++ {
					t.Rt.Data[j][i] += math.Log(0.25)
				}
				sl = floats.LogSumExp(t.Rt.Data[j])
				going = false
				break
			}
			cur = cur.Par
		}
		results <- sl
	}
}

// CalcLogLikeWorkMarked this is intended to calculate only on the marked nodes back to teh root
func CalcLogLikeWorkMarked(t *gophy.Tree, x *gophy.DNAModel, jobs <-chan int, results chan<- float64) {
	for j := range jobs {
		sl := 0.0
		for _, n := range t.Post {
			if len(n.Chs) > 0 {
				if n.Marked == true {
					CalcLogLikeNode(n, x, j)
				}
			}
			if t.Rt == n && n.Marked == true {
				for i := 0; i < 4; i++ {
					t.Rt.Data[j][i] += math.Log(0.25)
				}
				sl = floats.LogSumExp(t.Rt.Data[j])
			} else {
				sl = floats.LogSumExp(t.Rt.Data[j])
			}
		}
		results <- sl
	}
}

// CalcLogLikeNode calculates likelihood for node
func CalcLogLikeNode(nd *gophy.Node, model *gophy.DNAModel, site int) {
	for i := 0; i < 4; i++ {
		nd.Data[site][i] = 0.
	}
	x1 := 0.0
	x2 := []float64{0.0, 0.0, 0.0, 0.0}
	for _, c := range nd.Chs {
		P := model.GetPMap(c.Len)
		if len(c.Chs) == 0 {
			for i := 0; i < 4; i++ {
				x1 = 0.0
				for j := 0; j < 4; j++ {
					x1 += P.At(i, j) * c.Data[site][j]
				}
				nd.Data[site][i] += math.Log(x1)
			}
		} else {
			for i := 0; i < 4; i++ {
				for j := 0; j < 4; j++ {
					x2[j] = math.Log(P.At(i, j)) + c.Data[site][j]
				}
				nd.Data[site][i] += floats.LogSumExp(x2)
			}
		}
	}
}

func main() {
	rand.Seed(uint64(time.Now().UTC().UnixNano()))
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")
	flag.Parse()
	if len(*tfn) == 0 {
		fmt.Fprintln(os.Stderr, "need a filename")
		os.Exit(1)
	}
	if len(*afn) == 0 {
		fmt.Fprintln(os.Stderr, "need a filename")
		os.Exit(1)
	}
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}
	x := gophy.NewDNAModel()
	x.SetupQ()
	x.SetNucMap()

	//read a seq file
	nsites := 0
	seqs := map[string]string{}
	for _, i := range gophy.ReadSeqsFromFile(*afn) {
		seqs[i.NM] = i.SQ
		nsites = len(i.SQ)
	}
	//read a tree file
	f, err := os.Open(*tfn)
	if err != nil {
		fmt.Println(err)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	var rt *gophy.Node
	t := gophy.NewTree()
	for scanner.Scan() {
		ln := scanner.Text()
		if len(ln) < 2 {
			continue
		}
		rt = gophy.ReadNewickString(ln)
		t.Instantiate(rt)
	}
	//each site has a vector of 4 for DNA conditionals
	for _, n := range t.Post {
		n.Data = make([][]float64, nsites)
		for i := 0; i < nsites; i++ {
			n.Data[i] = []float64{0.0, 0.0, 0.0, 0.0}
		}
		if len(n.Chs) == 0 {
			for i := 0; i < nsites; i++ {
				n.Data[i][x.CharMap[string(seqs[n.Nam][i])]] = 1.0
			}
		}
	}

	// calc likelihood
	start := time.Now()
	w := 10
	if nsites < w {
		w = nsites
	}
	l := PCalcLogLike(t, x, nsites, 10)
	end := time.Now()
	fmt.Println("lnL:", l)
	fmt.Println(end.Sub(start))
	//fmt.Println(t.Rt.Newick(true))
	start = time.Now()
	/*
		for _, n := range t.Post {
			OptimizeBL(n, t, x, nsites, 10)
		}
		for _, n := range t.Post {
			OptimizeBL(n, t, x, nsites, 10)
		}
		for _, n := range t.Post {
			OptimizeBL(n, t, x, nsites, 10)
		}
		for _, n := range t.Post {
			OptimizeBL(n, t, x, nsites, 10)
		}
		for _, n := range t.Post {
			OptimizeBL(n, t, x, nsites, 10)
		}*/
	MCMC(t, x, nsites, 10, "temp.mcmc.tre")
	end = time.Now()
	fmt.Println(end.Sub(start))
}
