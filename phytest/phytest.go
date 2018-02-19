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

	"gonum.org/v1/gonum/floats"
)

// PCalcLogLike this will calculate log like in parallel
func PCalcLogLike(t gophy.Tree, x gophy.DNAModel, nsites int, wks int) (fl float64) {
	fl = 0.0
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	for i := 0; i < wks; i++ {
		go CalcLogLikeWork(t, x, jobs, results)
	}
	for i := 0; i < nsites; i++ {
		jobs <- i
	}
	close(jobs)
	for i := 0; i < nsites; i++ {
		fl += <-results
	}
	return
}

// CalcLogLikeWork this is intended for a worker that will be executing this per site
func CalcLogLikeWork(t gophy.Tree, x gophy.DNAModel, jobs <-chan int, results chan<- float64) {
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

// CalcLogLikeNode calculates likelihood for node
func CalcLogLikeNode(nd *gophy.Node, model gophy.DNAModel, site int) {
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
	x := gophy.DNAModel{}
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
	var t gophy.Tree
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
		//populate Pmatrices
		x.SetP(n.Len)
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
	//fmt.Println(math.Log(l))
}
