// phytest is just for testing different functions in the phylogenetics stuff
//
package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"runtime/pprof"
	"strconv"
	"strings"
	"time"

	"github.com/FePhyFoFum/gophy"

	"golang.org/x/exp/rand"
)

func slidingWindow(x, w float64) float64 {
	return math.Abs((rand.Float64() * ((x + w/2) - (x - w/2))) + (x - w/2.))
}

// MCMC simple MCMC for branch lengths
func MCMC(t *gophy.Tree, x *gophy.DiscreteModel, patternval []float64, wks int, outfilename string) {
	f, err := os.Create(outfilename)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	wr := bufio.NewWriter(f)
	w := 0.1
	printf := 1000
	iter := 100000
	curlike := gophy.PCalcLogLikePatterns(t, x, patternval, wks)
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
		//newlike = gophy.PCalcLogLikeMarked(t, x, nsites, wks)
		//make pattern marked
		newlike = gophy.PCalcLikePatternsMarked(t, x, patternval, wks)
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

func printModelDNA(modelparams []float64, basefreqs []float64) {
	fmt.Fprintln(os.Stderr, "basefreqs -- A:", basefreqs[0], " C:", basefreqs[1], " G:", basefreqs[2], " T:", basefreqs[3])
	fmt.Fprintln(os.Stderr, "modelparams --")
	fmt.Fprintln(os.Stderr, " - ", modelparams[0], modelparams[1], modelparams[2])
	fmt.Fprintln(os.Stderr, modelparams[0], " - ", modelparams[3], modelparams[4])
	fmt.Fprintln(os.Stderr, modelparams[1], modelparams[3], " - ", 1.0)
	fmt.Fprintln(os.Stderr, modelparams[2], modelparams[4], 1.0, "-")
}

// A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V
func printModelProt(d *gophy.DiscreteModel) {
	fmt.Fprintln(os.Stderr, "resfreqs --A:", d.BF[0], "resfreqs --R:", d.BF[1], "resfreqs --N:", d.BF[2], "resfreqs --D:", d.BF[3], "resfreqs --C:", d.BF[4])
	fmt.Fprintln(os.Stderr, "resfreqs --Q:", d.BF[5], "resfreqs --E:", d.BF[6], "resfreqs --G:", d.BF[7], "resfreqs --H:", d.BF[8], "resfreqs --I:", d.BF[9])
	fmt.Fprintln(os.Stderr, "resfreqs --L:", d.BF[10], "resfreqs --K:", d.BF[11], "resfreqs --M:", d.BF[12], "resfreqs --F:", d.BF[13], "resfreqs --P:", d.BF[14])
	fmt.Fprintln(os.Stderr, "resfreqs --S:", d.BF[15], "resfreqs --T:", d.BF[16], "resfreqs --W:", d.BF[17], "resfreqs --Y:", d.BF[18], "resfreqs --V:", d.BF[19])
	scanner := bufio.NewScanner(strings.NewReader(d.Ex))
	fmt.Fprintln(os.Stderr, "Exchangeabilties (S_ij, symmetric portion only)")
	i := 0
	for scanner.Scan() {
		if i < 20 {
			fmt.Fprintln(os.Stderr, scanner.Text())
			i++
		}
	}
}

func main() {
	rand.Seed(uint64(time.Now().UTC().UnixNano()))
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
	st := flag.String("st", "nuc", "sequence type [nuc/aa/mult]")
	mdr := flag.String("mdr", "1.0,1.0,1.0,1.0,1.0", "five params for GTR (if sequence type == nuc), or x params (if seq type == mult)")
	m := flag.String("m", "JTT", "empirical amino acid [JTT/WAG/LG] (if sequence type == aa)")
	mbf := flag.String("mbf", "emp", "model base frequencies [mod(el)/emp(irical)] (if sequence type == aa)")
	wks := flag.Int("w", 4, "number of threads")
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")
	flag.Parse()
	var dt gophy.DataType
	if len(*tfn) == 0 {
		flag.PrintDefaults()
		os.Exit(1)
	}
	if len(*afn) == 0 {
		flag.PrintDefaults()
		os.Exit(1)
	}
	if len(*st) == 0 {
		flag.PrintDefaults()
		os.Exit(1)
	} else {
		if *st == "nuc" {
			dt = gophy.Nucleotide
		} else if *st == "aa" {
			dt = gophy.AminoAcid
		} else if *st == "mult" {
			dt = gophy.MultiState
		} else {
			fmt.Fprintln(os.Stderr, "sequence type string is not a recognised datatype, please use [nuc/aa/mult]")
			os.Exit(1)
		}
	}
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
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

	fmt.Println(t.Rt.Newick(true))

	//read a seq file
	var nsites int
	var x gophy.DiscreteModel
	var patternval []float64

	if dt == gophy.Nucleotide {
		seqs, patternsint, ns, bf := gophy.ReadPatternsSeqsFromFile(*afn, true)
		nsites = ns
		patternval, _ = gophy.PreparePatternVecs(t, patternsint, seqs)
		y := gophy.NewDNAModel()
		y.M.SetBaseFreqs(bf)
		mds := strings.Split(*mdr, ",")
		modelparams := make([]float64, 5)
		if len(mds) != 5 {
			fmt.Fprintln(os.Stderr, "your model contains ", len(mds), " params, not 5")
			os.Exit(1)
		} else {
			for i, j := range mds {
				f, err := strconv.ParseFloat(j, 64)
				if err != nil {
					fmt.Fprintln(os.Stderr, "problem parsing ", j, " as float in model specs")
					os.Exit(1)
				}
				modelparams[i] = f
			}
		}
		y.M.SetRateMatrix(modelparams)
		y.M.SetupQGTR()
		x = y.M
	} else if dt == gophy.AminoAcid {
		seqs, patternsint, ns, bf := gophy.ReadPatternsSeqsFromFile(*afn, false)
		nsites = ns
		patternval, _ = gophy.PreparePatternVecsProt(t, patternsint, seqs)
		y := gophy.NewProteinModel()
		if *m == "JTT" {
			y.SetRateMatrixJTT()
		} else if *m == "WAG" {
			y.SetRateMatrixWAG()
		} else if *m == "LG" {
			y.SetRateMatrixLG()
		} else {
			fmt.Fprintln(os.Stderr, "amino acid model string not recognized, please use [JTT/WAG/LG]")
			os.Exit(1)
		}
		if *mbf == "mod" {
			y.M.SetModelBF()
		} else if *mbf == "emp" {
			y.M.SetBaseFreqs(bf)
		}
		y.M.SetupQGTR()
		x = y.M
	} else {
		seqs, patternsint, ns, bf, numstates := gophy.ReadPatternsMSeqsFromFile(*afn)
		nsites = ns
		patternval, _ = gophy.PreparePatternVecsMS(t, patternsint, seqs, gophy.GetMap(numstates), numstates)
		mds := strings.Split(*mdr, ",")
		modelparams := make([]float64, len(mds))
		if len(mds) == (((numstates*numstates)-numstates)/2)-1 {
			fmt.Fprintln(os.Stderr, "your model contains ", len(mds), " params, and is symmetric")
		} else if len(mds) == ((numstates*numstates)-numstates)-1 {
			fmt.Fprintln(os.Stderr, "your model contains ", len(mds), " params, and is asymmetric")
		} else {
			fmt.Fprintln(os.Stderr, "not enough parameter values for number of states!")
			os.Exit(1)
		}
		for i, j := range mds {
			f, err := strconv.ParseFloat(j, 64)
			if err != nil {
				fmt.Fprintln(os.Stderr, "problem parsing ", j, " as float in model specs")
				os.Exit(1)
			}
			modelparams[i] = f
		}
		y := gophy.NewMultStateModel(numstates)
		y.M.SetScaledRateMatrix(modelparams, true)
		y.M.SetBaseFreqs(bf)
		y.M.EBF = bf
		y.M.SetupQGTR()
		x = y.M
	}

	l := gophy.PCalcLikePatterns(t, &x, patternval, *wks)
	fmt.Println("starting lnL:", l)
	fmt.Println("getting starting branch lengths (parsimony)")
	s := gophy.PCalcSankParsPatternsMultState(t, x.NumStates, patternval, *wks)
	fmt.Println("sank:", s)
	gophy.EstParsBLMultState(t, x.NumStates, patternval, nsites)
	fmt.Println("start:\n" + t.Rt.Newick(true) + ";")
	//os.Exit(0)

	// calc likelihood
	start := time.Now()
	l = gophy.PCalcLikePatterns(t, &x, patternval, *wks)
	fmt.Println("lnL:", l)

	gophy.OptimizeBLNR(t, &x, patternval, 10)
	l = gophy.PCalcLikePatterns(t, &x, patternval, *wks)
	fmt.Println("final:\n" + t.Rt.Newick(true) + ";")
	fmt.Println("lnL:", l)
	end := time.Now()
	fmt.Fprintln(os.Stderr, end.Sub(start))

}
