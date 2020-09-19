package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime/pprof"

	"github.com/FePhyFoFum/gophy"
)

func main() {
	tfn := flag.String("t", "", "tree filename")
	wks := flag.Int("w", 4, "number of threads")
	v := flag.Bool("v", false, "verbose")
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")
	flag.Parse()
	fmt.Fprintln(os.Stderr, *tfn, *wks, *v)
	if len(os.Args) < 2 {
		flag.PrintDefaults()
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
	//reading trees
	f, err := os.Open(*tfn)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	defer f.Close()
	scanner := bufio.NewReader(f)
	fmt.Fprint(os.Stderr, "reading trees\n")
	trees := make([]gophy.Tree, 0)
	for {
		ln, err := scanner.ReadString('\n')
		if len(ln) > 0 {
			rt := gophy.ReadNewickString(ln)
			var t gophy.Tree
			t.Instantiate(rt)
			trees = append(trees, t)
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Printf("error: %v", err)
			break
		}
	}
	fmt.Fprint(os.Stderr, "\n")

	for _, t := range trees {
		minmap := make(map[*gophy.Node]float64, 1)
		maxmap := make(map[*gophy.Node]float64, 1)
		minmap[t.Rt] = 100.0
		maxmap[t.Rt] = 100.0
		p := gophy.PLObj{}
		p.Smoothing = 1.0
		numsites := 10000.
		p.SetValues(t, numsites, minmap, maxmap)
		x := p.RunLF(numsites / 20.)
		fmt.Println(p.PrintNewickDurations(t) + ";")
		nds := make([]*gophy.Node, 2)
		nds[0], _ = t.GetTipByName("taxon_1")
		nds[1], _ = t.GetTipByName("taxon_12")
		nd := gophy.GetMrca(nds, t.Rt)
		mrcagroups := make([]*gophy.Node, 1)
		mrcagroups[0] = nd
		x = p.RunMLF(x.X[0], mrcagroups, t)
		fmt.Println(p.PrintNewickDurations(t) + ";")
		//x = p.RunPL(x.X[0])
		//fmt.Println(p.PrintNewickDurations(t) + ";")
	}
}
