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
			log.Printf("read %d bytes: %v", ln, err)
			break
		}
	}
	fmt.Fprint(os.Stderr, "\n")

	for _, t := range trees {
		minmap := make(map[*gophy.Node]float64, 1)
		maxmap := make(map[*gophy.Node]float64, 1)
		minmap[t.Rt] = 10.0
		maxmap[t.Rt] = 10.0
		p := gophy.PLObj{}
		p.SetValues(t, 100., minmap, maxmap)
		p.Smoothing = 1
	}
}
