package main

import (
	"flag"
	"fmt"
	"math/rand"
	"time"

	"github.com/FePhyFoFum/gophy"
)

func main() {
	rand.Seed(time.Now().UTC().UnixNano())
	treefn := flag.String("t", "", "input tree")
	traitfn := flag.String("m", "", "continuous traits")
	flag.Parse()
	nwk := gophy.ReadLine(*treefn)[0]
	//rt := cophymaru.ReadTree(nwk)
	rt := gophy.ReadNewickString(nwk)
	t := gophy.NewTree()
	t.Instantiate(rt)
	for _, nd := range t.Pre {
		if nd != t.Rt {
			nd.BMLen = rand.Float64()
		}
	}
	gophy.MapContinuous(t, *traitfn)
	//gophy.IterateBMLengths(t, 100)
	gophy.InitMissingValues(t.Pre)
	var sites []int
	for i := range t.Rt.ContData {
		sites = append(sites, i)
	}
	start := time.Now()
	gophy.BMOptimBLEM(t, 100)
	fmt.Println(gophy.SubUnrootedLogLikeParallel(t.Rt, sites, 3))
	elapsed := time.Since(start)
	fmt.Println(elapsed)
	fmt.Println(t.Rt.BMPhylogram())
}
