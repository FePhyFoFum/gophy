package main

import (
	"flag"
	"fmt"
	"gophy"
	"math/rand"
	"time"
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
	fmt.Println(t)
	gophy.MapContinuous(t, *traitfn)
	//gophy.IterateBMLengths(t, 100)
	gophy.InitMissingValues(t.Pre)
	start := time.Now()
	for i := 0; i < 4; i++ {
		gophy.CalcExpectedTraits(t.Rt)
		gophy.BMCalcLensBackFront(t)
	}
	elapsed := time.Since(start)
	fmt.Println(elapsed)
	fmt.Println(t.Rt.BMPhylogram())
}
