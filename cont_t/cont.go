package main

import (
	"flag"
	"fmt"
	"math/rand"
	"time"

	"github.com/FePhyFoFum/gophy"
	"github.com/FePhyFoFum/gophy/clustering"
)

func main() {
	rand.Seed(time.Now().UTC().UnixNano())
	treefn := flag.String("t", "", "input tree")
	traitfn := flag.String("m", "", "continuous traits")
	algfn := flag.String("f", "", "0 optimize branch lengths and output loglikelihood\n1 output loglikelihoods for each trait")
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
	clustering.InitMissingValues(t.Pre)
	var sites []int
	for i := range t.Rt.ContData {
		sites = append(sites, i)
	}
	start := time.Now()
	clustering.BMOptimBLEM(t, 100)
	if *algfn == "0" {
		fmt.Println(clustering.SubUnrootedLogLikeParallel(t.Rt, sites, 3))
	} else if *algfn == "1" {
		for _, site := range sites {
			sll := clustering.SingleSiteLL(t.Rt, site)
			fmt.Println(site, sll)
		}
	}
	elapsed := time.Since(start)
	fmt.Println(elapsed)
	fmt.Println(t.Rt.BMPhylogram())
}
