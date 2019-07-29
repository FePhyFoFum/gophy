package gophy

import (
	"fmt"
	"math"
	"os"
)

func rootedDispSiteLL(n *Node, nlikes *float64, startFresh bool, site int) {
	for _, chld := range n.Chs {
		rootedDispSiteLL(chld, nlikes, startFresh, site)
	}
	nchld := len(n.Chs)
	if n.Marked == true {
		if startFresh == false {
			if nchld != 0 {
				*nlikes += n.LL[site]
			}
		}
	}
	if n.Marked == false || startFresh == true {
		n.ConPruneLen[site] = n.BMLen
		log2pi := 1.8378770664093453
		if nchld != 0 {
			if nchld != 2 {
				fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
				os.Exit(0)
			}
			c0 := n.Chs[0]
			c1 := n.Chs[1]
			curlike := float64(0.0)
			var tempChar float64
			curVar := (c0.ConPruneLen[site]) + (c1.ConPruneLen[site])
			contrast := c0.ContData[site] - c1.ContData[site]
			curlike += ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2) / (curVar))))
			tempChar = (((c0.ConPruneLen[site]) * c1.ContData[site]) + ((c1.ConPruneLen[site]) * c0.ContData[site])) / (curVar)
			n.ContData[site] = tempChar
			*nlikes += curlike
			tempBranchLength := n.ConPruneLen[site] + (((c0.ConPruneLen[site]) * (c1.ConPruneLen[site])) / ((c0.ConPruneLen[site]) + (c1.ConPruneLen[site]))) // need to calculate the prune length by adding the averaged lengths of the daughter nodes to the length
			n.ConPruneLen[site] = tempBranchLength                                                                                                            // need to calculate the "prune length" by adding the length to the uncertainty
			n.LL[site] = curlike
			//n.Marked = true
		}
	}
}

func rootedMissingSiteLL(n *Node, nlikes *float64, startFresh bool, site int) {
	for _, chld := range n.Chs {
		rootedDispSiteLL(chld, nlikes, startFresh, site)
	}
	nchld := len(n.Chs)
	if n.Marked == true {
		if startFresh == false {
			if nchld != 0 {
				*nlikes += n.LL[site]
			}
		}
	}
	if n.Marked == false || startFresh == true {
		n.ConPruneLen[site] = n.BMLen
		log2pi := 1.8378770664093453
		if nchld != 0 {
			if nchld != 2 {
				fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
				os.Exit(0)
			}
			c0 := n.Chs[0]
			c1 := n.Chs[1]
			if c0.Mis[site] == false && c1.Mis[site] == false {
				curlike := float64(0.0)
				var tempChar float64
				curVar := (c0.ConPruneLen[site]) + (c1.ConPruneLen[site])
				contrast := c0.ContData[site] - c1.ContData[site]
				curlike += ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2) / (curVar))))
				tempChar = (((c0.ConPruneLen[site]) * c1.ContData[site]) + ((c1.ConPruneLen[site]) * c0.ContData[site])) / (curVar)
				n.ContData[site] = tempChar
				*nlikes += curlike
				tempBranchLength := n.ConPruneLen[site] + (((c0.ConPruneLen[site]) * (c1.ConPruneLen[site])) / ((c0.ConPruneLen[site]) + (c1.ConPruneLen[site]))) // need to calculate the prune length by adding the averaged lengths of the daughter nodes to the length
				n.ConPruneLen[site] = tempBranchLength                                                                                                            // need to calculate the "prune length" by adding the length to the uncertainty
				n.LL[site] = curlike
				//n.Marked = true
			} else if c0.Mis[site] == true && c1.Mis[site] == false {
				n.ConPruneLen[site] += (c1.ConPruneLen[site])
				n.ContData[site] = c1.ContData[site]
				n.LL[site] = 0.0
			} else if c1.Mis[site] == true && c0.Mis[site] == false {
				n.ConPruneLen[site] += c0.ConPruneLen[site]
				n.ContData[site] = c0.ContData[site]
				n.LL[site] = 0.0
			}
		}
	}
}

func rootedTreeLike(tree *Node, startFresh bool, jobs <-chan int, results chan<- float64) {
	for site := range jobs {
		tmpll := 0.
		rootedMissingSiteLL(tree, &tmpll, startFresh, site)
		tmpll = tree.LL[site]
		results <- tmpll
	}
}

//PBMLogLikeRt will calculate the BM log like on a rooted tree
func PBMLogLikeRt(tree *Node, startFresh bool, workers int) (sitelikes float64) {
	nsites := len(tree.Chs[0].ContData)
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	for w := 0; w < workers; w++ {
		go rootedTreeLike(tree, startFresh, jobs, results)
	}

	for site := 0; site < nsites; site++ {
		jobs <- site
	}
	close(jobs)

	for site := 0; site < nsites; site++ {
		sitelikes += <-results
	}
	return
}
