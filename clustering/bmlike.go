package clustering

import (
	"fmt"
	"math"
	"os"

	"github.com/FePhyFoFum/gophy"
)

//BMPruneRooted will prune BM branch lens and PICs down to a rooted node
//root node should be a real (ie. bifurcating) root
func BMPruneRooted(n *gophy.Node) {
	for _, chld := range n.Chs {
		BMPruneRooted(chld)
	}
	n.PruneLen = n.BMLen
	nchld := len(n.Chs)
	if nchld != 0 { //&& n.Marked == false {
		var tempChar float64
		if nchld != 2 {
			fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
		}
		c0 := n.Chs[0]
		c1 := n.Chs[1]
		bot := ((1.0 / c0.PruneLen) + (1.0 / c1.PruneLen))
		//fmt.Println(n.Nam, n.BMLen, n.PruneLen)
		n.PruneLen += 1.0 / bot
		if n.PruneLen == 0 {
			fmt.Println("BMPruneRooted()", n.Nam, n.BMLen, n.PruneLen)
			os.Exit(0)
		}
		if math.IsNaN(n.PruneLen) {
			fmt.Println(c0.Nam, c0.BMLen, c1.Nam, c1.BMLen)
			os.Exit(0)
		}
		for i := range n.Chs[0].ContData {
			if math.IsNaN(c1.ContData[i]) || math.IsNaN(c0.ContData[i]) {
				fmt.Println("NaN encountered in BMPruneRooted()")
				fmt.Println(c0.Nam, c0.ContData[i], c0.Mis[i], c0.PruneLen, c1.Nam, c1.ContData[i], c1.Mis[i], c1.PruneLen)
				os.Exit(0)
			}
			tempChar = (((1 / c0.PruneLen) * c1.ContData[i]) + ((1 / c1.PruneLen) * c0.ContData[i])) / bot
			n.ContData[i] = tempChar
		}
	}
}

//BMPruneRootedSingle will prune BM branch lens and calculate PIC of a single trait down to a rooted node
//root node should be a real (ie. bifurcating) root
func BMPruneRootedSingle(n *gophy.Node, i int) {
	for _, chld := range n.Chs {
		BMPruneRootedSingle(chld, i)
	}
	n.PruneLen = n.BMLen
	nchld := len(n.Chs)
	if nchld != 0 { //&& n.Marked == false {
		if nchld != 2 {
			fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
		}
		c0 := n.Chs[0]
		c1 := n.Chs[1]
		bot := ((1.0 / c0.PruneLen) + (1.0 / c1.PruneLen))
		n.PruneLen += 1.0 / bot
		tempCharacter := (((1 / c0.PruneLen) * c1.ContData[i]) + ((1 / c1.PruneLen) * c0.ContData[i])) / bot
		n.ContData[i] = tempCharacter
	}
}

//AssertUnrootedTree is a quick check to make sure the tree passed is unrooted
func AssertUnrootedTree(tree *gophy.Node) {
	if len(tree.Chs) != 3 {
		fmt.Print("BRAncH LenGTHS MUST BE ITERATED ON AN UNROOTED TREE. THIS TREE IS ROOTED.")
		os.Exit(0)
	}
}

//CalcUnrootedLogLike will calculate the log-likelihood of an unrooted tree, while assuming that no sites have missing data.
func CalcUnrootedLogLike(tree *gophy.Node, startFresh bool) (chll float64) {
	chll = 0.0
	for _, ch := range tree.Chs {
		curlike := 0.0
		CalcRootedLogLike(ch, &curlike, startFresh)
		chll += curlike
	}
	sitelikes := 0.0
	var tmpll float64
	var contrast, curVar float64
	for i := range tree.Chs[0].ContData {
		tmpll = 0.
		if tree.Chs[0].Mis[i] == false && tree.Chs[1].Mis[i] == false && tree.Chs[2].Mis[i] == false { //do the standard calculation when no subtrees have missing traits
			contrast = tree.Chs[0].ContData[i] - tree.Chs[1].ContData[i]
			curVar = tree.Chs[0].PruneLen + tree.Chs[1].PruneLen
			tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
			tmpPruneLen := ((tree.Chs[0].PruneLen * tree.Chs[1].PruneLen) / (tree.Chs[0].PruneLen + tree.Chs[1].PruneLen))
			tmpChar := ((tree.Chs[0].PruneLen * tree.Chs[1].ContData[i]) + (tree.Chs[1].PruneLen * tree.Chs[0].ContData[i])) / curVar
			contrast = tmpChar - tree.Chs[2].ContData[i]
			curVar = tree.Chs[2].PruneLen + tmpPruneLen
			tmpll += ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
		}
		sitelikes += tmpll
	}
	chll += sitelikes
	return
}

//CalcRootedLogLike will return the BM likelihood of a tree assuming that no data are missing from the tips.
func CalcRootedLogLike(n *gophy.Node, nlikes *float64, startFresh bool) {
	for _, chld := range n.Chs {
		CalcRootedLogLike(chld, nlikes, startFresh)
	}
	nchld := len(n.Chs)
	n.PruneLen = n.BMLen
	if nchld != 0 {
		if nchld != 2 {
			fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
		}
		c0 := n.Chs[0]
		c1 := n.Chs[1]
		curlike := float64(0.0)
		var tempChar float64
		for i := range n.Chs[0].ContData {
			curVar := c0.PruneLen + c1.PruneLen
			contrast := c0.ContData[i] - c1.ContData[i]
			curlike += ((-0.5) * ((math.Log(2 * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2) / (curVar))))
			tempChar = ((c0.PruneLen * c1.ContData[i]) + (c1.PruneLen * c0.ContData[i])) / (curVar)
			n.ContData[i] = tempChar
		}
		*nlikes += curlike
		tempBranchLength := n.BMLen + ((c0.PruneLen * c1.PruneLen) / (c0.PruneLen + c1.PruneLen)) // need to calculate the prune length by adding the averaged lengths of the daughter nodes to the length
		n.PruneLen = tempBranchLength                                                             // need to calculate the "prune length" by adding the length to the uncertainty

	}

}

//calcUnrootedNodeLikes  will calculate the likelihood of an unrooted tree at each site (i) of the continuous character alignment
func calcUnrootedSiteLLParallel(tree *gophy.Node, i int) (tmpll float64) {
	var contrast, curVar float64
	log2pi := 1.8378770664093453
	if tree.Chs[0].Mis[i] == false && tree.Chs[1].Mis[i] == false && tree.Chs[2].Mis[i] == false { //do the standard calculation when no subtrees have missing traits
		contrast = tree.Chs[0].ContData[i] - tree.Chs[1].ContData[i]
		curVar = tree.Chs[0].ConPruneLen[i] + tree.Chs[1].ConPruneLen[i]
		tmpll = ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
		tmpConPruneLen := ((tree.Chs[0].ConPruneLen[i] * tree.Chs[1].ConPruneLen[i]) / (tree.Chs[0].ConPruneLen[i] + tree.Chs[1].ConPruneLen[i]))
		tmpChar := ((tree.Chs[0].ConPruneLen[i] * tree.Chs[1].ContData[i]) + (tree.Chs[1].ConPruneLen[i] * tree.Chs[0].ContData[i])) / curVar
		contrast = tmpChar - tree.Chs[2].ContData[i]
		curVar = tree.Chs[2].ConPruneLen[i] + tmpConPruneLen
		tmpll += ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	} else if tree.Chs[0].Mis[i] == false && tree.Chs[1].Mis[i] == false && tree.Chs[2].Mis[i] == true { // do standard "rooted" calculation on Chs[0] and Chs [1] if Chs[2] is missing
		contrast = tree.Chs[0].ContData[i] - tree.Chs[1].ContData[i]
		curVar = tree.Chs[0].ConPruneLen[i] + tree.Chs[1].ConPruneLen[i]
		tmpll = ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	} else if tree.Chs[0].Mis[i] == false && tree.Chs[2].Mis[i] == false && tree.Chs[1].Mis[i] == true { // do standard "rooted" calculation on Chs[0] and Chs [2] if Chs[1] is missing
		contrast = tree.Chs[0].ContData[i] - tree.Chs[2].ContData[i]
		curVar = tree.Chs[0].ConPruneLen[i] + tree.Chs[2].ConPruneLen[i]
		tmpll = ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	} else if tree.Chs[1].Mis[i] == false && tree.Chs[2].Mis[i] == false && tree.Chs[0].Mis[i] == true { // do standard "rooted" calculation on Chs[1] and Chs [2] if Chs[0] is missing
		contrast = tree.Chs[1].ContData[i] - tree.Chs[2].ContData[i]
		curVar = tree.Chs[1].ConPruneLen[i] + tree.Chs[2].ConPruneLen[i]
		tmpll = ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	}
	return
}

//calcRootedSiteLL will return the BM likelihood of a tree assuming that no data are missing from the tips.
func calcRootedSiteLLParallel(n *gophy.Node, nlikes *float64, startFresh bool, site int) {
	for _, chld := range n.Chs {
		calcRootedSiteLLParallel(chld, nlikes, startFresh, site)
	}
	nchld := len(n.Chs)
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
		} else if c0.Mis[site] == true && c1.Mis[site] == false {
			n.ConPruneLen[site] += (c1.ConPruneLen[site])
			n.ContData[site] = c1.ContData[site]
		} else if c1.Mis[site] == true && c0.Mis[site] == false {
			n.ConPruneLen[site] += c0.ConPruneLen[site]
			n.ContData[site] = c0.ContData[site]
		} else if c1.Mis[site] == true && c0.Mis[site] == true {
			n.Mis[site] = true
		}
	}
}

func siteTreeLikeParallel(tree, ch1, ch2, ch3 *gophy.Node, startFresh bool, weights []float64, jobs <-chan int, results chan<- float64) {
	for site := range jobs {
		tmpll := 0.
		rootedMissingSiteLL(ch1, &tmpll, true, site)
		rootedMissingSiteLL(ch2, &tmpll, true, site)
		rootedMissingSiteLL(ch3, &tmpll, true, site)
		tmpll += calcUnrootedSiteLLParallel(tree, site)
		tmpll = tmpll * weights[site]
		results <- tmpll
	}
}

//WeightedUnrootedLogLikeParallel will calculate the log-likelihood of an unrooted tree, while assuming that some sites have missing data. This can be used to calculate the likelihoods of trees that have complete trait sampling, but it will be slower than CalcRootedLogLike.
func WeightedUnrootedLogLikeParallel(tree *gophy.Node, startFresh bool, weights []float64, workers int) (sitelikes float64) {
	nsites := len(tree.Chs[0].ContData)
	ch1 := tree.Chs[0] //.PostorderArray()
	ch2 := tree.Chs[1] //.PostorderArray()
	ch3 := tree.Chs[2] //.PostorderArray()
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	for w := 0; w < workers; w++ {
		go siteTreeLikeParallel(tree, ch1, ch2, ch3, startFresh, weights, jobs, results)
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

func subSiteTreeLikeParallel(tree, ch1, ch2, ch3 *gophy.Node, startFresh bool, jobs <-chan int, results chan<- float64) {
	for site := range jobs {
		tmpll := 0.
		//calcRootedSiteLLParallel(ch1, &tmpll, startFresh, site)
		//calcRootedSiteLLParallel(ch2, &tmpll, startFresh, site)
		//calcRootedSiteLLParallel(ch3, &tmpll, startFresh, site)
		rootedMissingSiteLL(ch1, &tmpll, true, site)
		rootedMissingSiteLL(ch2, &tmpll, true, site)
		rootedMissingSiteLL(ch3, &tmpll, true, site)
		tmpll += calcUnrootedSiteLLParallel(tree, site)
		results <- tmpll
	}
}

//SubUnrootedLogLikeParallel will calculate the log-likelihood of an unrooted tree, while assuming that some sites have missing data. This can be used to calculate the likelihoods of trees that have complete trait sampling, but it will be slower than CalcRootedLogLike.
func SubUnrootedLogLikeParallel(tree *gophy.Node, sites []int, workers int) (sitelikes float64) {
	nsites := len(sites)
	ch1 := tree.Chs[0] //.PostorderArray()
	ch2 := tree.Chs[1] //.PostorderArray()
	ch3 := tree.Chs[2] //.PostorderArray()
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	for w := 0; w < workers; w++ {
		go subSiteTreeLikeParallel(tree, ch1, ch2, ch3, true, jobs, results)
	}
	//for site := 0; site < nsites; site++ {
	for _, site := range sites {
		jobs <- site
	}
	close(jobs)

	for range sites {
		sitelikes += <-results
	}
	return
}

//SingleSiteLL will return the likelihood of a single site
func SingleSiteLL(tree *gophy.Node, site int) (sitelike float64) {
	sitelike = 0.0
	var tmpll float64
	if len(tree.Chs) == 3 {
		ch1 := tree.Chs[0] //.PostorderArray()
		ch2 := tree.Chs[1] //.PostorderArray()
		ch3 := tree.Chs[2] //.PostorderArray()
		tmpll = 0.
		rootedMissingSiteLL(ch1, &tmpll, true, site)
		rootedMissingSiteLL(ch2, &tmpll, true, site)
		rootedMissingSiteLL(ch3, &tmpll, true, site)
		tmpll += calcUnrootedSiteLLParallel(tree, site)
		//fmt.Println(tmpll)
	} else if len(tree.Chs) == 2 {
		tmpll = 0.
		rootedMissingSiteLL(tree, &tmpll, true, site)
	}
	sitelike = tmpll
	return
}

func rootedMissingSiteLL(n *gophy.Node, nlikes *float64, startFresh bool, site int) {
	for _, chld := range n.Chs {
		rootedMissingSiteLL(chld, nlikes, startFresh, site)
	}
	nchld := len(n.Chs)
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
		} else if c0.Mis[site] == true && c1.Mis[site] == false {
			n.ConPruneLen[site] += (c1.ConPruneLen[site])
			n.ContData[site] = c1.ContData[site]
		} else if c1.Mis[site] == true && c0.Mis[site] == false {
			n.ConPruneLen[site] += c0.ConPruneLen[site]
			n.ContData[site] = c0.ContData[site]
		} else if c1.Mis[site] == true && c0.Mis[site] == true {
			n.Mis[site] = true
		}
	}

}

func rootedTreeLike(tree *gophy.Node, startFresh bool, jobs <-chan int, results chan<- float64) {
	for site := range jobs {
		tmpll := 0.
		rootedMissingSiteLL(tree, &tmpll, startFresh, site)
		results <- tmpll
	}
}

//PBMLogLikeRt will calculate the BM log like on a rooted tree
func PBMLogLikeRt(tree *gophy.Node, startFresh bool, workers int) (sitelikes float64) {
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
