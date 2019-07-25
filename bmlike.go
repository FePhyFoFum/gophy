package gophy

import (
	"fmt"
	"math"
	"os"
)

func postorder(curnode *Node) {
	for _, chld := range curnode.Chs {
		postorder(chld)
	}
	fmt.Println(curnode.Nam, curnode.ContData)
}

//BMPruneRooted will prune BM branch lens and PICs down to a rooted node
//root node should be a real (ie. bifurcating) root
func BMPruneRooted(n *Node) {
	for _, chld := range n.Chs {
		BMPruneRooted(chld)
	}
	n.PruneLen = n.Len
	nchld := len(n.Chs)
	if nchld != 0 { //&& n.Marked == false {
		var tempChar float64
		if nchld != 2 {
			fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
		}
		c0 := n.Chs[0]
		c1 := n.Chs[1]
		bot := ((1.0 / c0.PruneLen) + (1.0 / c1.PruneLen))
		n.PruneLen += 1.0 / bot
		for i := range n.Chs[0].ContData {
			tempChar = (((1 / c0.PruneLen) * c1.ContData[i]) + ((1 / c1.PruneLen) * c0.ContData[i])) / bot
			n.ContData[i] = tempChar
		}
	}
}

//AncBMPruneRooted will prune BM branch lens and PICs down to a rooted node
//root node should be a real (ie. bifurcating) root
func AncBMPruneRooted(n *Node) {
	for _, chld := range n.Chs {
		AncBMPruneRooted(chld)
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
		n.PruneLen += 1.0 / bot
		for i := range n.Chs[0].ContData {
			tempChar = (((1 / c0.PruneLen) * c1.ContData[i]) + ((1 / c1.PruneLen) * c0.ContData[i])) / bot
			n.ContData[i] = tempChar
		}
	}
}

//AncBMPruneRootedSingle will prune BM branch lens and calculate PIC of a single trait down to a rooted node
//root node should be a real (ie. bifurcating) root
func AncBMPruneRootedSingle(n *Node, i int) {
	for _, chld := range n.Chs {
		AncBMPruneRootedSingle(chld, i)
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

//BMPruneRootedSingle will prune BM branch lens and calculate PIC of a single trait down to a rooted node
//root node should be a real (ie. bifurcating) root
func BMPruneRootedSingle(n *Node, i int) {
	for _, chld := range n.Chs {
		BMPruneRootedSingle(chld, i)
	}
	n.PruneLen = n.Len
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

func debugParChld(tree *Node) {
	for _, ch := range tree.Chs {
		fmt.Println(ch.Nam, "\t")
	}
}

//AssertUnrootedTree is a quick check to make sure the tree passed is unrooted
func AssertUnrootedTree(tree *Node) {
	if len(tree.Chs) != 3 {
		fmt.Print("BRAncH LenGTHS MUST BE ITERATED ON AN UNROOTED TREE. THIS TREE IS ROOTED.")
		os.Exit(0)
	}
}

//IterateBMLengths will iteratively calculate the ML branch lengths for a particular topology, assuming that traits are fully sampled at the tips. should use MissingTraitsEM if there are missing sites.
func IterateBMLengths(tree *Tree, niter int) {
	AssertUnrootedTree(tree.Rt)
	itercnt := 0
	for {
		calcBMLengths(tree)
		itercnt++
		if itercnt == niter {
			break
		}
	}
}

//AncMissingTraitsEM will iteratively calculate the ML branch lengths for a particular topology
func AncMissingTraitsEM(tree *Tree, niter int) (nparam float64) {
	AssertUnrootedTree(tree.Rt)
	nodes := tree.Pre
	InitMissingValues(nodes)
	itercnt := 0
	for {
		AncCalcExpectedTraits(tree.Rt) //calculate Expected trait values
		ancCalcBMLengths(tree)         //maximize likelihood of branch lengths
		itercnt++
		if itercnt == niter {
			break
		}
	}
	nparam = 0.0
	for _, n := range nodes {
		if n.Anc == true && len(n.Chs) == 0 {
			continue
		}
		nparam += 1.0
	}
	return
}

//MissingTraitsEM will iteratively calculate the ML branch lengths for a particular topology
func MissingTraitsEM(tree *Tree, niter int) {
	AssertUnrootedTree(tree.Rt)
	nodes := tree.Pre
	InitMissingValues(nodes)
	itercnt := 0
	for {
		CalcExpectedTraits(tree.Rt) //calculate Expected trait values
		calcBMLengths(tree)         //maximize likelihood of branch lengths
		itercnt++
		if itercnt == niter {
			break
		}
	}
}

//ancCalcBMLengths will perform a single pass of the branch length ML estimation
func ancCalcBMLengths(tree *Tree) {
	rnodes := tree.Pre
	lnode := 0
	for ind, newroot := range rnodes {
		if len(newroot.Chs) == 0 {
			continue
		} else if newroot != rnodes[0] {
			newroot.RerootBM(rnodes[lnode])
			lnode = ind
		}
		for _, cn := range tree.Rt.Chs {
			AncBMPruneRooted(cn)
		}
		AncTritomyML(tree.Rt)
	}
	rnodes[0].RerootBM(tree.Rt)
	//fmt.Println(tree.Newick(true))
}

//calcBMLengths will perform a single pass of the branch length ML estimation
func calcBMLengths(tree *Tree) {
	rnodes := tree.Pre
	lnode := 0
	for ind, newroot := range rnodes {
		if len(newroot.Chs) == 0 {
			continue
		} else if newroot != rnodes[0] {
			newroot.Reroot(rnodes[lnode])
			lnode = ind
		}
		for _, cn := range tree.Rt.Chs {
			AncBMPruneRooted(cn)
		}
		TritomyML(tree.Rt)
	}
	rnodes[0].Reroot(tree.Rt)
	//fmt.Println(tree.Newick(true))
}

//PruneToStar will prune brlens and traits to a root
func PruneToStar(tree *Node) {
	for _, cn := range tree.Chs {
		BMPruneRooted(cn)
	}
}

//CalcUnrootedLogLike will calculate the log-likelihood of an unrooted tree, while assuming that no sites have missing data.
func CalcUnrootedLogLike(tree *Node, startFresh bool) (chll float64) {
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
func CalcRootedLogLike(n *Node, nlikes *float64, startFresh bool) {
	for _, chld := range n.Chs {
		CalcRootedLogLike(chld, nlikes, startFresh)
	}
	nchld := len(n.Chs)
	if n.Marked == true && startFresh == false {
		if nchld != 0 {
			for _, l := range n.LL {
				*nlikes += l
			}
		}
	} else if n.Marked == false || startFresh == true {
		n.PruneLen = n.Len
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
				n.LL[i] = curlike
				tempChar = ((c0.PruneLen * c1.ContData[i]) + (c1.PruneLen * c0.ContData[i])) / (curVar)
				n.ContData[i] = tempChar
			}
			*nlikes += curlike
			tempBranchLength := n.Len + ((c0.PruneLen * c1.PruneLen) / (c0.PruneLen + c1.PruneLen)) // need to calculate the prune length by adding the averaged lengths of the daughter nodes to the length
			n.PruneLen = tempBranchLength                                                           // need to calculate the "prune length" by adding the length to the uncertainty

			n.Marked = true
		}
	}
}

//SitewiseLogLike will calculate the log-likelihood of an unrooted tree, while assuming that some sites have missing data. This can be used to calculate the likelihoods of trees that have complete trait sampling, but it will be slower than CalcRootedLogLike.
func SitewiseLogLike(tree *Node) (sitelikes []float64) {
	var tmpll float64
	ch1 := tree.Chs[0]                       //.PostorderArray()
	ch2 := tree.Chs[1]                       //.PostorderArray()
	ch3 := tree.Chs[2]                       //.PostorderArray()
	for site := range tree.Chs[0].ContData { //calculate log likelihood at each site
		tmpll = 0.
		calcRootedSiteLL(ch1, &tmpll, true, site)
		calcRootedSiteLL(ch2, &tmpll, true, site)
		calcRootedSiteLL(ch3, &tmpll, true, site)
		tmpll += calcUnrootedSiteLL(tree, site)
		//fmt.Println(site, tmpll)
		sitelikes = append(sitelikes, tmpll)
	}
	return
}

/*/WeightedUnrootedLogLike will calculate the log-likelihood of an unrooted tree, while assuming that no sites have missing data.
func WeightedUnrootedLogLike(tree *Node, startFresh bool, weights []float64) (chll float64) {
	chll = 0.0
	for _, ch := range tree.Chs {
		n := ch.PostorderArray()
		for i := range ch.ContData {
			curlike := calcRootedSiteLL(n, startFresh, i)
			chll += curlike
		}
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
		sitelikes += tmpll * weights[i]
	}
	chll += sitelikes
	return
}
*/

//MarkAll will mark all of the nodes in a tree ARRAY
func MarkAll(nodes []*Node) {
	for _, n := range nodes {
		n.Marked = true
	}
}

//calcUnrootedNodeLikes  will calculate the likelihood of an unrooted tree at each site (i) of the continuous character alignment
func calcUnrootedSiteLLParallel(tree *Node, i int) (tmpll float64) {
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
func calcRootedSiteLLParallel(n *Node, nlikes *float64, startFresh bool, site int) {
	for _, chld := range n.Chs {
		calcRootedSiteLLParallel(chld, nlikes, startFresh, site)
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
		n.ConPruneLen[site] = n.Len
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

/*/calcRootedSiteLL will return the BM likelihood of a tree assuming that no data are missing from the tips.
func calcRootedSiteLLParallel(n *Node, nlikes *float64, startFresh bool, site int) {
	for _, chld := range n.Chs {
		calcRootedSiteLLParallel(chld, nlikes, startFresh, site)
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
		n.ConPruneLen[site] = n.Len
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
			curVar := c0.ConPruneLen[site] + c1.ConPruneLen[site]
			contrast := c0.ContData[site] - c1.ContData[site]
			curlike += ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2) / (curVar))))
			tempChar = ((c0.ConPruneLen[site] * c1.ContData[site]) + (c1.ConPruneLen[site] * c0.ContData[site])) / (curVar)
			n.ContData[site] = tempChar
			*nlikes += curlike
			tempBranchLength := n.ConPruneLen[site] + ((c0.ConPruneLen[site] * c1.ConPruneLen[site]) / (c0.ConPruneLen[site] + c1.ConPruneLen[site])) // need to calculate the prune length by adding the averaged lengths of the daughter nodes to the length
			n.ConPruneLen[site] = tempBranchLength                                                                                            // need to calculate the "prune length" by adding the length to the uncertainty
			n.LL[site] = curlike
			//n.Marked = true
		}
	}
}*/

func siteTreeLikeParallel(tree, ch1, ch2, ch3 *Node, startFresh bool, weights []float64, jobs <-chan int, results chan<- float64) {
	for site := range jobs {
		tmpll := 0.
		calcRootedSiteLLParallel(ch1, &tmpll, startFresh, site)
		calcRootedSiteLLParallel(ch2, &tmpll, startFresh, site)
		calcRootedSiteLLParallel(ch3, &tmpll, startFresh, site)
		tmpll += calcUnrootedSiteLLParallel(tree, site)
		tmpll = tmpll * weights[site]
		results <- tmpll
	}
}

//WeightedUnrootedLogLikeParallel will calculate the log-likelihood of an unrooted tree, while assuming that some sites have missing data. This can be used to calculate the likelihoods of trees that have complete trait sampling, but it will be slower than CalcRootedLogLike.
func WeightedUnrootedLogLikeParallel(tree *Node, startFresh bool, weights []float64, workers int) (sitelikes float64) {
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

//WeightedUnrootedLogLike will calculate the log-likelihood of an unrooted tree, while assuming that some sites have missing data. This can be used to calculate the likelihoods of trees that have complete trait sampling, but it will be slower than CalcRootedLogLike.
func WeightedUnrootedLogLike(tree *Node, startFresh bool, weights []float64) (sitelikes float64) {
	sitelikes = 0.0
	var tmpll float64
	ch1 := tree.Chs[0]                       //.PostorderArray()
	ch2 := tree.Chs[1]                       //.PostorderArray()
	ch3 := tree.Chs[2]                       //.PostorderArray()
	for site := range tree.Chs[0].ContData { //calculate log likelihood at each site
		tmpll = 0.
		calcRootedSiteLL(ch1, &tmpll, startFresh, site)
		calcRootedSiteLL(ch2, &tmpll, startFresh, site)
		calcRootedSiteLL(ch3, &tmpll, startFresh, site)
		tmpll += calcUnrootedSiteLL(tree, site)
		tmpll = tmpll * weights[site]
		//fmt.Println(tmpll)
		sitelikes += tmpll
	}
	return
}

//MissingUnrootedLogLike will calculate the log-likelihood of an unrooted tree. It is deprecated and is only a convienience for testing LL calculation.
func MissingUnrootedLogLike(tree *Node, startFresh bool) (sitelikes float64) {
	sitelikes = 0.0
	var tmpll float64
	ch1 := tree.Chs[0]                       //.PostorderArray()
	ch2 := tree.Chs[1]                       //.PostorderArray()
	ch3 := tree.Chs[2]                       //.PostorderArray()
	for site := range tree.Chs[0].ContData { //calculate log likelihood at each site
		tmpll = 0.
		calcRootedSiteLL(ch1, &tmpll, startFresh, site)
		calcRootedSiteLL(ch2, &tmpll, startFresh, site)
		calcRootedSiteLL(ch3, &tmpll, startFresh, site)
		tmpll += calcUnrootedSiteLL(tree, site)
		sitelikes += tmpll
	}
	return
}

//calcUnrootedNodeLikes  will calculate the likelihood of an unrooted tree at each site (i) of the continuous character alignment
func calcUnrootedSiteLL(tree *Node, i int) (tmpll float64) {
	var contrast, curVar float64
	if tree.Chs[0].Mis[i] == false && tree.Chs[1].Mis[i] == false && tree.Chs[2].Mis[i] == false { //do the standard calculation when no subtrees have missing traits
		contrast = tree.Chs[0].ContData[i] - tree.Chs[1].ContData[i]
		curVar = tree.Chs[0].PruneLen + tree.Chs[1].PruneLen
		tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
		tmpPruneLen := ((tree.Chs[0].PruneLen * tree.Chs[1].PruneLen) / (tree.Chs[0].PruneLen + tree.Chs[1].PruneLen))
		tmpChar := ((tree.Chs[0].PruneLen * tree.Chs[1].ContData[i]) + (tree.Chs[1].PruneLen * tree.Chs[0].ContData[i])) / curVar
		contrast = tmpChar - tree.Chs[2].ContData[i]
		curVar = tree.Chs[2].PruneLen + tmpPruneLen
		tmpll += ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	} else if tree.Chs[0].Mis[i] == false && tree.Chs[1].Mis[i] == false && tree.Chs[2].Mis[i] == true { // do standard "rooted" calculation on Chs[0] and Chs [1] if Chs[2] is missing
		contrast = tree.Chs[0].ContData[i] - tree.Chs[1].ContData[i]
		curVar = tree.Chs[0].PruneLen + tree.Chs[1].PruneLen
		tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	} else if tree.Chs[0].Mis[i] == false && tree.Chs[2].Mis[i] == false && tree.Chs[1].Mis[i] == true { // do standard "rooted" calculation on Chs[0] and Chs [2] if Chs[1] is missing
		contrast = tree.Chs[0].ContData[i] - tree.Chs[2].ContData[i]
		curVar = tree.Chs[0].PruneLen + tree.Chs[2].PruneLen
		tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	} else if tree.Chs[1].Mis[i] == false && tree.Chs[2].Mis[i] == false && tree.Chs[0].Mis[i] == true { // do standard "rooted" calculation on Chs[1] and Chs [2] if Chs[0] is missing
		contrast = tree.Chs[1].ContData[i] - tree.Chs[2].ContData[i]
		curVar = tree.Chs[1].PruneLen + tree.Chs[2].PruneLen
		tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	}
	return
}

//MissingRootedLogLike will return the BM log-likelihood of a tree for a single site, pruning tips that have missing data
func MissingRootedLogLike(n *Node, startFresh bool) (sitelikes float64) {
	//nodes := n.PostorderArray()
	sitelikes = 0.
	for site := range n.ContData {
		calcRootedSiteLL(n, &sitelikes, startFresh, site)
	}
	return
}

//calcRootedSiteLL will return the BM likelihood of a tree assuming that no data are missing from the tips.
func calcRootedSiteLL(n *Node, nlikes *float64, startFresh bool, site int) {
	for _, chld := range n.Chs {
		calcRootedSiteLL(chld, nlikes, startFresh, site)
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
		n.PruneLen = n.Len
		if nchld != 0 {
			if nchld != 2 {
				fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
				os.Exit(0)
			}
			c0 := n.Chs[0]
			c1 := n.Chs[1]
			curlike := float64(0.0)
			var tempChar float64
			curVar := c0.PruneLen + c1.PruneLen
			contrast := c0.ContData[site] - c1.ContData[site]
			curlike += ((-0.5) * ((math.Log(2 * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2) / (curVar))))
			tempChar = ((c0.PruneLen * c1.ContData[site]) + (c1.PruneLen * c0.ContData[site])) / (curVar)
			n.ContData[site] = tempChar
			*nlikes += curlike
			tempBranchLength := n.Len + ((c0.PruneLen * c1.PruneLen) / (c0.PruneLen + c1.PruneLen)) // need to calculate the prune length by adding the averaged lengths of the daughter nodes to the length
			n.PruneLen = tempBranchLength                                                           // need to calculate the "prune length" by adding the length to the uncertainty
			n.LL[site] = curlike
			//n.Marked = true
		}
	}
}

//TritomyML will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree
func TritomyML(tree *Node) {
	ntraits := len(tree.Chs[0].ContData)
	fntraits := float64(ntraits)
	var x1, x2, x3 float64
	sumV1 := 0.0
	sumV2 := 0.0
	sumV3 := 0.0
	for i := range tree.Chs[0].ContData {
		x1 = tree.Chs[0].ContData[i]
		x2 = tree.Chs[1].ContData[i]
		x3 = tree.Chs[2].ContData[i]
		sumV1 += ((x1 - x2) * (x1 - x3))
		sumV2 += ((x2 - x1) * (x2 - x3))
		sumV3 += ((x3 - x1) * (x3 - x2))
	}
	if sumV1 < 0.0 {
		sumV1 = 0.000001
		sumV2 = 0.0
		sumV3 = 0.0
		for i := range tree.Chs[0].ContData {
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.Chs[2].ContData[i]
			sumV2 += (x1 - x2) * (x1 - x2)
			sumV3 += (x1 - x3) * (x1 - x3)
		}
	} else if sumV2 < 0.0 {
		sumV1 = 0.0
		sumV2 = 0.00001
		sumV3 = 0.0
		for i := range tree.Chs[0].ContData {
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.Chs[2].ContData[i]
			sumV1 += (x2 - x1) * (x2 - x1)
			sumV3 += (x2 - x3) * (x2 - x3)
		}
	} else if sumV3 < 0.0 {
		sumV1 = 0.0
		sumV2 = 0.0
		sumV3 = 0.0001
		for i := range tree.Chs[0].ContData {
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.Chs[2].ContData[i]
			sumV1 += (x3 - x1) * (x3 - x1)
			sumV2 += (x3 - x2) * (x3 - x2)
		}
	}
	sumV1 = sumV1 / fntraits
	sumV2 = sumV2 / fntraits
	sumV3 = sumV3 / fntraits
	sumV1 = sumV1 - (tree.Chs[0].PruneLen - tree.Chs[0].Len)
	sumV2 = sumV2 - (tree.Chs[1].PruneLen - tree.Chs[1].Len)
	sumV3 = sumV3 - (tree.Chs[2].PruneLen - tree.Chs[2].Len)
	if sumV1 < 0. {
		sumV1 = 0.0001
	}
	if sumV2 < 0. {
		sumV2 = 0.0001
	}
	if sumV3 < 0. {
		sumV3 = 0.0001
	}
	tree.Chs[0].Len = sumV1
	tree.Chs[1].Len = sumV2
	tree.Chs[2].Len = sumV3
}

//AncTritomyML will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree assuming that direct ancestors may be in the tree
func AncTritomyML(tree *Node) {
	ntraits := len(tree.Chs[0].ContData)
	fntraits := float64(ntraits)
	var x1, x2, x3 float64
	sumV1 := 0.0
	sumV2 := 0.0
	sumV3 := 0.0
	for i := range tree.Chs[0].ContData {
		x1 = tree.Chs[0].ContData[i]
		x2 = tree.Chs[1].ContData[i]
		x3 = tree.Chs[2].ContData[i]
		sumV1 += ((x1 - x2) * (x1 - x3))
		sumV2 += ((x2 - x1) * (x2 - x3))
		sumV3 += ((x3 - x1) * (x3 - x2))
	}
	if sumV1 < 0.0 || tree.Chs[0].Anc == true && len(tree.Chs[0].Chs) == 0 {
		sumV1 = 0.0000000000001
		sumV2 = 0.0
		sumV3 = 0.0
		for i := range tree.Chs[0].ContData {
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.Chs[2].ContData[i]
			sumV2 += (x1 - x2) * (x1 - x2)
			sumV3 += (x1 - x3) * (x1 - x3)
		}
	} else if sumV2 < 0.0 || tree.Chs[1].Anc == true && len(tree.Chs[1].Chs) == 0 {
		sumV1 = 0.0
		sumV2 = 0.0000000000001 //0.0
		sumV3 = 0.0
		for i := range tree.Chs[0].ContData {
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.Chs[2].ContData[i]
			sumV1 += (x2 - x1) * (x2 - x1)
			sumV3 += (x2 - x3) * (x2 - x3)
		}
	} else if sumV3 < 0.0 || tree.Chs[2].Anc == true && len(tree.Chs[2].Chs) == 0 {
		sumV1 = 0.0
		sumV2 = 0.0
		sumV3 = 0.0000000000001 //0.0
		for i := range tree.Chs[0].ContData {
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.Chs[2].ContData[i]
			sumV1 += (x3 - x1) * (x3 - x1)
			sumV2 += (x3 - x2) * (x3 - x2)
		}
	}
	if sumV1 != 0.0 {
		sumV1 = sumV1 / fntraits
		sumV1 = sumV1 - (tree.Chs[0].PruneLen - tree.Chs[0].BMLen)
	}
	if sumV2 != 0.0 {
		sumV2 = sumV2 / fntraits
		sumV2 = sumV2 - (tree.Chs[1].PruneLen - tree.Chs[1].BMLen)
	}
	if sumV3 != 0.0 {
		sumV3 = sumV3 / fntraits
		sumV3 = sumV3 - (tree.Chs[2].PruneLen - tree.Chs[2].BMLen)
	}
	if sumV1 < 0. {
		sumV1 = 0.0
	}
	if sumV2 < 0. {
		sumV2 = 0.0
	}
	if sumV3 < 0. {
		sumV3 = 0.0
	}
	tree.Chs[0].BMLen = sumV1
	tree.Chs[1].BMLen = sumV2
	tree.Chs[2].BMLen = sumV3
}
