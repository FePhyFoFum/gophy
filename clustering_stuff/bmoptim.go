package clustering

import (
	"fmt"
	"math"
	"math/rand"
	"os"
)

//ClusterMissingTraitsEM will calculate the BM branch lengths using an iterative EM calculation that imputes missing data using PICs using the traits in a single cluster
func ClusterMissingTraitsEM(t *Tree, cluster *Cluster, niter int) {
	best := -1000000000.0
	var newlen []float64
	for it := 0; it < 1; it++ {
		/*for _, n := range t.Pre[1:] {
			r := rand.Float64()
			n.BMLen = r
		}*/
		for i := 0; i < niter; i++ {
			CalcExpectedTraits(t.Rt)
			BMCalcLensBackFront(t, cluster.Sites)
		}
		if len(t.Rt.Chs) == 3 {
			ll := SubUnrootedLogLikeParallel(t.Rt, cluster.Sites, 4)
			if ll > best {
				best = ll
				newlen = nil
				for _, n := range t.Pre {
					newlen = append(newlen, n.BMLen) //store the newly calculated branch lengths
				}
			}
		} else {
			fmt.Println("picard currently only runs on unrooted trees.")
		}
	}
	for i, n := range t.Pre {
		n.BMLen = newlen[i]
	}
	cluster.BranchLengths = newlen
}

//GreedyIterateLengthsMissing will calculate the BM branch lengths using an iterative EM calculation that imputes missing data using PICs using the traits in a single cluster
func GreedyIterateLengthsMissing(t *Tree, sites []int, niter int) {
	best := -1000000000.0
	var newlen []float64
	for it := 0; it < 1; it++ {
		/*for _, n := range t.Pre[1:] {
			r := rand.Float64()
			n.BMLen = r
		}*/
		for i := 0; i < niter; i++ {
			CalcExpectedTraits(t.Rt)
			BMCalcLensBackFront(t, sites)
		}
		if len(t.Rt.Chs) == 3 {
			ll := SubUnrootedLogLikeParallel(t.Rt, sites, 4)
			if ll > best {
				best = ll
				for _, n := range t.Pre {
					newlen = append(newlen, n.BMLen) //store the newly calculated branch lengths
				}
			}
		} else {
			ll := PBMLogLikeRt(t.Rt, true, 4)
			if ll > best {
				best = ll
				newlen = nil
				for _, n := range t.Pre {
					newlen = append(newlen, n.BMLen) //store the newly calculated branch lengths
				}

			}
		}
	}
	for i, n := range t.Pre {
		n.BMLen = newlen[i]
	}
}

//BMOptimBLEM will calculate the BM branch lengths using an iterative EM calculation that imputes missing data using PICs
func BMOptimBLEM(t *Tree, niter int) {
	var sites []int
	for i := range t.Rt.ContData {
		sites = append(sites, i)
	}
	best := -1000000000.0
	var newlen []float64
	for it := 0; it < 10; it++ {
		for _, n := range t.Pre[1:] {
			r := rand.Float64()
			n.BMLen = r
		}
		for i := 0; i < niter; i++ {
			CalcExpectedTraits(t.Rt)
			BMCalcLensBackFront(t, sites)
		}
		if len(t.Rt.Chs) == 3 {
			ll := SubUnrootedLogLikeParallel(t.Rt, sites, 4)
			if ll > best {
				best = ll
				newlen = nil
				for _, n := range t.Pre {
					newlen = append(newlen, n.BMLen) //store the newly calculated branch lengths
				}
			}
		} else {
			ll := PBMLogLikeRt(t.Rt, true, 4)
			if ll > best {
				best = ll
				newlen = nil
				for _, n := range t.Pre {
					newlen = append(newlen, n.BMLen) //store the newly calculated branch lengths
				}
				for i, n := range t.Pre {
					n.BMLen = newlen[i]
				}
				//fmt.Println("MID:", PBMLogLikeRt(t.Rt, true, 4), best, ll)
				//fmt.Println(t.Rt.BMPhylogram())
			}
		}

	}
	for i, n := range t.Pre {
		n.BMLen = newlen[i]
	}
	//CalcExpectedTraits(t.Rt)
	//ll := PBMLogLikeRt(t.Rt, true, 4)
	//fmt.Println("DONE:", PBMLogLikeRt(t.Rt, true, 4), best)

}

//BMCalcLensBackFront will do one pass of the EM branch length estimation
func BMCalcLensBackFront(t *Tree, sites []int) {
	if len(t.Rt.Chs) == 3 { // the tree is unrooted
		//parSubtreePruneLen := make(map[*Node]float64)
		for _, c := range t.Rt.Chs {
			BMPruneRooted(c)
		}
		AncTritomyML(t.Rt, sites)
		for _, c := range t.Rt.Chs {
			BMPruneRooted(c)
		}
	} else {
		BMPruneRooted(t.Rt)
	}
	for _, c := range t.Rt.Chs {
		for _, nd := range c.PreorderArray() {
			//fmt.Println(nd.BMLen, nd.PruneLen, nd.Nam)
			if len(nd.Chs) == 0 { //don't do this at tips
				continue
			}
			var c0, c1 *Node
			var bot float64
			if nd.Par == t.Rt && len(t.Rt.Chs) == 3 {
				var sibs []*Node
				for _, c := range t.Rt.Chs {
					if c != nd {
						sibs = append(sibs, c)
					}
				}
				c0 = sibs[0]
				c1 = sibs[1]
			} else if nd.Par == t.Rt && len(t.Rt.Chs) == 2 {
				sib := nd.GetSib()
				if len(sib.Chs) == 0 { // special case if the tree is rooted on a single taxon
					nd.BMLen += sib.BMLen
					nd.ContData = sib.ContData
					nd.PruneLen = nd.BMLen
					virtualTritomyML(nd, sites)
					fmt.Println("WARNING: BM Branch length optimization is currently wonky on trees rooted on a single taxon")
					//os.Exit(0)
					nd.PruneLen = nd.BMLen
					sib.BMLen = nd.BMLen / 2
					nd.BMLen = nd.BMLen / 2
					continue // don't do the stuff below; move on to the next node
				} else {
					c0 = sib.Chs[0]
					c1 = sib.Chs[1]
					nd.BMLen += sib.BMLen
				}
			} else if nd.Par != t.Rt {
				c0 = nd.GetSib()
				c1 = nd.Par
				/*if c0.PruneLen == 0 {
					c0.PruneLen = 0.00001
				}*/
			}
			bot = ((1.0 / c0.PruneLen) + (1.0 / c1.PruneLen))
			//fmt.Println("HERE", c0.PruneLen, c1.PruneLen, bot)
			if c0.PruneLen == 0 {
				fmt.Println(c0.BMLen, c0.PruneLen, c0.Nam)
			}
			nd.PruneLen = nd.BMLen + 1.0/bot
			//parSubtreePruneLen[nd] = 1.0 / bot
			for i := range sites {
				nd.ContData[i] = (((1 / c0.PruneLen) * c1.ContData[i]) + ((1 / c1.PruneLen) * c0.ContData[i])) / bot
			}
			virtualTritomyML(nd, sites)
			//nd.PruneLen = nd.BMLen + (1.0 / bot)
		}
		BMPruneRooted(c)
	}
	if len(t.Rt.Chs) == 3 {
		AncTritomyML(t.Rt, sites)
	} else if len(t.Rt.Chs) == 2 { // average between the estimated branch lengths descending from the root if the tree is rooted
		sumChsLens := (t.Rt.Chs[0].BMLen / 2) + (t.Rt.Chs[1].BMLen / 2)
		for _, c := range t.Rt.Chs {
			c.BMLen = sumChsLens / 2.0
		}
	}
}

//virtualTritomyML will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree
func virtualTritomyML(tree *Node, sites []int) {
	ntraits := len(sites)
	fntraits := float64(ntraits)
	var x1, x2, x3 float64
	sumV1 := 0.0
	sumV2 := 0.0
	sumV3 := 0.0
	for i := range sites {
		x1 = tree.Chs[0].ContData[i]
		x2 = tree.Chs[1].ContData[i]
		x3 = tree.ContData[i]
		sumV1 += ((x1 - x2) * (x1 - x3))
		sumV2 += ((x2 - x1) * (x2 - x3))
		sumV3 += ((x3 - x1) * (x3 - x2))
	}
	if sumV1 < 0.0 {
		sumV1 = 0.000001
		sumV2 = 0.0
		sumV3 = 0.0
		for i := range sites {
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.ContData[i]
			sumV2 += (x1 - x2) * (x1 - x2)
			sumV3 += (x1 - x3) * (x1 - x3)
		}
	} else if sumV2 < 0.0 {
		sumV1 = 0.0
		sumV2 = 0.00001
		sumV3 = 0.0
		for i := range sites {
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.ContData[i]
			sumV1 += (x2 - x1) * (x2 - x1)
			sumV3 += (x2 - x3) * (x2 - x3)
		}
	} else if sumV3 < 0.0 {
		sumV1 = 0.0
		sumV2 = 0.0
		sumV3 = 0.0001
		for i := range sites {
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.ContData[i]
			sumV1 += (x3 - x1) * (x3 - x1)
			sumV2 += (x3 - x2) * (x3 - x2)
		}
	}
	pruneDiff := tree.Chs[0].PruneLen - tree.Chs[0].BMLen
	tree.Chs[0].BMLen = (sumV1 / fntraits) - (tree.Chs[0].PruneLen - tree.Chs[0].BMLen)
	tree.Chs[0].PruneLen = tree.Chs[0].BMLen + pruneDiff
	pruneDiff = tree.Chs[1].PruneLen - tree.Chs[1].BMLen
	tree.Chs[1].BMLen = (sumV2 / fntraits) - (tree.Chs[1].PruneLen - tree.Chs[1].BMLen)
	tree.Chs[1].PruneLen = tree.Chs[1].BMLen + pruneDiff
	pruneDiff = (tree.PruneLen - tree.BMLen)
	tree.BMLen = (sumV3 / fntraits) - pruneDiff
	tree.PruneLen = tree.BMLen + pruneDiff
	if tree.Chs[0].BMLen <= 0. {
		tree.Chs[0].BMLen = 0.0001
		tree.Chs[0].PruneLen = tree.Chs[0].BMLen + pruneDiff
	}
	if tree.Chs[1].BMLen <= 0. {
		tree.Chs[1].BMLen = 0.0001
		tree.Chs[1].PruneLen = tree.Chs[1].BMLen + pruneDiff
	}
	if tree.BMLen <= 0. {
		tree.BMLen = 0.0001
		tree.PruneLen = tree.BMLen + pruneDiff
	}
}

//AncTritomyML will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree assuming that direct ancestors may be in the tree
func AncTritomyML(tree *Node, sites []int) {
	ntraits := len(sites)
	fntraits := float64(ntraits)
	var x1, x2, x3 float64
	sumV1 := 0.0
	sumV2 := 0.0
	sumV3 := 0.0
	for i := range sites {
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
		for i := range sites {
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
		for i := range sites {
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
		for i := range sites {
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
	if sumV1 <= 0. {
		sumV1 = 0.0
	}
	if sumV2 <= 0. {
		sumV2 = 0.0
	}
	if sumV3 <= 0. {
		sumV3 = 0.0
	}
	tree.Chs[0].BMLen = sumV1
	tree.Chs[1].BMLen = sumV2
	tree.Chs[2].BMLen = sumV3
	if math.IsNaN(sumV1) || math.IsNaN(sumV2) || math.IsNaN(sumV3) {
		fmt.Println("Error in AncTritomyML() in bmoptim.go")
		fmt.Println(tree.Chs[0].Nam, sumV1, tree.Chs[1].Nam, sumV2, tree.Chs[2].Nam, sumV3)
		os.Exit(0)
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
////These next functions are for calculating branch lengths when doing branch length modularity clustering
//////////////////////////////////////////////////////////////////////////////////////////////

//IterateLengthsWeighted will iteratively calculate the ML branch lengths for a particular topology and cluster when doing the greedy site clustering procedure.
func IterateLengthsWeighted(tree *Tree, cluster *Cluster, niter int) {
	AssertUnrootedTree(tree.Rt)
	//nodes := tree.PreorderArray()
	//InitMissingValues(nodes)
	for i := 0; i < niter; i++ {
		CalcExpectedTraits(tree.Rt)                //calculate Expected trait values
		calcBMLensBackFrontWeighted(tree, cluster) //maximize likelihood of branch lengths
	}
	var newlen []float64
	for _, n := range tree.Pre {
		newlen = append(newlen, n.BMLen) //store the newly calculated branch lengths
	}
	cluster.BranchLengths = newlen
}

//BMCalcLensBackFront will do one pass of the EM branch length estimation
func calcBMLensBackFrontWeighted(t *Tree, cluster *Cluster) {
	if len(t.Rt.Chs) == 3 { // the tree is unrooted
		//parSubtreePruneLen := make(map[*Node]float64)
		for _, c := range t.Rt.Chs {
			BMPruneRooted(c)
		}
		tritomyWeightedML(t.Rt, cluster.SiteWeights)
		for _, c := range t.Rt.Chs {
			BMPruneRooted(c)
		}
	} else {
		BMPruneRooted(t.Rt)
	}
	for _, c := range t.Rt.Chs {
		for _, nd := range c.PreorderArray() {
			if len(nd.Chs) == 0 { //don't do this at tips
				continue
			}
			var c0, c1 *Node
			var bot float64
			if nd.Par == t.Rt && len(t.Rt.Chs) == 3 {
				var sibs []*Node
				for _, c := range t.Rt.Chs {
					if c != nd {
						sibs = append(sibs, c)
					}
				}
				c0 = sibs[0]
				c1 = sibs[1]
			} else if nd.Par == t.Rt && len(t.Rt.Chs) == 2 {
				sib := nd.GetSib()
				if len(sib.Chs) == 0 { // special case if the tree is rooted on a single taxon
					nd.BMLen += sib.BMLen
					nd.ContData = sib.ContData
					nd.PruneLen = nd.BMLen
					virtualTritomyMLWeights(nd, cluster.SiteWeights)
					fmt.Println("WARNING: BM Branch length optimization is currently wonky on trees rooted on a single taxon")
					//os.Exit(0)
					nd.PruneLen = nd.BMLen
					sib.BMLen = nd.BMLen / 2
					nd.BMLen = nd.BMLen / 2
					continue // don't do the stuff below; move on to the next node
				} else {
					c0 = sib.Chs[0]
					c1 = sib.Chs[1]
					nd.BMLen += sib.BMLen
				}
			} else if nd.Par != t.Rt {
				c0 = nd.GetSib()
				c1 = nd.Par
			}
			bot = ((1.0 / c0.PruneLen) + (1.0 / c1.PruneLen))
			nd.PruneLen = nd.BMLen + 1.0/bot
			//parSubtreePruneLen[nd] = 1.0 / bot
			for i := range nd.ContData {
				nd.ContData[i] = (((1 / c0.PruneLen) * c1.ContData[i]) + ((1 / c1.PruneLen) * c0.ContData[i])) / bot
			}
			virtualTritomyMLWeights(nd, cluster.SiteWeights)
			nd.PruneLen = nd.BMLen + 1.0/bot
		}
		BMPruneRooted(c)
	}
	if len(t.Rt.Chs) == 3 {
		tritomyWeightedML(t.Rt, cluster.SiteWeights)
	} else if len(t.Rt.Chs) == 2 { // average between the estimated branch lengths descending from the root if the tree is rooted
		sumChsLens := (t.Rt.Chs[0].BMLen / 2) + (t.Rt.Chs[1].BMLen / 2)
		for _, c := range t.Rt.Chs {
			c.BMLen = sumChsLens / 2.0
		}
	}
}

//virtualTritomyMLWeights will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree
func virtualTritomyMLWeights(tree *Node, weights map[int]float64) {
	ntraits := len(tree.Chs[0].ContData)
	var x1, x2, x3, wt float64
	sumV1 := 0.0
	sumV2 := 0.0
	sumV3 := 0.0
	for i := 0; i < ntraits; i++ {
		wt = weights[i]
		if wt == 0.0 {
			continue
		}
		x1 = tree.Chs[0].ContData[i]
		x2 = tree.Chs[1].ContData[i]
		x3 = tree.ContData[i]
		sumV1 += ((x1 - x2) * (x1 - x3)) * wt
		sumV2 += ((x2 - x1) * (x2 - x3)) * wt
		sumV3 += ((x3 - x1) * (x3 - x2)) * wt
	}
	if sumV1 < 0.0 {
		sumV1 = 0.000001
		sumV2 = 0.0
		sumV3 = 0.0
		for i := 0; i < ntraits; i++ {
			wt = weights[i]
			if wt == 0.0 {
				continue
			}
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.ContData[i]
			sumV2 += ((x1 - x2) * (x1 - x2)) * wt
			sumV3 += ((x1 - x3) * (x1 - x3)) * wt
		}
	} else if sumV2 < 0.0 {
		sumV1 = 0.0
		sumV2 = 0.00001
		sumV3 = 0.0
		for i := 0; i < ntraits; i++ {
			wt = weights[i]
			if wt == 0.0 {
				continue
			}
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.ContData[i]
			sumV1 += ((x2 - x1) * (x2 - x1)) * wt
			sumV3 += ((x2 - x3) * (x2 - x3)) * wt
		}
	} else if sumV3 < 0.0 {
		sumV1 = 0.0
		sumV2 = 0.0
		sumV3 = 0.0001
		for i := 0; i < ntraits; i++ {
			wt = weights[i]
			if wt == 0.0 {
				continue
			}
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.ContData[i]
			sumV1 += ((x3 - x1) * (x3 - x1)) * wt
			sumV2 += ((x3 - x2) * (x3 - x2)) * wt
		}
	}
	tree.Chs[0].BMLen = (sumV1) - (tree.Chs[0].PruneLen - tree.Chs[0].BMLen)
	tree.Chs[1].BMLen = (sumV2) - (tree.Chs[1].PruneLen - tree.Chs[1].BMLen)
	parPruneDiff := (tree.PruneLen - tree.BMLen)
	tree.BMLen = (sumV3) - parPruneDiff
	tree.PruneLen = tree.BMLen + parPruneDiff

	if tree.Chs[0].BMLen <= 0. {
		tree.Chs[0].BMLen = 0.0001
	}
	if tree.Chs[1].BMLen <= 0. {
		tree.Chs[1].BMLen = 0.0001
	}
	if tree.BMLen <= 0. {
		tree.BMLen = 0.0001
	}
	if math.IsNaN(sumV1) || math.IsNaN(sumV2) || math.IsNaN(sumV3) {
		fmt.Println("error in virtualTritomyMLWeights() in bmoptim.go")
		fmt.Println(tree.Nam, sumV1, tree.Chs[0].Nam, sumV2, tree.Chs[1].Nam, sumV3)
		os.Exit(0)
	}

}

//tritomyWeightedML will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree
func tritomyWeightedML(tree *Node, weights map[int]float64) {
	//ntraits := len(tree.Chs[0].ContData)
	/*fntraits := 0.
	for _, w := range weights {
		fntraits += w
	}
	//fmt.Println(fntraits)*/
	var x1, x2, x3 float64
	sumV1 := 0.0
	sumV2 := 0.0
	sumV3 := 0.0
	wt := 0.0
	for i := range tree.Chs[0].ContData {
		wt = weights[i]
		if wt == 0.0 {
			continue
		}
		x1 = tree.Chs[0].ContData[i]
		x2 = tree.Chs[1].ContData[i]
		x3 = tree.Chs[2].ContData[i]

		sumV1 += ((x1 - x2) * (x1 - x3)) * wt
		sumV2 += ((x2 - x1) * (x2 - x3)) * wt
		sumV3 += ((x3 - x1) * (x3 - x2)) * wt
	}
	if sumV1 <= 0.0 {
		sumV1 = 0.000001
		sumV2 = 0.0
		sumV3 = 0.0
		for i := range tree.Chs[0].ContData {
			wt = weights[i]
			if wt == 0.0 {
				continue
			}
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.Chs[2].ContData[i]
			sumV2 += (x1 - x2) * (x1 - x2) * wt
			sumV3 += (x1 - x3) * (x1 - x3) * wt
		}
	} else if sumV2 <= 0.0 {
		sumV1 = 0.0
		sumV2 = 0.00001
		sumV3 = 0.0
		for i := range tree.Chs[0].ContData {
			wt = weights[i]
			if wt == 0.0 {
				continue
			}

			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.Chs[2].ContData[i]
			sumV1 += (x2 - x1) * (x2 - x1) * wt
			sumV3 += (x2 - x3) * (x2 - x3) * wt
		}
	} else if sumV3 <= 0.0 {
		sumV1 = 0.0
		sumV2 = 0.0
		sumV3 = 0.0001
		for i := range tree.Chs[0].ContData {
			wt = weights[i]
			if wt == 0.0 {
				continue
			}

			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.Chs[2].ContData[i]
			sumV1 += (x3 - x1) * (x3 - x1) * wt
			sumV2 += (x3 - x2) * (x3 - x2) * wt
		}
	}
	//sumV1 = sumV1 // fntraits
	//sumV2 = sumV2 // fntraits
	//sumV3 = sumV3 // fntraits
	sumV1 = sumV1 - (tree.Chs[0].PruneLen - tree.Chs[0].BMLen)
	sumV2 = sumV2 - (tree.Chs[1].PruneLen - tree.Chs[1].BMLen)
	sumV3 = sumV3 - (tree.Chs[2].PruneLen - tree.Chs[2].BMLen)
	if sumV1 <= 0. {
		sumV1 = 0.0001
	}
	if sumV2 <= 0. {
		sumV2 = 0.0001
	}
	if sumV3 <= 0. {
		sumV3 = 0.0001
	}
	//fmt.Println(tree.Chs[0].Nam, sumV1, tree.Chs[1].Nam, sumV2, tree.Chs[2].Nam, sumV3)
	if math.IsNaN(sumV1) || math.IsNaN(sumV2) || math.IsNaN(sumV3) {
		fmt.Println("error in virtualTritomyML() in bmoptim.go")
		fmt.Println(tree.Chs[0].Nam, sumV1, tree.Chs[1].Nam, sumV2, tree.Chs[2].Nam, sumV3)
		os.Exit(0)
	}
	tree.Chs[0].BMLen = sumV1
	tree.Chs[1].BMLen = sumV2
	tree.Chs[2].BMLen = sumV3
}
