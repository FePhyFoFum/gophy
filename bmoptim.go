package gophy

import "fmt"

//BMOptimBLEM will calculate the BM branch lengths using an iterative EM calculation that imputes missing data using PICs
func BMOptimBLEM(t *Tree, niter int) {
	CalcExpectedTraits(t.Rt)
	for i := 0; i < niter; i++ {
		BMCalcLensBackFront(t)
		CalcExpectedTraits(t.Rt)
	}
	BMCalcLensBackFront(t)
}

//BMCalcLensBackFront will do one pass of the EM branch length estimation
func BMCalcLensBackFront(t *Tree) {
	if len(t.Rt.Chs) == 3 { // the tree is unrooted
		//parSubtreePruneLen := make(map[*Node]float64)
		for _, c := range t.Rt.Chs {
			BMPruneRooted(c)
		}
		AncTritomyML(t.Rt)
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
					virtualTritomyML(nd)
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
			virtualTritomyML(nd)
			nd.PruneLen = nd.BMLen + 1.0/bot
		}
		BMPruneRooted(c)
	}
	if len(t.Rt.Chs) == 3 {
		AncTritomyML(t.Rt)
	} else if len(t.Rt.Chs) == 2 { // average between the estimated branch lengths descending from the root if the tree is rooted
		sumChsLens := (t.Rt.Chs[0].BMLen / 2) + (t.Rt.Chs[1].BMLen / 2)
		for _, c := range t.Rt.Chs {
			c.BMLen = sumChsLens / 2.0
		}
	}
}

//virtualTritomyML will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree
func virtualTritomyML(tree *Node) {
	ntraits := len(tree.Chs[0].ContData)
	fntraits := float64(ntraits)
	var x1, x2, x3 float64
	sumV1 := 0.0
	sumV2 := 0.0
	sumV3 := 0.0
	for i := 0; i < ntraits; i++ {
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
		for i := 0; i < ntraits; i++ {
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
		for i := 0; i < ntraits; i++ {
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
		for i := 0; i < ntraits; i++ {
			x1 = tree.Chs[0].ContData[i]
			x2 = tree.Chs[1].ContData[i]
			x3 = tree.ContData[i]
			sumV1 += (x3 - x1) * (x3 - x1)
			sumV2 += (x3 - x2) * (x3 - x2)
		}
	}
	tree.Chs[0].BMLen = (sumV1 / fntraits) - (tree.Chs[0].PruneLen - tree.Chs[0].BMLen)
	tree.Chs[1].BMLen = (sumV2 / fntraits) - (tree.Chs[1].PruneLen - tree.Chs[1].BMLen)
	tree.BMLen = (sumV3 / fntraits) - (tree.PruneLen - tree.BMLen)
	if tree.Chs[0].BMLen < 0. {
		tree.Chs[0].BMLen = 0.0001
	}
	if tree.Chs[1].BMLen < 0. {
		tree.Chs[1].BMLen = 0.0001
	}
	if tree.BMLen < 0. {
		tree.BMLen = 0.0001
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
