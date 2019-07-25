package gophy

//InitMissingValues will find the missing sites in a data matrix and plug in values corresponding to the mean of the remaining sites
func InitMissingValues(tree []*Node) {
	means := CalcSiteMeans(tree)
	for _, n := range tree {
		if len(n.Chs) == 0 {
			MakeMissingMeansTip(n, means)
		}
	}
}

//MakeMissingMeansTip will replace missing values with the mean across all tips for a single tip
func MakeMissingMeansTip(n *Node, means []float64) {
	for i := range n.ContData {
		if n.Mis[i] == true {
			n.ContData[i] = means[i]
		}
	}
}

//CalcSiteMeans will calculate the mean value for all the sites in the matrix for which the site is not missing
func CalcSiteMeans(nodes []*Node) (siteSum []float64) {
	var ntraits []int
	for range nodes[0].ContData {
		siteSum = append(siteSum, 0.0)
		ntraits = append(ntraits, 0)
	}
	for _, n := range nodes {
		if len(n.Chs) != 0 {
			continue
		}
		for i, tr := range n.ContData {
			if n.Mis[i] != true {
				siteSum[i] += tr
				ntraits[i]++
			}
		}
	}
	for i := range siteSum {
		siteSum[i] = siteSum[i] / float64(ntraits[i])
	}
	return
}

//CalcExpectedTraits will plug in the expected values for missing traits under BM using the pruning/PIC ancestral state estimation approach
func CalcExpectedTraits(tree *Node) {
	rnodes := tree.PreorderArray()
	lnode := 0
	var expect float64
	var bot float64
	var top float64
	var ltrait float64
	var llen float64
	for ind, newroot := range rnodes { //visit each internal node and check any leaves for missing data. if found, calculate input expected value as the PIC at newroot
		if len(newroot.Chs) == 0 {
			continue
		} else if newroot != rnodes[0] {
			tree = newroot.Reroot(rnodes[lnode])
			lnode = ind
		}
		for _, cn := range tree.Chs {
			if len(cn.Chs) == 0 { // visit any leaves subtending from newroot
				for traitIndex := range cn.ContData {
					if cn.Mis[traitIndex] == true { // check if trait is missing and calculate expectation for each missing value
						bot = 0.
						childCount := 0
						for _, cn2 := range tree.Chs { //calculate PIC at newroot
							if cn2 != cn {
								childCount++
								BMPruneRootedSingle(cn2, traitIndex) // prune root to 3-tip tree
								bot += 1. / cn2.PruneLen
								if childCount == 2 {
									top = ((1. / cn2.PruneLen) * ltrait) + ((1. / llen) * cn2.ContData[traitIndex])
									break
								}
								ltrait = cn2.ContData[traitIndex]
								llen = cn2.PruneLen
							}
						}
						expect = top / bot
						cn.ContData[traitIndex] = expect
					}
				}
			}
		}
	}
	tree = rnodes[0].Reroot(tree)
	//fmt.Println(tree.Newick(true))
}

//AncCalcExpectedTraits will plug in the expected values for missing traits under BM using the pruning/PIC ancestral state estimation approach
func AncCalcExpectedTraits(tree *Node) {
	rnodes := tree.PreorderArray()
	lnode := 0
	var expect float64
	var bot float64
	var top float64
	var ltrait float64
	var llen float64
	for ind, newroot := range rnodes { //visit each internal node and check any leaves for missing data. if found, calculate input expected value as the PIC at newroot
		if len(newroot.Chs) == 0 {
			continue
		} else if newroot != rnodes[0] {
			tree = newroot.RerootBM(rnodes[lnode])
			lnode = ind
		}
		for _, cn := range tree.Chs {
			if len(cn.Chs) == 0 { // visit any leaves subtending from newroot
				for traitIndex := range cn.ContData {
					if cn.Mis[traitIndex] == true { // check if trait is missing and calculate expectation for each missing value
						bot = 0.
						childCount := 0
						for _, cn2 := range tree.Chs { //calculate PIC at newroot
							if cn2 != cn {
								childCount++
								AncBMPruneRootedSingle(cn2, traitIndex) // prune root to 3-tip tree
								bot += 1. / cn2.PruneLen
								if childCount == 2 {
									top = ((1. / cn2.PruneLen) * ltrait) + ((1. / llen) * cn2.ContData[traitIndex])
									break
								}
								ltrait = cn2.ContData[traitIndex]
								llen = cn2.PruneLen
							}
						}
						expect = top / bot
						cn.ContData[traitIndex] = expect
					}
				}
			}
		}
	}
	tree = rnodes[0].RerootBM(tree)
	//fmt.Println(tree.Newick(true))
}
