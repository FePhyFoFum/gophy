package clustering

import (
	"fmt"
	"math"
	"os"

	"github.com/FePhyFoFum/gophy"
)

//InitMissingValues will find the missing sites in a data matrix and plug in values corresponding to the mean of the remaining sites
func InitMissingValues(tree []*gophy.Node) {
	means := CalcSiteMeans(tree)
	for _, n := range tree {
		if len(n.Chs) == 0 {
			MakeMissingMeansTip(n, means)
		}
	}
}

//MakeMissingMeansTip will replace missing values with the mean across all tips for a single tip
func MakeMissingMeansTip(n *gophy.Node, means []float64) {
	for i := range n.ContData {
		if n.Mis[i] == true {
			n.ContData[i] = means[i]
		}
	}
}

//CalcSiteMeans will calculate the mean value for all the sites in the matrix for which the site is not missing
func CalcSiteMeans(nodes []*gophy.Node) (siteSum []float64) {
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
func CalcExpectedTraits(tree *gophy.Node) {
	if len(tree.Chs) == 3 {
		for _, c := range tree.Chs {
			BMPruneRooted(c)
		}
	} else if len(tree.Chs) == 2 {
		BMPruneRooted(tree)
	}
	for _, c := range tree.Chs {
		for _, nd := range c.PreorderArray() {
			if nd.Par == tree {
				if len(tree.Chs) == 3 {
					var sibs []*gophy.Node
					for _, c := range tree.Chs {
						if c != nd {
							sibs = append(sibs, c)
						}
					}
					bot := ((1.0 / sibs[0].PruneLen) + (1.0 / sibs[1].PruneLen))
					nd.PruneLen = nd.BMLen + (1.0 / bot)
					for ind := range nd.ContData {
						if nd.Mis[ind] == true || len(nd.Chs) == 2 {
							nd.ContData[ind] = (((1 / sibs[0].PruneLen) * sibs[1].ContData[ind]) + ((1 / sibs[1].PruneLen) * sibs[0].ContData[ind])) / bot
						} else if len(nd.Chs) > 2 {
							fmt.Println("This tree has a polytomy somewhere other than the root. Please fix.")
							os.Exit(0)
						}
					}
				} else if len(tree.Chs) == 2 {
					if len(nd.GetSib().Chs) > 0 {
						sib := nd.GetSib()
						sibs := sib.Chs //if the tree is rooted and you are imputing data at a taxon descending from the root, calc average from sibling's children (to get other subtree)
						bot := ((1.0 / sibs[0].PruneLen) + (1.0 / sibs[1].PruneLen))
						//nd.PruneLen = (sib.BMLen + nd.BMLen) + (1.0 / bot)
						nd.PruneLen = nd.BMLen + sib.PruneLen
						for ind := range nd.ContData {
							nd.ContData[ind] = (((1 / sibs[0].PruneLen) * sibs[1].ContData[ind]) + ((1 / sibs[1].PruneLen) * sibs[0].ContData[ind])) / bot
						}
					} else { //if the subtree on the other side of the root only has a single tip
						sib := nd.GetSib()
						nd.PruneLen += sib.PruneLen
						nd.ContData = sib.ContData
					}
				}
			} else {
				sib := nd.GetSib()
				for ind := range nd.ContData {
					if len(nd.Chs) > 2 {
						fmt.Println("This tree has a polytomy somewhere other than the root. Please fix.")
						os.Exit(0)
					}
					if nd.Mis[ind] == true || len(nd.Chs) != 0 {
						bot := ((1.0 / sib.PruneLen) + (1.0 / nd.Par.PruneLen))
						//fmt.Println(sib.Nam, sib.PruneLen, nd.Par.Nam, nd.Par.PruneLen)

						if math.IsNaN(bot) || math.IsInf(bot, 0) {
							fmt.Println(sib.Nam, sib.PruneLen, nd.Par.Nam, nd.Par.PruneLen)
							os.Exit(0)
						}

						trt := (((1 / sib.PruneLen) * nd.Par.ContData[ind]) + ((1 / nd.Par.PruneLen) * sib.ContData[ind])) / bot
						if math.IsNaN(trt) {
							fmt.Println(sib.Nam, sib.ContData[ind], nd.Par.Nam, nd.Par.ContData[ind])
							os.Exit(0)
						}
						nd.ContData[ind] = trt
						nd.PruneLen = nd.BMLen + 1.0/bot
					}
				}
			}
		}
		BMPruneRooted(c)
	}
}

//ClusterCalcExpectedTraits will plug in the expected values for missing traits under BM using the pruning/PIC ancestral state estimation approach
func ClusterCalcExpectedTraits(tree *gophy.Node, sites []int) {
	if len(tree.Chs) == 3 {
		for _, c := range tree.Chs {
			BMPruneRooted(c)
		}
	} else if len(tree.Chs) == 2 {
		BMPruneRooted(tree)
	}
	for _, c := range tree.Chs {
		for _, nd := range c.PreorderArray() {
			if nd.Par == tree {
				if len(tree.Chs) == 3 {
					var sibs []*gophy.Node
					for _, c := range tree.Chs {
						if c != nd {
							sibs = append(sibs, c)
						}
					}
					bot := ((1.0 / sibs[0].PruneLen) + (1.0 / sibs[1].PruneLen))
					nd.PruneLen = nd.BMLen + (1.0 / bot)
					for ind := range sites {
						if nd.Mis[ind] == true || len(nd.Chs) == 2 {
							nd.ContData[ind] = (((1 / sibs[0].PruneLen) * sibs[1].ContData[ind]) + ((1 / sibs[1].PruneLen) * sibs[0].ContData[ind])) / bot
						} else if len(nd.Chs) > 2 {
							fmt.Println("This tree has a polytomy somewhere other than the root. Please fix.")
							os.Exit(0)
						}
					}
				} else if len(tree.Chs) == 2 {
					if len(nd.GetSib().Chs) > 0 {
						sib := nd.GetSib()
						sibs := sib.Chs //if the tree is rooted and you are imputing data at a taxon descending from the root, calc average from sibling's children (to get other subtree)
						bot := ((1.0 / sibs[0].PruneLen) + (1.0 / sibs[1].PruneLen))
						//nd.PruneLen = (sib.BMLen + nd.BMLen) + (1.0 / bot)
						nd.PruneLen = nd.BMLen + sib.PruneLen
						for ind := range sites {
							nd.ContData[ind] = (((1 / sibs[0].PruneLen) * sibs[1].ContData[ind]) + ((1 / sibs[1].PruneLen) * sibs[0].ContData[ind])) / bot
						}
					} else { //if the subtree on the other side of the root only has a single tip
						sib := nd.GetSib()
						nd.PruneLen += sib.PruneLen
						nd.ContData = sib.ContData
					}
				}
			} else {
				sib := nd.GetSib()
				for ind := range sites {
					if len(nd.Chs) > 2 {
						fmt.Println("This tree has a polytomy somewhere other than the root. Please fix.")
						os.Exit(0)
					}
					if nd.Mis[ind] == true || len(nd.Chs) != 0 {
						bot := ((1.0 / sib.PruneLen) + (1.0 / nd.Par.PruneLen))
						//fmt.Println(sib.Nam, sib.PruneLen, nd.Par.Nam, nd.Par.PruneLen)

						if math.IsNaN(bot) || math.IsInf(bot, 0) {
							fmt.Println("HERE")
							fmt.Println(sib.Nam, sib.PruneLen, nd.Par.Nam, nd.Par.PruneLen)
							os.Exit(0)
						}

						trt := (((1 / sib.PruneLen) * nd.Par.ContData[ind]) + ((1 / nd.Par.PruneLen) * sib.ContData[ind])) / bot
						if math.IsNaN(trt) {
							fmt.Println(sib.Nam, sib.ContData[ind], nd.Par.Nam, nd.Par.ContData[ind])
							os.Exit(0)
						}
						nd.ContData[ind] = trt
						nd.PruneLen = nd.BMLen + 1.0/bot
					}
				}
			}
		}
		BMPruneRooted(c)
	}
}

func calcSingleExpTrait(c1, c2 *gophy.Node, ind int) (val float64) {
	bot := ((1.0 / c1.PruneLen) + (1.0 / c2.PruneLen))
	val = (((1 / c1.PruneLen) * c2.ContData[ind]) + ((1 / c2.PruneLen) * c1.ContData[ind])) / bot
	return
}
