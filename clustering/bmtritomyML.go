package clustering

import (
	"fmt"
	"math"
	"os"

	"github.com/FePhyFoFum/gophy"
)

//TritomySubML will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree using only the sites indicated in sites
func TritomySubML(tree *gophy.Node, sites []int) {
	ntraits := len(sites)
	fntraits := float64(ntraits)
	var x1, x2, x3 float64
	sumV1 := 0.0
	sumV2 := 0.0
	sumV3 := 0.0
	for _, i := range sites {
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
		for _, i := range sites {
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
		for _, i := range sites {
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
		for _, i := range sites {
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

	if sumV1 <= 0. {
		sumV1 = 0.0001
	}
	if sumV2 <= 0. {
		sumV2 = 0.0001
	}
	if sumV3 <= 0. {
		sumV3 = 0.0001
	}
	tree.Chs[0].Len = sumV1
	tree.Chs[1].Len = sumV2
	tree.Chs[2].Len = sumV3
}

//TritomyWeightedML will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree
func TritomyWeightedML(tree *gophy.Node, weights map[int]float64) {
	//ntraits := len(tree.Chs[0].ContData)
	fntraits := 0.
	for _, w := range weights {
		fntraits += w
	}
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
	sumV1 = sumV1 / fntraits
	sumV2 = sumV2 / fntraits
	sumV3 = sumV3 / fntraits
	sumV1 = sumV1 - (tree.Chs[0].PruneLen - tree.Chs[0].Len)
	sumV2 = sumV2 - (tree.Chs[1].PruneLen - tree.Chs[1].Len)
	sumV3 = sumV3 - (tree.Chs[2].PruneLen - tree.Chs[2].Len)
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
		fmt.Println(tree.Chs[0].Nam, sumV1, tree.Chs[1].Nam, sumV2, tree.Chs[2].Nam, sumV3)
		os.Exit(0)
	}
	tree.Chs[0].Len = sumV1
	tree.Chs[1].Len = sumV2
	tree.Chs[2].Len = sumV3
}
