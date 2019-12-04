package gophy

// Calculate stratigraphic likelihoods with trees

import (
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

// PoissonTreeLoglike calculates the Poisson LogLike based on
//  stratigraphic ranges
func PoissonTreeLoglike(tree *Tree) float64 {
	nodels := tree.Pre
	lam := 1.0
	c := 1.0
	treelik := 0.0
	for _, node := range nodels {
		if node.Nam == "root" {
			continue
		}
		tf := node.Par.Height * c
		tl := node.Height * c
		brlik := 0.0
		if node.FData["FINDS"] > 1 {
			f := node.FData["FAD"] * c
			l := node.FData["LAD"] * c
			a := math.Log(math.Pow(f-l, (node.FData["FINDS"] - 2.0)))
			b := math.Log(math.Pow(lam, node.FData["FINDS"]))
			c := -lam * (tf - tl)
			brlik = (a + b + c) - float64(LogFactorial(int(node.FData["FINDS"]-2.0)))
		} else if node.FData["FINDS"] == 1 {
			brlik = math.Log(lam) + (-lam * (tf - tl))
		} else if node.FData["FINDS"] == 0 {
			brlik = -lam * (tf - tl)
		}
		treelik += brlik
	}
	return treelik
}

// ADPoissonTreeLoglike calculates the ancestor - descendent
//  for a set of stratigraphic ranges
func ADPoissonTreeLoglike(nodels []*Node, lam float64) float64 {
	//lam := 10.0
	c := 1.0
	treelik := 0.0
	for _, node := range nodels {
		if node.Nam == "root" {
			continue
		} else if node.Anc == true && len(node.Chs) == 0 { // will calc likelihood of ancestors on the corresponding tip
			continue
		}
		var tf float64
		if node.Anc == false {
			tf = node.Par.Height * c
		} else if node.Anc == true && len(node.Chs) != 0 {
			tf = node.Par.Par.Height * c
		}
		tl := node.Height * c
		brlik := 0.0
		if node.FData["FINDS"] > 1 {
			f := node.FData["FAD"] * c
			l := node.FData["LAD"] * c
			a := math.Log(math.Pow(f-l, (node.FData["FINDS"] - 2.0)))
			b := math.Log(math.Pow(lam, node.FData["FINDS"]))
			c := -lam * (tf - tl)
			brlik = (a + b + c) - float64(LogFactorial(int(node.FData["FINDS"]-2.0)))
		} else if node.FData["FINDS"] == 1 {
			brlik = math.Log(lam) + (-lam * (tf - tl))
		} else if node.FData["FINDS"] == 0 {
			brlik = -lam * (tf - tl)
		}
		treelik += brlik
	}
	return treelik
}

// ReadStrat reads stratigraphic ranges from a file and assigns those
//  data to a tree
func ReadStrat(stratfl string, t *Tree) {
	nodels := t.Pre
	lines := ReadLine(stratfl)
	matchcount := 0
	linecount := 0
	matched := make(map[string]bool)
	for _, line := range lines {
		if line == "" {
			continue
		}
		ss := strings.Split(line, "\t")
		if len(ss) != 4 {
			fmt.Println("this line:")
			fmt.Println(line)
			fmt.Println("does not have 4 columns. should be <taxon\tFAD\tLAD\tn_occurrences>")
			os.Exit(0)
		}
		curtax := ss[0]
		var fad, lad, n float64
		var err error
		fad, err = strconv.ParseFloat(ss[1], 64)
		if err != nil {
			fmt.Println("couldn't convert FAD on this line:")
			fmt.Println(line)
			fmt.Println("to float")
			os.Exit(0)
		}
		lad, err = strconv.ParseFloat(ss[2], 64)
		if err != nil {
			fmt.Println("couldn't convert LAD on this line:")
			fmt.Println(line)
			fmt.Println("to float")
			os.Exit(0)
		}
		n, err = strconv.ParseFloat(ss[3], 64)
		if err != nil {
			fmt.Println("couldn't convert n_occurrences on this line:")
			fmt.Println(line)
			fmt.Println("to float")
			os.Exit(0)
		}
		linecount++
		for _, node := range nodels {
			if node.Nam == curtax {
				if _, ok := matched[curtax]; !ok {
					node.FData["FAD"] = fad
					node.FData["LAD"] = lad
					node.FData["FINDS"] = n
					matched[curtax] = true
					matchcount++
				} else {
					fmt.Println("range file contains duplicate entries for taxon ", curtax)
					os.Exit(0)
				}
			}
		}
	}
	if matchcount != linecount {
		fmt.Println("WARNING: the range file did not contain a range for all of the taxa in the tree")
		for _, node := range nodels {
			if node.Nam == "" {
				continue
			}
			if _, ok := matched[node.Nam]; !ok {
				fmt.Println(node.Nam)
			}
		}
		fmt.Println("is/are not in the provided range file")
	}
}

func assignHeights(node *Node) {
	for _, chld := range node.Chs {
		assignHeights(chld)
	}
	if len(node.Chs) == 0 {
		if node.FData["LAD"] == 0.0 {
			node.Height = node.FData["LAD"]
		} else {
			node.Height = node.FData["LAD"] - 0.000001
		}
	} else {
		oldestChildHeight := OldestChildAge(node)
		node.Height = oldestChildHeight + 0.1
		node.FData["FAD"] = node.Height
	}

}

// MakeStratHeights assigns the strat heights
func MakeStratHeights(tree *Tree) {
	assignHeights(tree.Rt)
	postTree := tree.Post
	for _, node := range postTree {
		if node.Par == nil {
			continue
		} else {
			node.FData["TimeLen"] = node.Par.Height - node.Height
		}
	}
}

//TimeTraverse will visit all descendant nodes in order of their heights (earliest -> latest)
func TimeTraverse(preNodes []*Node, internalOnly bool) (ret []*Node) {
	var unsortNodes []*Node
	if internalOnly == true {
		for _, n := range preNodes {
			if len(n.Chs) != 0 {
				unsortNodes = append(unsortNodes, n)
			}
		}
	} else {
		unsortNodes = preNodes
	}
	ss := unsortNodes
	sort.Slice(ss, func(i, j int) bool {
		return ss[i].Height > ss[j].Height
	})
	ret = ss
	return
}

/*	var orderedChs []*Node
	added := make(map[*Node]bool)
	for _, cn := range unsortNodes {
		if cn.Height > unsortNodes[0].Height {
			orderedChs = append(orderedChs, cn) //put all Chs with height >  the first element in the new array first
			added[cn] = true
		}
	}
	for _, cn := range
	ret = orderedChs
	return
}
*/
// OldestChildAge returns the oldest Child
func OldestChildAge(node *Node) float64 {
	oldestChildHeight := 0.0
	for _, c := range node.Chs {
		if c.FData["FAD"] > oldestChildHeight && c.Nam+"_ancestral" != node.Nam {
			oldestChildHeight = c.FData["FAD"]
		}
	}
	return oldestChildHeight
}
