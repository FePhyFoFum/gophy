package gophy

import (
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)

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
		if node.FINDS > 1 {
			f := node.FAD * c
			l := node.LAD * c
			a := math.Log(math.Pow(f-l, (node.FINDS - 2.0)))
			b := math.Log(math.Pow(lam, node.FINDS))
			c := -lam * (tf - tl)
			brlik = (a + b + c) - float64(LogFactorial(int(node.FINDS-2.0)))
		} else if node.FINDS == 1 {
			brlik = math.Log(lam) + (-lam * (tf - tl))
		} else if node.FINDS == 0 {
			brlik = -lam * (tf - tl)
		}
		treelik += brlik
	}
	return treelik
}

func ADPoissonTreeLoglike(nodels []*Node, lam float64) float64 {
	//lam := 10.0
	c := 1.0
	treelik := 0.0
	for _, node := range nodels {
		if node.Nam == "root" {
			continue
		} else if node.ANC == true && node.ISTIP == false { // will calc likelihood of ancestors on the corresponding tip
			continue
		}
		var tf float64
		if node.ANC == false {
			tf = node.Par.Height * c
		} else if node.ANC == true && node.ISTIP == true {
			tf = node.Par.Par.Height * c
		}
		tl := node.Height * c
		brlik := 0.0
		if node.FINDS > 1 {
			f := node.FAD * c
			l := node.LAD * c
			a := math.Log(math.Pow(f-l, (node.FINDS - 2.0)))
			b := math.Log(math.Pow(lam, node.FINDS))
			c := -lam * (tf - tl)
			brlik = (a + b + c) - float64(LogFactorial(int(node.FINDS-2.0)))
		} else if node.FINDS == 1 {
			brlik = math.Log(lam) + (-lam * (tf - tl))
		} else if node.FINDS == 0 {
			brlik = -lam * (tf - tl)
		}
		treelik += brlik
	}
	return treelik
}

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
					node.FAD = fad
					node.LAD = lad
					node.FINDS = n
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
	if node.ISTIP {
		if node.LAD == 0.0 {
			node.Height = node.LAD
		} else {
			node.Height = node.LAD - 0.000001
		}
	} else {
		oldestChildHeight := OldestChildAge(node)
		node.Height = oldestChildHeight + 0.001
		node.FAD = node.Height
	}

}

func MakeStratHeights(tree *Tree) {
	assignHeights(tree.Rt)
	postTree := tree.Post
	for _, node := range postTree {
		if node.Par == nil {
			continue
		} else {
			node.TimeLen = node.Par.Height - node.Height
		}
	}
}

func OldestChildAge(node *Node) float64 {
	oldestChildHeight := 0.0
	for _, c := range node.Chs {
		if c.FAD > oldestChildHeight && c.Nam+"_ancestral" != node.Nam {
			oldestChildHeight = c.FAD
		}
	}
	return oldestChildHeight
}
