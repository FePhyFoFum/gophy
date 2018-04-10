package main

import (
	"bufio"
	"flag"
	"fmt"
	"gophy"
	"io"
	"os"
	"strconv"
	"strings"
)

//this is for scaling a tree to be ultrametric
// just splitting the difference and all that

//mrca filename should be
//name1,name2 date

func main() {
	tfn := flag.String("t", "", "tree filename")
	mfn := flag.String("m", "", "mrca filename")
	flag.Parse()
	if len(*tfn) == 0 {
		os.Exit(0)
	}
	if len(*mfn) == 0 {
		os.Exit(0)
	}
	fmt.Fprintln(os.Stderr, "treefile:", *tfn)
	fmt.Fprintln(os.Stderr, "mrcafile:", *mfn)

	// read tree file
	f, err := os.Open(*tfn)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	defer f.Close()
	scanner := bufio.NewReader(f)
	fmt.Fprintln(os.Stderr, "reading trees")
	var t gophy.Tree
	var rt *gophy.Node
	nmsnds := make(map[string]*gophy.Node)
	for {
		ln, err := scanner.ReadString('\n')
		if len(ln) > 0 {
			rt = gophy.ReadNewickString(ln)
			t.Instantiate(rt)
			for _, i := range t.Post {
				if len(i.Nam) > 0 {
					nmsnds[i.Nam] = i
				}
			}
			break
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			break
		}
	}
	//end tree file reading
	//read mrca file
	f, err = os.Open(*mfn)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	defer f.Close()
	mrcas := make(map[*gophy.Node]float64) //map is node and float is date
	scanner = bufio.NewReader(f)
	fmt.Fprintln(os.Stderr, "reading mrcas")
	for {
		ln, err := scanner.ReadString('\n')
		if len(ln) > 0 {
			fmt.Fprintln(os.Stderr, strings.Trim(ln, "\n"))
			spls1 := strings.Split(strings.Trim(ln, "\n"), " ") //mrca = 0, date = 1
			spls2 := strings.Split(spls1[0], ",")
			nds := make([]*gophy.Node, 2)
			nds[0] = nmsnds[spls2[0]]
			nds[1] = nmsnds[spls2[1]]
			ff, ferr := strconv.ParseFloat(spls1[1], 64)
			if ferr != nil {
				fmt.Fprintln(os.Stderr, "problem parsing", spls1[1], "as float64")
			}
			mrcas[gophy.GetMrca(nds, rt)] = ff
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			break
		}
	}

	// scale the tree
	for _, i := range t.Pre {
		if _, ok := mrcas[i]; ok {
			if i != rt {
				ih := i.Height
				i.Len = i.Len + ih - mrcas[i]
			}
			scaleSubTree(i, mrcas[i])
			gophy.SetHeights(&t)
		}
	}
	fmt.Println(t.Rt.Newick(true) + ";")
}

func getMaxNodes(node *gophy.Node) float64 {
	v := 0.
	if len(node.Chs) == 0 {
		return 0.
	}
	for _, i := range node.GetTips() {
		vl := 1.
		going := true
		cur := i
		for going {
			par := cur.Par
			if par == node {
				going = false
				break
			} else {
				vl += 1.
				cur = par
			}
		}
		if vl > v {
			v = vl
		}
	}
	return v
}

func scaleSubTree(node *gophy.Node, value float64) {
	for _, i := range node.Chs {
		mx := getMaxNodes(i)
		cv := mx + 1.
		v := value / cv
		i.Len = v
		scaleSubTree(i, value-v)
	}
}
