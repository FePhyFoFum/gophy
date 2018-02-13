package main

import (
	"flag"
	"fmt"
	"gophy"
	"math/rand"
	"os"
	"strconv"
	"time"
)

/*
 This just makes completely random trees
*/

func main() {
	if len(os.Args) < 2 {
		fmt.Fprint(os.Stderr, "rantree -n numtips\n")
		os.Exit(1)
	}
	nt := flag.Int("nt", 3, "how many tips")
	flag.Parse()
	nodeints := []int{}
	nodemap := make(map[int]*gophy.Node)
	for i := 0; i < *nt; i++ {
		nodeints = append(nodeints, i)
		nd := new(gophy.Node)
		nd.Nam = "nd_" + strconv.Itoa(i)
		nodemap[i] = nd
	}

	start := time.Now()
	curnodenum := len(nodemap)
	rand.Seed(time.Now().Unix())
	for len(nodeints) > 1 {
		n1 := rand.Int() % len(nodeints)
		tnd1 := nodemap[nodeints[n1]]
		nodeints = append(nodeints[:n1], nodeints[n1+1:]...)
		n1 = rand.Int() % len(nodeints)
		tnd2 := nodemap[nodeints[n1]]
		nodeints = append(nodeints[:n1], nodeints[n1+1:]...)
		nd := new(gophy.Node)
		nd.Nam = "nd_" + strconv.Itoa(curnodenum)
		tnd1.Par = nd
		tnd2.Par = nd
		nd.Chs = append(nd.Chs, tnd1)
		nd.Chs = append(nd.Chs, tnd2)
		nodeints = append(nodeints, curnodenum)
		nodemap[curnodenum] = nd
		curnodenum++
	}
	fmt.Println(nodemap[nodeints[0]].Newick(false) + ";")
	end := time.Now()
	fmt.Fprintln(os.Stderr, end.Sub(start))
}
