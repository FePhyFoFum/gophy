package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"time"

	"github.com/FePhyFoFum/gophy"

	"golang.org/x/exp/rand"
)

func main() {
	rand.Seed(uint64(time.Now().UTC().UnixNano()))
	tfn := flag.String("t", "", "tree filename")
	flag.Parse()
	if len(*tfn) == 0 {
		fmt.Fprintln(os.Stderr, "need a filename")
		os.Exit(1)
	}

	//read a tree file
	f, err := os.Open(*tfn)
	if err != nil {
		fmt.Println(err)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	var rt *gophy.Node
	t := gophy.NewTree()
	for scanner.Scan() {
		ln := scanner.Text()
		if len(ln) < 2 {
			continue
		}
		rt = gophy.ReadNewickString(ln)
		t.Instantiate(rt)
		break // read the first tree
	}
	fmt.Println(t.Rt.Newick(false) + ";")

	gophy.TritomyRoot(t)
	fmt.Println(t.Rt.Newick(false) + ";")

	moves := gophy.NNIMoves(t)
	for i, j := range moves {
		gophy.SwapBranch(j[0], j[1])
		fmt.Println(i, t.Rt.Newick(false)+";")
		gophy.SwapBranch(j[0], j[1])
		//fmt.Println("original", t.Rt.Newick(false)+";")
	}
}
