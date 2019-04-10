package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/FePhyFoFum/gophy"
)

func MakeBDTree(showDead bool) (root *gophy.Node) {
	rt := gophy.Node{}

	root = &rt
	return
}

func main() {
	if len(os.Args) < 2 {
		fmt.Fprint(os.Stderr, "bdsim -e numtips\n")
		os.Exit(1)
	}
	nt := flag.Int("e", 3, "how many extant tips")
	flag.Parse()
	fmt.Fprintln(os.Stderr, *nt)
	MakeBDTree(false)
}
