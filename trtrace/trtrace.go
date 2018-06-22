package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"math/rand"
	"os"
	"strings"
	"time"

	"github.com/FePhyFoFum/gophy"
)

func main() {
	nt := flag.String("t", "", "treefile name")
	nl := flag.String("n", "", "list of ids")
	nf := flag.String("f", "", "filename with ids")
	flag.Parse()
	if len(os.Args) < 2 {
		flag.PrintDefaults()
		os.Exit(1)
	}
	var ids []string
	if len(*nl) > 0 {
		ids = strings.Split(*nl, ",")
		fmt.Fprint(os.Stderr, ids, "\n")
	} else if len(*nf) > 0 {
		ids = make([]string, 0)
		f, err := os.Open(*nf)
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		defer f.Close()
		scanner := bufio.NewReader(f)
		for {
			ln, err := scanner.ReadString('\n')
			if len(ln) > 0 {
				ids = append(ids, strings.Trim(ln, "\n"))
			}
			if err == io.EOF {
				break
			}
			if err != nil {
				log.Printf("read %d bytes: %v", ln, err)
				break
			}
		}
	}
	fmt.Fprint(os.Stderr, "reading tree\n")
	f, err := os.Open(*nt)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	defer f.Close()
	scanner := bufio.NewReader(f)
	var tree gophy.Tree
	for {
		ln, err := scanner.ReadString('\n')
		if len(ln) > 0 {
			rt := gophy.ReadNewickString(ln)
			tree.Instantiate(rt)
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Printf("read %d bytes: %v", ln, err)
			break
		}
	}
	namedNodeMap := make(map[string]*gophy.Node)
	for _, n := range tree.Pre {
		if len(n.Nam) < 1 {
			continue
		}
		namedNodeMap[n.Nam] = n
	}
	nds := make([]*gophy.Node, 0)
	unmatched := make([]string, 0)
	for _, i := range ids {
		if _, ok := namedNodeMap[i]; ok {
			nds = append(nds, namedNodeMap[i])
		} else {
			unmatched = append(unmatched, i)
		}
	}
	//tracing
	r := rand.New(rand.NewSource(time.Now().UnixNano()))
	rid := r.Float64()
	traceTree(nds, rid)
	addLen := make(map[string]float64)
	for _, i := range nds {
		if len(i.Chs) > 0 {
			if len(getMarkedChs(i, rid)) == 0 {
				addLen[i.Nam] = getSubTendLen(i)
			}
		}
	}
	var newick string
	if len(unmatched) > 0 {
		fmt.Fprint(os.Stderr, "unmatched: ", unmatched, "\n")
	}
	curRt := tree.Rt
	going := true
	for going {
		x := getMarkedChs(curRt, rid)
		if len(x) == 1 {
			curRt = x[0]
		} else {
			going = false
			break
		}
	}
	x := curRt.NewickPaint(true, rid) + ";"
	//handle the tips thjat are not tips
	//this takes time but necessary
	nwt := gophy.ReadNewickString(x)
	var ntt gophy.Tree
	ntt.Instantiate(nwt)
	for _, i := range ntt.Tips {
		if _, ok := addLen[i.Nam]; ok {
			i.Len += addLen[i.Nam]
		}
	}
	newick = nwt.Newick(true) + ";"
	//end the handle
	untraceTree(nds, rid)
	fmt.Fprintf(os.Stdout, newick)
}

func getSubTendLen(nd *gophy.Node) float64 {
	x := 0.0
	going := true
	curNd := nd.Chs[0]
	for going {
		x += curNd.Len
		if len(curNd.Chs) > 0 {
			curNd = curNd.Chs[0]
		} else {
			going = false
			break
		}
	}
	return x
}

func getMarkedChs(nd *gophy.Node, rid float64) []*gophy.Node {
	x := make([]*gophy.Node, 0)
	for _, i := range nd.Chs {
		if _, ok := i.MarkedMap[rid]; ok {
			x = append(x, i)
		}
	}
	return x
}

func traceTree(nds []*gophy.Node, rid float64) {
	for _, n := range nds {
		n.MarkedMap[rid] = true
		going := true
		cur := n.Par
		for going {
			if _, ok := cur.MarkedMap[rid]; ok {
				break
			}
			cur.MarkedMap[rid] = true
			if cur.Par == nil {
				break
			}
			cur = cur.Par
		}
	}
}

func untraceTree(nds []*gophy.Node, rid float64) {
	for _, n := range nds {
		if _, ok := n.MarkedMap[rid]; ok {
			delete(n.MarkedMap, rid)
		}
		going := true
		cur := n.Par
		for going {
			if _, ok := cur.MarkedMap[rid]; !ok {
				break
			}
			delete(cur.MarkedMap, rid)
			if cur.Par == nil {
				break
			}
			cur = cur.Par
		}
	}
}
