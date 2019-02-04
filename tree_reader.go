package gophy

import (
	"bufio"
	"bytes"
	"fmt"
	"os"
	"strconv"
)

//ReadTreeFromFile read a single tree from a file
func ReadTreeFromFile(tfn string) (tree *Tree) {
	f, err := os.Open(tfn)
	if err != nil {
		fmt.Println(err)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	var rt *Node
	tree = NewTree()
	for scanner.Scan() {
		ln := scanner.Text()
		if len(ln) < 2 {
			continue
		}
		rt = ReadNewickString(ln)
		tree.Instantiate(rt)
		break
	}
	return
}

//ReadTreesFromFile read multi single tree from a file
func ReadTreesFromFile(tfn string) (trees []*Tree) {
	f, err := os.Open(tfn)
	if err != nil {
		fmt.Println(err)
	}
	defer f.Close()
	trees = make([]*Tree, 0)
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		var rt *Node
		tree := NewTree()
		ln := scanner.Text()
		if len(ln) < 2 {
			continue
		}
		rt = ReadNewickString(ln)
		tree.Instantiate(rt)
		trees = append(trees, tree)
	}
	return
}

// ReadNewickString given a string it will return a pointer to the root node
func ReadNewickString(ts string) (root *Node) {
	rt := Node{nil, nil, "root", map[string]string{}, map[string]float64{},
		map[string]int{}, 0, 0., nil, false, 0., map[float64]bool{}}
	x := 0
	nc := string(ts[x : x+1])
	start := true
	cn := new(Node)
	for {
		if nc == "(" {
			if start == true {
				cn = &rt
				start = false
			} else {
				nn := Node{cn, nil, "", map[string]string{}, map[string]float64{},
					map[string]int{}, 0, 0., nil, false, 0., map[float64]bool{}}
				cn.addChild(&nn)
				cn = &nn
			}
		} else if nc == "," {
			cn = cn.Par
		} else if nc == ")" {
			cn = cn.Par
			x++
			nc = ts[x : x+1]
			if nc == "," || nc == ")" || nc == ":" || nc == "[" || nc == ";" {
				continue
			}
			var nm bytes.Buffer
			for {
				nm.WriteString(nc)
				x++
				nc = ts[x : x+1]
				if nc == "," || nc == ")" || nc == ":" || nc == "[" || nc == ";" {
					break
				}
			}
			cn.Nam = nm.String()
			x--
		} else if nc == ";" {
			break
		} else if nc == ":" {
			x++
			nc = ts[x : x+1]
			var bl bytes.Buffer
			for {
				bl.WriteString(nc)
				x++
				nc = ts[x : x+1]
				if nc == "," || nc == ")" || nc == ":" || nc == "[" || nc == ";" {
					break
				}
			}
			b, err := strconv.ParseFloat(bl.String(), 64)
			if err != nil {
				fmt.Fprintf(os.Stderr, "There is an error in branch length processing\n")
			}
			cn.Len = b
			x--
		} else {
			nn := Node{cn, nil, "", map[string]string{}, map[string]float64{},
				map[string]int{}, 0, 0., nil, false, 0., map[float64]bool{}}
			cn.addChild(&nn)
			cn = &nn
			var nm bytes.Buffer
			for {
				nm.WriteString(nc)
				x++
				nc = ts[x : x+1]
				if nc == "," || nc == ")" || nc == ":" || nc == "[" {
					break
				}
			}
			x--
			nn.Nam = nm.String()
		}
		x++
		nc = ts[x : x+1]
	}
	root = &rt
	return
}
