package gophy

import (
	"bytes"
	"fmt"
	"os"
	"strconv"
)

// ReadNewickString given a string it will return a pointer to the root node
func ReadNewickString(ts string) (root *Node) {
	rt := Node{nil, nil, "root", map[string]string{}, 0., nil, false, 0., map[float64]bool{}}
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
				nn := Node{cn, nil, "", map[string]string{}, 0., nil, false, 0., map[float64]bool{}}
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
			nn := Node{cn, nil, "", map[string]string{}, 0., nil, false, 0., map[float64]bool{}}
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
