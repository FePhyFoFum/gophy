package gophy

import "math"

// GetMrca get the mrca
func GetMrca(nds []*Node, root *Node) (mrca *Node) {
	mrca = nil
	traceback := make([]*Node, 0)
	first := nds[0]
	for first != root {
		first = first.Par
		traceback = append(traceback, first)
		if first.Par == nil {
			break
		}
	}
	mrca = nds[0].Par
	for _, i := range nds[1:] {
		mrca = mrcaRecurse(mrca, i, traceback)
	}
	return
}

func mrcaRecurse(nd1 *Node, nd2 *Node, path1 []*Node) (mrca *Node) {
	mrca = nil
	path := path1[NodeSlicePosition(path1, nd1):]
	parent := nd2
	for parent != nil {
		if NodeSlicePosition(path, parent) != -1 {
			mrca = parent
			return
		}
		parent = parent.Par
	}
	return
}

// SetHeights set tree height
func SetHeights(tree *Tree) {
	for _, i := range tree.Tips {
		cur := i
		h := 0.0
		going := true
		for going {
			h += cur.Len
			cur = cur.Par
			if cur == nil {
				going = false
				break
			} else {
				if h > cur.Height {
					cur.Height = h
				}
			}
		}
	}
	for _, i := range tree.Pre {
		if i != tree.Rt {
			i.Height = math.Abs(Round(i.Par.Height-i.Len, 0.5, 5)) //weird rounding thing on some machines
		}
	}
}

// NodeNamesSliceIntersects checks to see whether two node slices intersect by name
func NodeNamesSliceIntersects(a, b []*Node) (rb bool) {
	rb = false
	for _, k := range a {
		for _, l := range b {
			if k.Nam == l.Nam {
				rb = true
				return
			}
		}
	}
	return
}
