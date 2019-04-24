package gophy

import (
	"fmt"
	"math"
)

// GetMrca get the mrca
func GetMrca(nds []*Node, root *Node) (mrca *Node) {
	if len(nds) == 1 {
		return nds[0]
	}
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

// StochasticNNI
func StochasticNNI() {

}

// to get the branches that would be NNIs, use the function and
// then send the relevant nodes here
// SwapBranch
func SwapBranch(nd1 *Node, nd2 *Node) bool {
	if nd1.Par == nil || nd2.Par == nil {
		return false
	}
	par1 := nd1.Par
	par2 := nd2.Par
	nd1.Par = par2
	nd2.Par = par1
	par1.removeChild(nd1)
	par2.removeChild(nd2)
	par1.addChild(nd2)
	par2.addChild(nd1)
	return true
}

func TritomyRoot(tr *Tree) {
	curroot := tr.Rt
	if len(curroot.Chs) > 2 {
		return
	}
	if len(curroot.Chs[0].Chs) > 0 { //internal
		currootCH := curroot.Chs[0]
		nbl := currootCH.Len
		curroot.Chs[1].Len = curroot.Chs[1].Len + nbl
		curroot.removeChild(currootCH)
		for _, i := range currootCH.Chs {
			curroot.addChild(i)
			i.Par = curroot
		}
	} else {
		currootCH := curroot.Chs[1]
		nbl := currootCH.Len
		curroot.Chs[0].Len = curroot.Chs[0].Len + nbl
		curroot.removeChild(currootCH)
		for _, i := range currootCH.Chs {
			curroot.addChild(i)
			i.Par = curroot
		}
	}
}

func Reroot(inroot *Node, tr *Tree) {
	tempParent := inroot.Par
	newRoot := new(Node)
	newRoot.addChild(inroot)
	inroot.Par = newRoot
	tempParent.removeChild(inroot)
	tempParent.addChild(newRoot)
	newRoot.Par = tempParent
	newRoot.Len = inroot.Len / 2.
	inroot.Len = inroot.Len / 2.
	processReRoot(newRoot)
	tr.Rt = newRoot
}

func processReRoot(node *Node) {
	if node.Par == nil || len(node.Chs) == 0 {
		return
	}
	if node.Par != nil {
		processReRoot(node.Par)
	}
	// Exchange branch label, length et cetera
	exchangeInfo(node.Par, node)
	// Rearrange topology
	parent := node.Par
	node.addChild(parent)
	parent.removeChild(node)
	parent.Par = node
}

func exchangeInfo(node1 *Node, node2 *Node) {
	swaps := ""
	swapd := 0.0
	swaps = node1.Nam
	node1.Nam = node2.Nam
	node2.Nam = swaps
	swapd = node1.Len
	node1.Len = node2.Len
	node2.Len = swapd
}

//NNIMoves looks at the root and returns the NNIs
func NNIMoves(tr *Tree) [][]*Node {
	if len(tr.Rt.Chs) > 3 {
		fmt.Errorf("needs to be a tritomy root.")
		return nil
	}
	x0 := len(tr.Rt.Chs[0].Chs)
	x1 := len(tr.Rt.Chs[1].Chs)
	x2 := len(tr.Rt.Chs[2].Chs)
	//  1  2 34
	//   \ | x
	//    \|/
	//need to label these
	//x2 := len(tr.Rt.Chs[2].Chs)
	var nd1 *Node
	var nd2 *Node
	var nd3 *Node
	var nd4 *Node
	moves := [][]*Node{}
	if x2 > 1 {
		nd1 = tr.Rt.Chs[0]
		nd2 = tr.Rt.Chs[1]
		nd3 = tr.Rt.Chs[2].Chs[0]
		nd4 = tr.Rt.Chs[2].Chs[1]
		tm := []*Node{nd1, nd3}
		moves = append(moves, tm)
		tm2 := []*Node{nd2, nd4}
		moves = append(moves, tm2)
	}
	if x1 > 1 {
		nd1 = tr.Rt.Chs[0]
		nd2 = tr.Rt.Chs[2]
		nd3 = tr.Rt.Chs[1].Chs[0]
		nd4 = tr.Rt.Chs[1].Chs[1]
		tm := []*Node{nd1, nd3}
		moves = append(moves, tm)
		tm2 := []*Node{nd2, nd4}
		moves = append(moves, tm2)
	}
	if x0 > 1 {
		nd1 = tr.Rt.Chs[1]
		nd2 = tr.Rt.Chs[2]
		nd3 = tr.Rt.Chs[0].Chs[0]
		nd4 = tr.Rt.Chs[0].Chs[1]
		tm := []*Node{nd1, nd3}
		moves = append(moves, tm)
		tm2 := []*Node{nd2, nd4}
		moves = append(moves, tm2)
	}
	fmt.Println(moves)
	return moves
}
