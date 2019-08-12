package gophy

import (
	"bytes"
	"fmt"
	"os"
	"strconv"
)

// Node minimal node struct
type Node struct {
	Par       *Node   //parent
	Chs       []*Node //childs
	Nam       string  //name
	SData     map[string]string
	FData     map[string]float64
	IData     map[string]int
	Num       int
	Len       float64     //branch length
	Data      [][]float64 // [site][states]
	Marked    bool        //just for like calculations
	Height    float64
	MarkedMap map[float64]bool //the float is for an id for the query
	//anc+bl
	// X------>--X
	// TP        RT
	//
	// X--<------X
	// RvTp      RV
	TpConds   [][]float64 //[site][states]
	RtConds   [][]float64
	RvConds   [][]float64
	RvTpConds [][]float64
	//TODO: need comments for these below
	FAD         float64
	LAD         float64
	FINDS       float64
	TimeLen     float64
	ContData    []float64
	Mis         []bool
	PruneLen    float64
	ConPruneLen []float64 // prevent race condition when calculating BM likelihood
	//BMPruneLen  []float64
	BMLen    float64
	LL       []float64
	Anc      bool
	ClustLen map[int]float64
	//ClustPruneLen map[int]float64 // [][]float64
	//
}

// Walk just a simple walker with chans
// can use this with this
/*
ch := make(chan *Node)
go func() {
	rt.Walk(ch)
	close(ch)
}()
for n := range ch {
	t.Post = append(t.Post, n)
	if len(n.Chs) == 0 {
		t.Tips = append(t.Tips, n)
	}
}*/
func (n *Node) Walk(ch chan *Node) {
	if n == nil {
		return
	}
	for i := range n.Chs {
		n.Chs[i].Walk(ch)
	}
	ch <- n
}

// GetTips returns a slice with node pointers
func (n Node) GetTips() (tips []*Node) {
	x := NewNodeStack()
	x.Push(&n)
	for x.Empty() == false {
		c, _ := x.Pop()
		if len(c.Chs) == 0 {
			tips = append(tips, c)
		} else {
			for _, h := range c.Chs {
				x.Push(h)
			}
		}
	}
	return
}

// GetTipNames returns a slice with node pointers
func (n Node) GetTipNames() (tips []string) {
	x := NewNodeStack()
	x.Push(&n)
	for x.Empty() == false {
		c, _ := x.Pop()
		if len(c.Chs) == 0 {
			tips = append(tips, c.Nam)
		} else {
			for _, h := range c.Chs {
				x.Push(h)
			}
		}
	}
	return
}

// NewickChronogram returns a string newick with the branch lengths scaled to time
func (n Node) NewickChronogram() (ret string) {
	bl := true
	var buffer bytes.Buffer
	for in, cn := range n.Chs {
		if in == 0 {
			buffer.WriteString("(")
		}
		buffer.WriteString(cn.NewickChronogram())
		if bl == true {
			s := strconv.FormatFloat(cn.TimeLen, 'f', -1, 64)
			buffer.WriteString(":")
			buffer.WriteString(s)
		}
		if in == len(n.Chs)-1 {
			buffer.WriteString(")")
		} else {
			buffer.WriteString(",")
		}
	}
	buffer.WriteString(n.Nam)
	ret = buffer.String()
	return
}

// Newick returns a string newick
func (n Node) Newick(bl bool) (ret string) {
	var buffer bytes.Buffer
	for in, cn := range n.Chs {
		if in == 0 {
			buffer.WriteString("(")
		}
		buffer.WriteString(cn.Newick(bl))
		if bl == true {
			s := strconv.FormatFloat(cn.Len, 'f', -1, 64)
			buffer.WriteString(":")
			buffer.WriteString(s)
		}
		if in == len(n.Chs)-1 {
			buffer.WriteString(")")
		} else {
			buffer.WriteString(",")
		}
	}
	buffer.WriteString(n.Nam)
	ret = buffer.String()
	return
}

// BMPhylogram returns a string newick with brownian motion branch lengths
//TODO can this and the TimeLen merge with just a boolean?
func (n Node) BMPhylogram() (ret string) {
	bl := true
	var buffer bytes.Buffer
	for in, cn := range n.Chs {
		if in == 0 {
			buffer.WriteString("(")
		}
		buffer.WriteString(cn.BMPhylogram())
		if bl == true {
			s := strconv.FormatFloat(cn.BMLen, 'f', -1, 64)
			buffer.WriteString(":")
			buffer.WriteString(s)
		}
		if in == len(n.Chs)-1 {
			buffer.WriteString(")")
		} else {
			buffer.WriteString(",")
		}
	}
	buffer.WriteString(n.Nam)
	ret = buffer.String()
	return
}

// NewickFloatBL returns a string newick with branch lengths of the data in FData[fl]
func (n Node) NewickFloatBL(fl string) (ret string) {
	var buffer bytes.Buffer
	for in, cn := range n.Chs {
		if in == 0 {
			buffer.WriteString("(")
		}
		buffer.WriteString(cn.NewickFloatBL(fl))
		if _, ok := cn.FData[fl]; ok {
			s := strconv.FormatFloat(cn.FData[fl], 'f', -1, 64)
			buffer.WriteString(":")
			buffer.WriteString(s)
		}
		if in == len(n.Chs)-1 {
			buffer.WriteString(")")
		} else {
			buffer.WriteString(",")
		}
	}
	buffer.WriteString(n.Nam)
	ret = buffer.String()
	return
}

// NewickPaint returns a string newick
func (n Node) NewickPaint(bl bool, rid float64) (ret string) {
	var buffer bytes.Buffer
	painted := make([]*Node, 0)
	for _, cn := range n.Chs {
		if _, ok := cn.MarkedMap[rid]; ok {
			painted = append(painted, cn)
		}
	}
	for in, cn := range painted {
		if in == 0 {
			buffer.WriteString("(")
		}
		buffer.WriteString(cn.NewickPaint(bl, rid))
		if bl == true {
			s := strconv.FormatFloat(cn.Len, 'f', -1, 64)
			buffer.WriteString(":")
			buffer.WriteString(s)
		}
		if in == len(painted)-1 {
			buffer.WriteString(")")
		} else {
			buffer.WriteString(",")
		}
	}
	if _, ok := n.MarkedMap[rid]; ok {
		buffer.WriteString(n.Nam)
	}
	ret = buffer.String()
	return
}

func (n *Node) addChild(c *Node) {
	n.Chs = append(n.Chs, c)
}

func (n *Node) removeChild(c *Node) {
	s := -1
	for i, j := range n.Chs {
		if j == c {
			s = i
		}
	}
	if s == -1 {
		return
	}
	n.Chs = append(n.Chs[:s], n.Chs[s+1:]...)
}

func (n Node) String() string {
	return n.Newick(false)
}

//RerootLS reroots all the nodes represented in a graph on n
func (n *Node) RerootBM(oldroot *Node) *Node {
	if n == oldroot {
		fmt.Println("you are trying to reroot on the current root!")
	}
	nnodes := 0
	oldroot.NNodes(&nnodes)
	var pathnodes = make([]*Node, nnodes)
	//var pathnodes []*Node
	curnode := n
	pathlen := 0 //this will count the number of nodes between the newroot and the oldroot
	for ind := range pathnodes {
		pathnodes[ind] = curnode
		//pathnodes = append(pathnodes,curnode)
		if curnode == oldroot {
			break
		}
		pathlen++
		curnode = curnode.Par
	}
	var newpar *Node
	for i := pathlen; i >= 1; i-- {
		newpar = pathnodes[i-1]
		curnode = pathnodes[i]
		curnode.removeChild(newpar)
		newpar.addChild(curnode)
		curnode.BMLen = newpar.BMLen
	}
	//curnode = nil
	//newpar = nil
	n.BMLen = 0.0
	return n
}

//Reroot reroots all the nodes represented in a graph on n
func (n *Node) Reroot(oldroot *Node) *Node {
	if n == oldroot {
		fmt.Println("you are trying to reroot on the current root!")
	}
	nnodes := 0
	oldroot.NNodes(&nnodes)
	var pathnodes = make([]*Node, nnodes)
	//var pathnodes []*Node
	curnode := n
	pathlen := 0 //this will count the number of nodes between the newroot and the oldroot
	for ind := range pathnodes {
		pathnodes[ind] = curnode
		//pathnodes = append(pathnodes,curnode)
		if curnode == oldroot {
			break
		}
		pathlen++
		curnode = curnode.Par
	}
	var newpar *Node
	for i := pathlen; i >= 1; i-- {
		newpar = pathnodes[i-1]
		curnode = pathnodes[i]
		curnode.removeChild(newpar)
		newpar.addChild(curnode)
		curnode.Len = newpar.Len
	}
	//curnode = nil
	//newpar = nil
	n.Len = 0.0
	return n
}

//NNodes is a helper method that will return the number of internal nodes descending from n (including n)
func (n *Node) NNodes(count *int) {
	*count++
	for _, ch := range n.Chs {
		ch.NNodes(count)
	}
}

//PreorderArray will return a preordered array of all the nodes in a tree
//TODO: Can remove after some refactoring but here while getting off the ground
func (n *Node) PreorderArray() (ret []*Node) {
	var buffer []*Node
	buffer = append(buffer, n)
	for _, cn := range n.Chs {
		for _, cret := range cn.PreorderArray() {
			buffer = append(buffer, cret)
		}
	}
	ret = buffer
	return
}

func (n *Node) GetBackbone(higherNode *Node) (backbone []*Node) {
	cur := n
	for {
		if cur.Par == nil && cur != higherNode {
			fmt.Println("failed at getting backbone. higher node is probably not actually above the lower node")
			break
		}
		backbone = append(backbone, cur)
		if cur.Par == higherNode {
			break
		}
	}
	return
}

//GetSib returns the sibling of a node
func (n *Node) GetSib() *Node {
	if n.Nam == "root" {
		fmt.Println("Root node has no sibling")
		os.Exit(0)
	}
	par := n.Par
	var sib *Node
	if len(par.Chs) != 2 {
		if len(par.Chs) == 1 {
			fmt.Println("Singleton encountered in tree")
			os.Exit(0)
		} else {
			fmt.Println("Multifurcation found in tree")
			os.Exit(0)
		}
	}
	for _, c := range par.Chs {
		if c != n {
			sib = c
		}
	}
	if sib == nil {
		fmt.Println("something is messed up with the tree. can't find a sister node for node", n)
		os.Exit(0)
	}
	return sib
}
