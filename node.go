package gophy

import (
	"bytes"
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
