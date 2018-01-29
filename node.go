package gophy

import (
    "bytes"
    "strconv"
)

type Node struct {
	Par  *Node   //parent
	Chs []*Node  //childs
	Nam string    //name
	Len float64  //branch length
}

func (n Node) Tips() (tips []*Node) {
    x := NewNodeStack()
    x.Push(&n)
    for x.Empty() == false {
        _, c := x.Pop()
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

func (n *Node) addChild(c *Node) {
	n.Chs = append(n.Chs, c)
}

func (n Node) String() string {
    return n.Newick(false)
}
