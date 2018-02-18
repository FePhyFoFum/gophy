package gophy

// Tree minimal phylogenetic tree struct
// don't need this, can use root but here if you want these convienence functions below
type Tree struct {
	Rt    *Node
	Post  []*Node
	Pre   []*Node
	Tips  []*Node
	Index int // if there is an identifying index
}

func (t *Tree) populatePrePostIt(rt *Node) {
	stk := NewNodeStack()
	stk.Push(rt)

	for stk.Empty() == false {
		cur, _ := stk.Pop()
		t.Pre = append(t.Pre, cur)
		t.Post = append(t.Post, cur) // will be reversed
		for _, n := range cur.Chs {
			stk.Push(n)
		}
		if len(cur.Chs) == 0 {
			t.Tips = append(t.Tips, cur)
		}
	}
	for i, j := 0, len(t.Post)-1; i < j; i, j = i+1, j-1 {
		t.Post[i], t.Post[j] = t.Post[j], t.Post[i]
	}
}

func (t *Tree) populatePostorder(rt *Node) {
	for _, n := range rt.Chs {
		t.populatePostorder(n)
	}
	t.Post = append(t.Post, rt)
	if len(rt.Chs) == 0 {
		t.Tips = append(t.Tips, rt)
	}
}

func (t *Tree) populatePreorder(rt *Node) {
	t.Pre = append(t.Pre, rt)
	for _, n := range rt.Chs {
		t.populatePreorder(n)
	}
}

// Instantiate will preorder and postorder
func (t *Tree) Instantiate(rt *Node) {
	t.Rt = rt
	t.populatePrePostIt(t.Rt)
	//t.populatePreorder(t.Rt)
	//t.populatePostorder(t.Rt)
}
