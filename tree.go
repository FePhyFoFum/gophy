package gophy

// Tree minimal phylogenetic tree struct
// don't need this, can use root but here if you want these convienence functions below
type Tree struct {
	Rt   *Node
	Post []*Node
	Pre  []*Node
	Tips []*Node
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
	t.populatePreorder(t.Rt)
	t.populatePostorder(t.Rt)
}
