package gophy

import (
	"errors"
	"sync"
)

// NodeStack makes a node stack for pushing and pulling.
type NodeStack struct {
	lock sync.Mutex // you don't have to do this if you don't want thread safety
	s    []*Node
}

// NewNodeStack returns a pointer to a node stack
func NewNodeStack() *NodeStack {
	return &NodeStack{sync.Mutex{}, make([]*Node, 0)}
}

// Push push the node into the stack
func (s *NodeStack) Push(v *Node) {
	s.lock.Lock()
	defer s.lock.Unlock()

	s.s = append(s.s, v)
}

// Empty is the node stack empty?
func (s *NodeStack) Empty() (ret bool) {
	ret = true
	if len(s.s) > 0 {
		ret = false
	}
	return
}

// Pop a node off the stack
func (s *NodeStack) Pop() (error, *Node) {
	s.lock.Lock()
	defer s.lock.Unlock()

	l := len(s.s)
	if l == 0 {
		return errors.New("Empty NodeStack"), nil
	}

	res := s.s[l-1]
	s.s = s.s[:l-1]
	return nil, res
}
