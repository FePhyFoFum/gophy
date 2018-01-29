package gophy

import (
    "sync"
    "errors"
)

type NodeStack struct {
     lock sync.Mutex // you don't have to do this if you don't want thread safety
     s []*Node
}

func NewNodeStack() *NodeStack {
    return &NodeStack {sync.Mutex{}, make([]*Node,0), }
}

func (s *NodeStack) Push(v * Node) {
    s.lock.Lock()
    defer s.lock.Unlock()

    s.s = append(s.s, v)
}

func (s *NodeStack) Empty() (ret bool) {
    ret = true
    if len(s.s) > 0 {
        ret = false
    }
    return
}

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


