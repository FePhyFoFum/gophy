package gophy_test

import (
	"gophy"
	"testing"
)

func TestConcordantWith(t *testing.T) {
	blt := make(map[int]bool)
	blt[0] = true
	blt[1] = true
	blt[2] = true
	brt := make(map[int]bool)
	brt[3] = true
	brt[4] = true
	brt[5] = true
	glt := make(map[int]bool)
	glt[0] = true
	glt[1] = true
	grt := make(map[int]bool)
	grt[3] = true
	grt[4] = true
	grt[5] = true
	b := gophy.Bipart{blt, brt}
	g := gophy.Bipart{glt, grt}
	if b.ConcordantWith(g) == false {
		t.Fail()
	}
	if g.ConcordantWith(b) == false {
		t.Fail()
	}
}
