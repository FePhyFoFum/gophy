package gophy_test

import (
	"fmt"
	"testing"

	"github.com/FePhyFoFum/gophy"
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
	glt[3] = true
	glt[4] = true
	glt[5] = true
	grt := make(map[int]bool)
	grt[0] = true
	grt[1] = true
	grt[2] = true
	b := gophy.Bipart{Lt: blt, Rt: brt}
	g := gophy.Bipart{Lt: glt, Rt: grt}
	fmt.Println(b)
	fmt.Println(g)
	fmt.Println(b.ConcordantWith(g))
	if b.ConcordantWith(g) == false {
		t.Fail()
	}
	if g.ConcordantWith(b) == false {
		t.Fail()
	}
}

func TestConflictsWith(t *testing.T) {
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
	glt[3] = true
	glt[5] = true
	grt := make(map[int]bool)
	grt[4] = true
	grt[1] = true
	grt[2] = true
	b := gophy.Bipart{Lt: blt, Rt: brt}
	g := gophy.Bipart{Lt: glt, Rt: grt}
	if b.ConflictsWith(g) == false {
		t.Fail()
	}
	if g.ConflictsWith(b) == false {
		t.Fail()
	}
}

func TestEquals(t *testing.T) {
	blt := make(map[int]bool)
	blt[0] = true
	blt[1] = true
	blt[2] = true
	brt := make(map[int]bool)
	brt[3] = true
	brt[4] = true
	brt[5] = true
	glt := make(map[int]bool)
	glt[3] = true
	glt[4] = true
	glt[5] = true
	grt := make(map[int]bool)
	grt[0] = true
	grt[1] = true
	grt[2] = true
	b := gophy.Bipart{Lt: blt, Rt: brt}
	g := gophy.Bipart{Lt: grt, Rt: glt}
	if b.Equals(g) == false {
		t.Fail()
	}
	if g.Equals(b) == false {
		t.Fail()
	}
}
