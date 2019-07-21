package gophy

import (
	"errors"
	"strconv"
)

// Quartet are represented as map[int]bools, one for the left and one for the right
type Quartet struct {
	Lt    map[int]bool
	Rt    map[int]bool
	Lts   []map[int]bool
	Rts   []map[int]bool
	Len   float64
	Index int
}

// ConflictsWith checks whether two biparts conflict
func (b Quartet) ConflictsWith(ib Quartet) (con bool) {
	con = false
	if IntMapIntersects(ib.Rt, b.Rt) && IntMapIntersects(ib.Rt, b.Lt) {
		if IntMapIntersects(ib.Lt, b.Rt) && IntMapIntersects(ib.Lt, b.Lt) {
			con = true
			return
		}
	}
	return
}

// Match for matching quartets
func (b Quartet) Match(inq Quartet) (mat bool) {
	mat = true
	//need to match at least 2 for each
	if b.ConflictsWith(inq) {
		mat = false
		return
	}
	lmatchedl := map[int]bool{}
	lmatchedr := map[int]bool{}
	for _, i := range b.Lts {
		tl := 0
		for j := range inq.Lts {
			if IntMapIntersects(i, inq.Lts[j]) {
				tl++
				lmatchedl[j] = true
			}
		}
		for j := range inq.Rts {
			if IntMapIntersects(i, inq.Rts[j]) {
				tl++
				lmatchedr[j] = true
			}
		}
		if tl > 1 {
			mat = false
			return
		}
	}
	samed := true
	lmatched := lmatchedl
	if len(lmatchedl) > 0 && len(lmatchedr) > 0 {
		mat = false
		return
	}
	if len(lmatchedr) > 0 {
		samed = false
		lmatched = lmatchedr
	}
	if len(lmatched) < 2 {
		mat = false
		return
	}
	rmatchedl := map[int]bool{}
	rmatchedr := map[int]bool{}
	for _, i := range b.Rts {
		tl := 0
		for j := range inq.Lts {
			if IntMapIntersects(i, inq.Lts[j]) {
				tl++
				rmatchedl[j] = true
			}
		}
		for j := range inq.Rts {
			if IntMapIntersects(i, inq.Rts[j]) {
				tl++
				rmatchedr[j] = true
			}
		}
		if tl > 1 {
			mat = false
			return
		}
	}
	if len(rmatchedl) > 0 && len(rmatchedr) > 0 {
		mat = false
		return
	}
	rmatched := rmatchedr
	if samed && len(rmatchedl) > 0 {
		mat = false
		return
	}
	if samed == false && len(rmatchedr) > 0 {
		mat = false
		return
	}
	if samed == false {
		rmatched = rmatchedl
	}
	if len(rmatched) < 2 {
		mat = false
		return
	}
	mat = true
	return
}

// StringWithNames converts the ints to the the strings from nmmap
func (b Quartet) StringWithNames(nmmap map[int]string) (ret string) {
	for cb, i := range b.Lts {
		c := 0
		for j := range i {
			ret += nmmap[j]
			if c < len(i)-1 {
				ret += ","
			}
			c++
		}
		if cb < len(b.Lts)-1 {
			ret += " | "
		}
	}
	ret += " === "
	for cb, i := range b.Rts {
		c := 0
		for j := range i {
			ret += nmmap[j]
			if c < len(i)-1 {
				ret += ","
			}
			c++
		}
		if cb < len(b.Lts)-1 {
			ret += " | "
		}
	}
	s := strconv.FormatFloat(b.Len, 'f', 1, 64)
	ret += " (" + s + ")"
	return
}

//GetQuartet for node
func GetQuartet(nd *Node, tree Tree, maptips map[string]int) (Quartet, error) {
	if len(nd.Chs) == 0 || nd == tree.Rt {
		bp := Quartet{}
		return bp, errors.New("root or child")
	}
	rights := []map[int]bool{}
	lefts := []map[int]bool{}
	right := map[int]bool{}
	for _, i := range nd.Chs {
		r := map[int]bool{}
		for _, j := range i.GetTipNames() {
			r[maptips[j]] = true
			right[maptips[j]] = true
		}
		rights = append(rights, r)
	}
	p := nd.Par
	for _, i := range p.Chs {
		ilvs := map[int]bool{}
		for _, j := range i.GetTipNames() {
			ilvs[maptips[j]] = true
		}
		if IntMapIntersects(ilvs, right) {
			continue
		} else {
			lefts = append(lefts, ilvs)
		}
	}
	if p != tree.Rt {
		out := map[int]bool{}
		plvs := map[int]bool{}
		for _, j := range tree.Rt.GetTipNames() {
			out[maptips[j]] = true
		}
		for _, j := range p.GetTipNames() {
			plvs[maptips[j]] = true
		}
		inter := IntMapDifferenceRet(out, plvs)
		interm := map[int]bool{}
		for _, j := range inter {
			interm[j] = true
		}
		lefts = append(lefts, interm)
	}
	bp := Quartet{Lts: lefts, Rts: rights, Len: nd.Len}
	return bp, nil
}

//GetQuartets return slice of quartets
func GetQuartets(tree Tree, maptips map[string]int) (bps []Quartet) {
	for _, i := range tree.Pre {
		q, err := GetQuartet(i, tree, maptips)
		if err == nil {
			bps = append(bps, q)
		}
	}
	return bps
}
