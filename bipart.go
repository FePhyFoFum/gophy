package gophy

// Bipart are represented as map[int]bools, one for the left and one for the right
type Bipart struct {
	Lt map[int]bool
	Rt map[int]bool
}

// StringWithNames converts the ints to the the strings from nmmap
func (b Bipart) StringWithNames(nmmap map[int]string) (ret string) {
	for n := range b.Lt {
		ret += nmmap[n] + " "
	}
	ret += "|"
	for n := range b.Rt {
		ret += " " + nmmap[n]
	}
	return
}

// NewickWithNames does similar things to StringWithNames but sends a newick back
func (b Bipart) NewickWithNames(nmmap map[int]string) (ret string) {
	ret += "(("
	count := 0
	for n := range b.Lt {
		ret += nmmap[n]
		if count < len(b.Lt)-1 {
			ret += ","
		}
		count++
	}
	ret += ")"
	for n := range b.Rt {
		ret += "," + nmmap[n]
	}
	ret += ");"
	return
}

// Equals checks whether two biparts are the same
func (b Bipart) Equals(ib Bipart) (eq bool) {
	eq = false
	if len(b.Rt) == len(ib.Rt) && len(b.Lt) == len(ib.Lt) {
		for m := range b.Rt {
			if _, ok := ib.Rt[m]; !ok {
				return
			}
		}
		for m := range b.Lt {
			if _, ok := ib.Lt[m]; !ok {
				return
			}
		}
		eq = true
		return
	} else if len(b.Rt) == len(ib.Lt) && len(b.Lt) == len(ib.Rt) {
		for m := range b.Rt {
			if _, ok := ib.Lt[m]; !ok {
				return
			}
		}
		for m := range b.Lt {
			if _, ok := ib.Rt[m]; !ok {
				return
			}
		}
		eq = true
		return
	}
	return
}

// ConflictsWith checks whether two biparts conflict
func (b Bipart) ConflictsWith(ib Bipart) (con bool) {
	con = false
	if IntMapIntersects(ib.Rt, b.Rt) && IntMapIntersects(ib.Rt, b.Lt) {
		if IntMapIntersects(ib.Lt, b.Rt) && IntMapIntersects(ib.Lt, b.Lt) {
			con = true
			return
		}
	}
	return
}

// ConcordantWith tests whether something is concordant (not conflicting or nested, etc)
func (b Bipart) ConcordantWith(ib Bipart) (con bool) {
	con = false
	if IntMapIntersects2(ib.Rt, b.Rt) && IntMapIntersects2(ib.Lt, b.Lt) {
		if IntMapIntersects(ib.Rt, b.Lt) == false {
			if IntMapIntersects(ib.Lt, b.Rt) == false {
				con = true
				return
			}
		} else {
			return
		}
	}
	if IntMapIntersects2(ib.Lt, b.Rt) && IntMapIntersects2(ib.Rt, b.Lt) {
		if IntMapIntersects(ib.Rt, b.Lt) == false {
			if IntMapIntersects(ib.Lt, b.Rt) == false {
				con = true
				return
			}
		} else {
			return
		}
	}
	return
}

// CompatibleWith checks that it isn't conflicting but can be nested
func (b Bipart) CompatibleWith(ib Bipart) (con bool) {
	con = true
	return
}

// IntMapIntersects checks to see if the two map[int]bool intersect (in the set sense)
func IntMapIntersects(s1 map[int]bool, s2 map[int]bool) (in bool) {
	in = false
	for k := range s1 {
		if s2[k] {
			in = true
			return
		}
	}
	return
}

// IntMapIntersects2 checks to see if the two map[int]bool intersect (in the set sense)
// with at least 2 matches
func IntMapIntersects2(s1 map[int]bool, s2 map[int]bool) (in bool) {
	in = false
	count := 0
	for k := range s1 {
		if s2[k] {
			count++
			if count >= 2 {
				in = true
				return
			}
		}
	}
	return
}

// BipartSliceContains checks to see if the bipart slice contains the bipart and returns the index
func BipartSliceContains(bps []Bipart, bp Bipart) (ind int) {
	ind = -1
	for i, value := range bps {
		if value.Equals(bp) {
			ind = i
			return
		}
	}
	return
}

// IntSliceContains checks to see if the int slice contains an int and returns the bool
func IntSliceContains(is []int, s int) (rb bool) {
	for _, a := range is {
		if a == s {
			return true
		}
	}
	return false
}

// PConflicts is a parallel conflict check. The slice is sent. The jobs are the two indices to check.
// The results are the two indicies and an int 1 for conflict 0 for no conflict
func PConflicts(bps []Bipart, jobs <-chan []int, results chan<- []int) {
	for j := range jobs {
		in1, in2 := j[0], j[1]
		b := 0
		if bps[in1].ConflictsWith(bps[in2]) {
			b = 1
		}
		results <- []int{in1, in2, b}
	}
}

// PConflictsCompTree is similar to PConflict but the first index refers to the first Bipart slice and the second referts
// to the second Bipart slice
func PConflictsCompTree(bps []Bipart, comptreebps []Bipart, jobs <-chan []int, results chan<- []int) {
	for j := range jobs {
		in1, in2 := j[0], j[1]
		b := 0
		if bps[in1].ConflictsWith(comptreebps[in2]) {
			b = 1
		}
		results <- []int{in1, in2, b}
	}
}

// PConcordance is similar to the other comparison code but for concordance. The input jobs are the i, j for the bipart
// comparisons. The results are the i, j, and 0 for not concordant and 1 for concordant
func PConcordance(bps []Bipart, jobs <-chan []int, results chan<- []int) {
	for j := range jobs {
		in1, in2 := j[0], j[1]
		b := 0
		if bps[in1].ConcordantWith(bps[in2]) {
			b = 1
		}
		results <- []int{in1, in2, b}
	}
}
