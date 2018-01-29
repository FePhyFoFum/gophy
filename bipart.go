package gophy

/*
Biparts are represented as map[int]bools, one for the left and one for the right
*/
type Bipart struct {
	Lt map[int]bool
	Rt map[int]bool
}

func (b Bipart) StringWithNames(nmmap map[int]string) (ret string) {
	for n, _ := range b.Lt {
		ret += nmmap[n] + " "
	}
	ret += "|"
	for n, _ := range b.Rt {
		ret += " " + nmmap[n]
	}
	return
}

func (b Bipart) Equals(ib Bipart) (eq bool) {
	eq = false
	if len(b.Rt) == len(ib.Rt) && len(b.Lt) == len(ib.Lt) {
		for m, _ := range b.Rt {
			if _, ok := ib.Rt[m]; !ok {
				return
			}
		}
		for m, _ := range b.Lt {
			if _, ok := ib.Lt[m]; !ok {
				return
			}
		}
		eq = true
		return
	} else if len(b.Rt) == len(ib.Lt) && len(b.Lt) == len(ib.Rt) {
		for m, _ := range b.Rt {
			if _, ok := ib.Lt[m]; !ok {
				return
			}
		}
		for m, _ := range b.Lt {
			if _, ok := ib.Rt[m]; !ok {
				return
			}
		}
		eq = true
		return
	}
	return
}

func (b Bipart) Conflicts(ib Bipart) (con bool) {
	con = false
	if IntMapIntersects(ib.Rt, b.Rt) && IntMapIntersects(ib.Rt, b.Lt) {
		if IntMapIntersects(ib.Lt, b.Rt) && IntMapIntersects(ib.Lt, b.Lt) {
			con = true
			return
		}
	}
	return
}

func (b Bipart) Concord(ib Bipart) (con bool) {
	return
}

func IntMapIntersects(s1 map[int]bool, s2 map[int]bool) (in bool) {
	in = false
	for k, _ := range s1 {
		if s2[k] {
			in = true
			return
		}
	}
	return
}

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

func IntSliceContains(is []int, s int) (rb bool) {
	for _, a := range is {
		if a == s {
			return true
		}
	}
	return false
}

func PConflicts(bps []Bipart, jobs <-chan []int, results chan<- []int) {
	for j := range jobs {
		in1, in2 := j[0], j[1]
		b := 0
		if bps[in1].Conflicts(bps[in2]) {
			b = 1
		}
		results <- []int{in1, in2, b}
	}
}


