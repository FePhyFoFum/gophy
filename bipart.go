package gophy

import (
	"fmt"
	"sort"
)

// Bipart are represented as map[int]bools, one for the left and one for the right
type Bipart struct {
	Lt          map[int]bool
	Rt          map[int]bool
	Ct          int   // counts
	TreeIndices []int //index of which trees this is in
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

// Equals trying to be faster
func (b Bipart) Equals(ib Bipart) (eq bool) {
	eq = false
	if len(b.Rt) == len(ib.Rt) && len(b.Lt) == len(ib.Lt) {
		// if the lengths of left and right are the same, we have to check a special case
		// where they could be reversed
		if len(b.Rt) == len(b.Lt) {
			reverse := false
			for m := range b.Rt {
				if _, ok := ib.Rt[m]; !ok {
					reverse = true
				}
				break
			}
			if reverse == false {
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
			}
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
		if IntMapIntersects(ib.Rt, b.Rt) == false {
			if IntMapIntersects(ib.Lt, b.Lt) == false {
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

// OutputEdges just print the edges
// bpbpts bipart int to tree index int list
// mapints are int to string names for the taxa
// bps list of biparts
// ntrees number of trees
func OutputEdges(bpbpts map[int][]int, mapints map[int]string, bps []Bipart, ntrees int) {
	for i, b := range bps {
		fmt.Println(b.StringWithNames(mapints), len(bpbpts[i]), float64(len(bpbpts[i]))/float64(ntrees))
	}
}

// CompareTreeToBiparts take biparts from a set , comparetreebps, and compre them to another set bps
func CompareTreeToBiparts(bps []Bipart, comptreebps []Bipart, workers int, mapints map[int]string) {
	jobs := make(chan []int, len(bps)*len(comptreebps))
	results := make(chan []int, len(bps)*len(comptreebps))
	for w := 1; w <= workers; w++ {
		go PConflictsCompTree(bps, comptreebps, jobs, results)
	}
	for j := range comptreebps {
		for i := range bps {
			jobs <- []int{i, j}
		}
	}
	close(jobs)
	confs := make(map[int][]int)   // key is bipart and value are the conflicts
	allconfs := make(map[int]bool) // list of all the conflicting biparts
	for range comptreebps {
		for range bps {
			x := <-results
			if x[2] == 1 {
				confs[x[1]] = append(confs[x[1]], x[0])
				allconfs[x[0]] = true
			}
		}
	}
	/*
	 going to make a set of concordance biparts of the set of conflicting biparts
	*/
	jobs = make(chan []int, len(allconfs)*len(allconfs))
	results = make(chan []int, len(allconfs)*len(allconfs))

	for w := 1; w <= workers; w++ {
		go PConcordance(bps, jobs, results)
	}

	for i := range allconfs {
		for j := range allconfs {
			jobs <- []int{i, j}
		}
	}
	close(jobs)
	bpsConcCounts := make(map[int]int) // key bipart index, value number of concordant bps
	for i := range allconfs {
		for j := range allconfs {
			if i < j {
				x := <-results
				if x[2] == 1 {
					if _, ok := bpsConcCounts[x[0]]; ok {
						bpsConcCounts[x[0]]++
					} else {
						bpsConcCounts[x[0]] = 0
					}
					if _, ok := bpsConcCounts[x[1]]; ok {
						bpsConcCounts[x[1]]++
					} else {
						bpsConcCounts[x[1]] = 0
					}
				}
			}
		}
	}
	/*
	   sorting the results so that the larger bps are listed first. stop printing after a few.
	   add a sys command for listing all the results
	*/
	for x, y := range confs {
		fmt.Print(comptreebps[x].NewickWithNames(mapints) + "\n")
		n := map[int][]int{}
		var a []int
		for _, v := range y {
			//n[bpsCounts[v]] = append(n[bpsCounts[v]], v)
			n[bpsConcCounts[v]] = append(n[bpsConcCounts[v]], v)
		}
		for k := range n {
			a = append(a, k)
		}
		sort.Sort(sort.Reverse(sort.IntSlice(a)))
		count := 0
		for _, k := range a {
			for _, s := range n[k] {
				//s is the bps index, k is the count
				//fmt.Print(" ", bpsCounts[s], " "+bps[s].NewickWithNames(mapints)+"\n")
				fmt.Print("  ", bps[s].Ct, " ", bpsConcCounts[s], " "+bps[s].NewickWithNames(mapints)+"\n")
				// fmt.Println(s, k)
				if count >= 10 {
					break
				}
				count++
			}
			if count >= 5 {
				break
			}
		}
	}
}
