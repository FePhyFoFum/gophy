package gophy

import (
	"fmt"
	"sort"
	"strconv"
)

// Bipart are represented as map[int]bools, one for the left and one for the right
type Bipart struct {
	Lt          map[int]bool
	Rt          map[int]bool
	Ct          int     // counts
	TreeIndices []int   // index of which trees this is in
	Nds         []*Node // nodes associated with the bipart
	Index       int     // just a unique id
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
		// there must be at least one difference in the Trees so it isn't just the same tree
		if CalcSliceIntDifferenceInt(bps[in1].TreeIndices, bps[in2].TreeIndices) > 0 {
			if bps[in1].ConcordantWith(bps[in2]) {
				b = 1
			}
		}
		results <- []int{in1, in2, b}
	}
}

// PConcordanceTwoSets same as the one above but where there are two sets
func PConcordanceTwoSets(comp []Bipart, bps []Bipart, jobs <-chan []int, results chan<- []int) {
	for j := range jobs {
		in1, in2 := j[0], j[1]
		b := 0
		// there must be at least one difference in the Trees so it isn't just the same tree
		if comp[in1].ConcordantWith(bps[in2]) {
			b = 1
		}
		results <- []int{in1, in2, b}
	}
}

// OutputEdges just print the edges
// mapints are int to string names for the taxa
// bps list of biparts
// ntrees number of trees
func OutputEdges(mapints map[int]string, bps []Bipart, ntrees int) {
	//sorted
	nn := map[int][]int{}
	var sortedCounts []int
	for v := range bps {
		nn[bps[v].Ct] = append(nn[bps[v].Ct], v)
	}
	for k := range nn {
		sortedCounts = append(sortedCounts, k)
	}
	sort.Sort(sort.Reverse(sort.IntSlice(sortedCounts)))
	var sortedBps []int
	for _, m := range sortedCounts {
		for _, k := range nn[m] {
			sortedBps = append(sortedBps, k)
		}
	}
	fmt.Println("numintrees percintrees bipart")
	for _, x := range sortedBps {
		i := x
		b := bps[x]
		fmt.Println(len(bps[i].TreeIndices), float64(len(bps[i].TreeIndices))/float64(ntrees), b.NewickWithNames(mapints))
	}
}

// CompareTreeToBiparts take biparts from a set , comparetreebps, and compre them to another set bps
// this one is complicated so keep with it
func CompareTreeToBiparts(bps []Bipart, comptreebps []Bipart, workers int, mapints map[int]string, verbose bool) {
	jobs := make(chan []int, len(bps)*len(comptreebps))
	results := make(chan []int, len(bps)*len(comptreebps))
	for w := 1; w <= workers; w++ {
		go PConflictsCompTree(bps, comptreebps, jobs, results)
	}
	njobs := 0
	for j := range comptreebps {
		for i := range bps {
			jobs <- []int{i, j}
			njobs++
		}
	}
	close(jobs)
	compconfs := make(map[int][]int)             // key is compbipart and value are the conflicts
	allconfs := make(map[int]bool)               // list of all the conflicting biparts from bps
	compconfstrees := make(map[int]map[int]bool) //key is compbipart and value are the conflicting trees
	for i := 0; i < njobs; i++ {
		x := <-results
		if x[2] == 1 {
			compconfs[x[1]] = append(compconfs[x[1]], x[0])
			allconfs[x[0]] = true
			if _, ok := compconfstrees[x[1]]; !ok {
				compconfstrees[x[1]] = make(map[int]bool)
			}
			for _, m := range bps[x[0]].TreeIndices {
				compconfstrees[x[1]][m] = true
			}
		}
	}
	for x := range compconfs {
		for _, n := range comptreebps[x].Nds {
			n.SData["conf"] = strconv.Itoa(len(compconfstrees[x]))
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
	njobs = 0
	for i := range allconfs {
		for j := range allconfs {
			jobs <- []int{i, j}
			njobs++
		}
	}
	close(jobs)
	bpsConcCounts := make(map[int]int)             // key bipart index, value number of concordant bps
	bpsConcTrees := make(map[int]map[int]bool)     // key bipart index, value is list of concordant tree
	compbpsConcTrees := make(map[int]map[int]bool) // key bipart index, value is list of concordant tree
	for i := 0; i < njobs; i++ {
		x := <-results
		if x[2] == 1 {
			if _, ok := bpsConcCounts[x[0]]; ok {
				bpsConcCounts[x[0]]++

			} else {
				bpsConcCounts[x[0]] = 0
				bpsConcTrees[x[0]] = make(map[int]bool)
			}
			if _, ok := bpsConcCounts[x[1]]; ok {
				bpsConcCounts[x[1]]++
			} else {
				bpsConcCounts[x[1]] = 0
				bpsConcTrees[x[1]] = make(map[int]bool)
			}
			if verbose {
				for _, m := range bps[x[0]].TreeIndices {
					bpsConcTrees[x[0]][m] = true
					bpsConcTrees[x[1]][m] = true
				}
				for _, m := range bps[x[1]].TreeIndices {
					bpsConcTrees[x[0]][m] = true
					bpsConcTrees[x[1]][m] = true
				}
			}
		}
	}

	// verbose comp concordance with bps
	if verbose {
		jobs = make(chan []int, len(comptreebps)*len(bps))
		results = make(chan []int, len(comptreebps)*len(bps))

		for w := 1; w <= workers; w++ {
			go PConcordanceTwoSets(comptreebps, bps, jobs, results)
		}
		njobs = 0
		for i := range comptreebps {
			for j := range bps {
				jobs <- []int{i, j}
				njobs++
			}
		}
		close(jobs)
		for i := 0; i < njobs; i++ {
			x := <-results // x[0] is compbpsindex, x[1] is bpsindex
			if x[2] == 1 {
				if _, ok := compbpsConcTrees[x[0]]; !ok {
					compbpsConcTrees[x[0]] = make(map[int]bool)
				}
				for _, m := range bps[x[1]].TreeIndices {
					compbpsConcTrees[x[0]][m] = true
				}
			}
		}
	}
	/*
	   sorting the results so that the larger bps are listed first. stop printing after a few.
	   add a sys command for listing all the results
	*/
	minout := 100
	// add things that don't conflict so that we can get concordance
	if verbose {
		for x := range comptreebps {
			if _, ok := compconfs[x]; !ok {
				fmt.Print("(", comptreebps[x].Index, ") ", comptreebps[x].NewickWithNames(mapints)+"\n")
				fmt.Print("  conctrees [", len(compbpsConcTrees[x]), "]: ", IntMapSetString(compbpsConcTrees[x])+"\n")
			}
		}
	}
	for x, y := range compconfs {
		fmt.Print("(", comptreebps[x].Index, ") ", comptreebps[x].NewickWithNames(mapints)+"\n")
		if verbose {
			fmt.Print("  conctrees [", len(compbpsConcTrees[x]), "]: ", IntMapSetString(compbpsConcTrees[x])+"\n")
			fmt.Print("  conftrees [", len(compconfstrees[x]), "]: ", IntMapSetString(compconfstrees[x])+"\n")
		}
		// put the number of conc at the internal nodes
		for _, n := range comptreebps[x].Nds {
			n.SData["conc"] = strconv.Itoa(len(compbpsConcTrees[x]))
		}
		n := map[int][]int{}
		var a []int
		for _, v := range y {
			//n[bpsCounts[v]] = append(n[bpsCounts[v]], v)
			n[bpsConcCounts[v]] = append(n[bpsConcCounts[v]], v)
			if verbose {
				if _, ok := bpsConcTrees[v]; !ok {
					bpsConcTrees[v] = make(map[int]bool)
				}
				for _, m := range bps[v].TreeIndices {
					bpsConcTrees[v][m] = true
				}
			}
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
				fmt.Print("  ", "(", bps[s].Index, ") ", bps[s].Ct, " ", len(bpsConcTrees[s]), " ", bpsConcCounts[s], " "+bps[s].NewickWithNames(mapints)+"\n")
				if verbose {
					fmt.Print("    trees [", len(bpsConcTrees[s]), "]:", IntMapSetString(bpsConcTrees[s]), "\n")
				}
				if count >= 10 {
					break
				}
				count++
			}
			if count >= minout {
				break
			}
		}
	}
}
