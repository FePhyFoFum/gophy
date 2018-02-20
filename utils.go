package gophy

import (
	"math"
	"strconv"
)

// CalcSliceIntDifferenceInt calculate the size of the difference (set) between two int slices
func CalcSliceIntDifferenceInt(a, b []int) int {
	mb := map[int]bool{}
	for _, x := range b {
		mb[x] = true
	}
	ab := 0
	for _, x := range a {
		if _, ok := mb[x]; !ok {
			ab++
		}
	}
	return ab
}

// CalcSliceIntDifference calculate the difference (set) between two int slices
func CalcSliceIntDifference(a, b []int) []int {
	mb := map[int]bool{}
	for _, x := range b {
		mb[x] = true
	}
	ab := []int{}
	for _, x := range a {
		if _, ok := mb[x]; !ok {
			ab = append(ab, x)
		}
	}
	return ab
}

// PCalcSliceIntDifferenceInt calculate the size of the difference (set) between two int slices in parallel
// used for: RF distance where the bpts is the tree index -> bipart index list map
func PCalcSliceIntDifferenceInt(bpts map[int][]int, jobs <-chan []int, results chan<- []int) {
	for j := range jobs {
		in1, in2 := j[0], j[1]
		mb := map[int]bool{}
		for _, x := range bpts[in2] {
			mb[x] = true
		}
		ab := 0
		for _, x := range bpts[in1] {
			if _, ok := mb[x]; !ok {
				ab++
			}
		}
		results <- []int{in1, in2, ab}
	}
}

//PCalcRFDistancesPartial calculates the partial rf, bpts is the tree index, bipart list,
// bps is the list of biparts
func PCalcRFDistancesPartial(bpts map[int][]int, bps []Bipart, jobs <-chan []int, results chan<- []int) {
	for j := range jobs {
		in1, in2 := j[0], j[1]
		mb := map[string]bool{}
		for _, x := range bpts[in1] {
			for _, y := range bpts[in2] {
				if bps[x].ConflictsWith(bps[y]) {
					mb["t1"+string(x)] = true
					mb["t2"+string(y)] = true
				}
			}
		}
		ab := 0
		for _, x := range mb {
			if x == true {
				ab++
			}
		}
		results <- []int{in1, in2, ab}
	}
}

// IntSliceIntersects checks to see whether two int slices intersect
func IntSliceIntersects(a, b []int) (rb bool) {
	rb = false
	for _, k := range a {
		for _, l := range b {
			if k == l {
				rb = true
				return
			}
		}
	}
	return
}

// IntSliceContains checks to see if the int slice contains an int and returns the bool
func IntSliceContains(is []int, s int) (rb bool) {
	rb = false
	for _, a := range is {
		if a == s {
			rb = true
			return
		}
	}
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

// IntMapSetString get a string for printing off a set
func IntMapSetString(intmap map[int]bool) (s string) {
	s = ""
	for m := range intmap {
		s += strconv.Itoa(m) + " "
	}
	return
}

// StringSliceContains tells you whether the e string is in the slice
func StringSliceContains(s []string, e string) bool {
	for _, a := range s {
		if a == e {
			return true
		}
	}
	return false
}

// SumFloatVec sum the float vectors
func SumFloatVec(x []float64) (s float64) {
	for _, a := range x {
		s += a
	}
	return
}

// SumLogExp sum log of exps
func SumLogExp(a, b float64) float64 {
	return a + log1exp(b-a)
}

func log1exp(x float64) float64 {
	if x > 35 {
		return x
	}
	if x < -10 {
		return math.Exp(x)
	}
	return math.Log1p(math.Exp(x))
}
