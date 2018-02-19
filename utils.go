package gophy

import "math"

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

//SliceStringContains tells you whether the e string is in the slice
func SliceStringContains(s []string, e string) bool {
	for _, a := range s {
		if a == e {
			return true
		}
	}
	return false
}

func SumFloatVec(x []float64) (s float64) {
	for _, a := range x {
		s += a
	}
	return
}

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
