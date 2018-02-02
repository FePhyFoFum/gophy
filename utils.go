package gophy

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
