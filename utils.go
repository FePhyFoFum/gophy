package gophy

import (
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/stat/distuv"

	"gonum.org/v1/gonum/stat"
)

func LogFactorial(val int) (x float64) {
	x = 0.
	for i := 1; i <= val; i++ {
		x += math.Log(float64(i))
	}
	return
}

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

//Rfwresult struct for holding info from rfwp analysis
type Rfwresult struct {
	Tree1  int
	Tree2  int
	Weight float64
	MaxDev float64
}

//PCalcRFDistancesPartialWeighted includes the branch lengths
func PCalcRFDistancesPartialWeighted(bpts map[int][]int, tippenalty bool, bps []Bipart, jobs <-chan []int, results chan<- Rfwresult) {
	for j := range jobs {
		in1, in2 := j[0], j[1]
		ab := 0.0
		mb := map[string]float64{}
		used := make([]int, 0)
		maxdev := 0.0
		for _, x := range bpts[in1] {
			for _, y := range bpts[in2] {
				if bps[x].ConflictsWith(bps[y]) {
					mb["t1"+string(x)] = bps[x].NdsM[in1].Len
					mb["t2"+string(y)] = bps[y].NdsM[in2].Len
					// record the max dev
					if bps[x].NdsM[in1].Len > maxdev {
						maxdev = bps[x].NdsM[in1].Len
					}
					if bps[y].NdsM[in2].Len > maxdev {
						maxdev = bps[y].NdsM[in2].Len
					}
					if tippenalty {
						used = append(used, x)
						used = append(used, y)
					}
				} else if bps[x].ConcordantWith(bps[y]) {
					v := math.Abs(bps[x].NdsM[in1].Len - bps[y].NdsM[in2].Len)
					ab += v
					if v > maxdev {
						maxdev = v
					}
					if tippenalty {
						used = append(used, x)
						used = append(used, y)
					}
				} else if bps[x].Equals(bps[y]) {
					v := math.Abs(bps[x].NdsM[in1].Len - bps[y].NdsM[in2].Len)
					ab += v
					if v > maxdev {
						maxdev = v
					}
					if tippenalty {
						used = append(used, x)
						used = append(used, y)
					}
				}
			}
		}
		// just adding a penalty for all the taxa that are missing.
		if tippenalty {
			x := CalcSliceIntDifference(bpts[in1], used)
			y := CalcSliceIntDifference(bpts[in2], used)
			for _, m := range x {
				if len(bps[m].Lt) == 1 {
					ab += bps[m].NdsM[in1].Len
					if bps[m].NdsM[in1].Len > maxdev {
						maxdev = bps[m].NdsM[in1].Len
					}
				} else {
					ab += bps[m].NdsM[in1].Len
					if bps[m].NdsM[in1].Len > maxdev {
						maxdev = bps[m].NdsM[in1].Len
					}
				}
			}
			for _, m := range y {
				if len(bps[m].Lt) == 1 {
					ab += bps[m].NdsM[in2].Len
					if bps[m].NdsM[in2].Len > maxdev {
						maxdev = bps[m].NdsM[in2].Len
					}
				} else {
					ab += bps[m].NdsM[in2].Len
					if bps[m].NdsM[in2].Len > maxdev {
						maxdev = bps[m].NdsM[in2].Len
					}
				}
			}
		}
		for _, x := range mb {
			ab += x
		}
		results <- Rfwresult{Tree1: in1, Tree2: in2, Weight: ab, MaxDev: maxdev}
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

// IntMapIntersectsRet checks to see if the two map[int]bool intersect and returns the intersection (in the set sense)
func IntMapIntersectsRet(s1, s2 map[int]bool) (r []int) {
	for k := range s1 {
		if s2[k] {
			r = append(r, k)
		}
	}
	return
}

// IntMapDifferenceRet calculate the difference (set) between two int slices
func IntMapDifferenceRet(a, b map[int]bool) []int {
	mb := map[int]bool{}
	for x := range b {
		mb[x] = true
	}
	ab := []int{}
	for x := range a {
		if _, ok := mb[x]; !ok {
			ab = append(ab, x)
		}
	}
	return ab
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

// NodeSliceContains tells you whether the e string is in the slice
func NodeSliceContains(s []*Node, e *Node) bool {
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

// NodeSlicePosition take a *[]Node slice and teh get the index of the element node
func NodeSlicePosition(sl []*Node, nd *Node) (x int) {
	x = -1
	for p, v := range sl {
		if v == nd {
			x = p
			return
		}
	}
	return
}

// Round to the nearest place probably val=num, roundOn = 0.5 places = 5
func Round(val float64, roundOn float64, places int) (newVal float64) {
	var round float64
	pow := math.Pow(10, float64(places))
	digit := pow * val
	_, div := math.Modf(digit)
	if div >= roundOn {
		round = math.Ceil(digit)
	} else {
		round = math.Floor(digit)
	}
	newVal = round / pow
	return
}

/*
  Bits below are for getting sorted indices from a list
*/

//SortedIntIdxSlice for sorting indices of ints
type SortedIntIdxSlice struct {
	sort.IntSlice
	Idx []int
}

//Swap for sorting indices
func (s SortedIntIdxSlice) Swap(i, j int) {
	s.IntSlice.Swap(i, j)
	s.Idx[i], s.Idx[j] = s.Idx[j], s.Idx[i]
}

// NewSortedIdxSliceD usage
/* s := NewSlice(1, 25, 3, 5, 4)
sort.Sort(s)
will give s.IntSlice = [1 3 4 5 25]
s.idx = [0 2 4 3 1]*/
func NewSortedIdxSliceD(n ...int) *SortedIntIdxSlice {
	s := &SortedIntIdxSlice{IntSlice: sort.IntSlice(n), Idx: make([]int, len(n))}
	for i := range s.Idx {
		s.Idx[i] = i
	}
	return s
}

// NewSortedIdxSlice usage
func NewSortedIdxSlice(n []int) *SortedIntIdxSlice {
	s := &SortedIntIdxSlice{IntSlice: sort.IntSlice(n), Idx: make([]int, len(n))}
	for i := range s.Idx {
		s.Idx[i] = i
	}
	return s
}

// LogFact calculate the log factorial
func LogFact(k float64) float64 {
	return k*(math.Log(k)-1) + math.Log(math.Sqrt(6.28318*k))
}

// MedianF calculate the "median" value
func MedianF(n []float64) float64 {
	sf := sort.Float64Slice(n)
	sf.Sort()
	m := len(sf) / 2
	if len(sf)%2 != 0 { //odd
		return sf[m]
	}
	return (sf[m-1] + sf[m]) / 2
}

// MaxF max
func MaxF(n []float64) float64 {
	sf := sort.Float64Slice(n)
	sf.Sort()
	return sf[len(sf)-1]
}

// MinF max
func MinF(n []float64) float64 {
	sf := sort.Float64Slice(n)
	sf.Sort()
	return sf[0]
}

// ConfInt95TF returns 95% conf int
// t-stat
func ConfInt95TF(nums []float64) (float64, float64) {
	mean := stat.Mean(nums, nil)
	SE := stat.StdDev(nums, nil) / math.Sqrt(float64(len(nums)))
	V := float64(len(nums) - 1)
	ST := distuv.StudentsT{Mu: 0, Nu: V, Sigma: 1}
	QL := mean + ST.Quantile(0.025)*SE
	QH := mean + ST.Quantile(0.975)*SE
	return QL, QH
}

//ReadLine is like the Python readline() and readlines()
func ReadLine(path string) (ln []string) {
	b, err := ioutil.ReadFile(path)
	if err != nil {
		fmt.Println(err)
		fmt.Println("There was an error when reading in the file:", path, ". Are you sure that it exists?")
		os.Exit(0)
	}
	ss := string(b)
	ln = strings.Split(ss, "\n")
	return
}
