package gophy_test

import (
	"fmt"
	"math"
	"sort"
	"testing"

	"github.com/FePhyFoFum/gophy"
)

// helper function to compare two int slices
func Equal(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i, v := range a {
		if v != b[i] {
			return false
		}
	}
	return true
}

func TestLogFactorial(t *testing.T) {
	if gophy.LogFactorial(10) != 15.104412573075518 {
		fmt.Println(gophy.LogFactorial(10))
		t.Fail()
	}
}

func TestCalcSliceIntDifferenceInt(t *testing.T) {
	a := []int{0, 1, 2}
	b := []int{2, 3, 4}
	x := gophy.CalcSliceIntDifferenceInt(a, b)
	y := gophy.CalcSliceIntDifferenceInt(b, a)
	if x != 2 || y != 2 {
		fmt.Println(x)
		fmt.Println(y)
		t.Fail()
	}
}

func TestCalcSlicIntDifference(t *testing.T) {
	a := []int{0, 1, 2}
	b := []int{2, 3, 4}
	x := gophy.CalcSliceIntDifference(a, b)
	y := gophy.CalcSliceIntDifference(b, a)
	if !Equal([]int{0, 1}, x) || !Equal([]int{3, 4}, y) {
		fmt.Println(x)
		fmt.Println(y)
		t.Fail()
	}
}

// TestPCalcSliceIntDifferenceInt goes here

func TestIntSliceContains(t *testing.T) {
	a := 1
	b := []int{0, 1, 2}
	if !(gophy.IntSliceContains(b, a)) {
		fmt.Println(gophy.IntSliceContains(b, a))
		t.Fail()
	}
}

func TestIntMapIntersects(t *testing.T) {
	a := map[int]bool{0: true, 1: false, 2: false}
	b := map[int]bool{0: false, 1: false, 2: true}
	if !(gophy.IntMapIntersects(a, b)) {
		fmt.Println(gophy.IntMapIntersects(a, b))
		t.Fail()
	}
}

func TestIntMapIntersects2(t *testing.T) {
	a := map[int]bool{0: true, 1: true, 2: false}
	b := map[int]bool{0: true, 1: false, 2: true}
	if !(gophy.IntMapIntersects2(a, b)) {
		fmt.Println(gophy.IntMapIntersects2(a, b))
		t.Fail()
	}
}

func TestIntMapIntersectsRet(t *testing.T) {
	a := map[int]bool{0: true, 1: false, 2: false}
	b := map[int]bool{2: true, 3: true, 4: true}
	if !(Equal(gophy.IntMapIntersectsRet(a, b), []int{2})) {
		fmt.Println(gophy.IntMapIntersectsRet(a, b))
		t.Fail()
	}
}

func TestIntMapDifferenceRet(t *testing.T) {
	a := map[int]bool{0: true, 1: false, 2: false}
	b := map[int]bool{2: true, 3: true, 4: true}
	if !(Equal(gophy.IntMapDifferenceRet(a, b), []int{0, 1})) || !(Equal(gophy.IntMapDifferenceRet(b, a), []int{3, 4})) {
		fmt.Println(gophy.IntMapDifferenceRet(a, b))
		fmt.Println(gophy.IntMapDifferenceRet(b, a))
		t.Fail()
	}
}

func TestIntMapSetString(t *testing.T) {
	a := map[int]bool{0: true, 1: true, 2: true}
	if !(gophy.IntMapSetString(a) == "0 1 2 ") {
		fmt.Println(gophy.IntMapSetString(a))
		t.Fail()
	}
}

func TestStringSliceContains(t *testing.T) {
	e := "test"
	s := []string{"this", "is", "a", "test"}
	if !(gophy.StringSliceContains(s, e)) {
		fmt.Println(gophy.StringSliceContains(s, e))
		t.Fail()
	}
}

func TestNodeSliceContains(t *testing.T) {
	tree := gophy.ReadTreeFromFile("test_files/10tips.tre")
	n := tree.Rt
	a := tree.Post
	if !(gophy.NodeSliceContains(a, n)) {
		fmt.Println(gophy.NodeSliceContains(a, n))
		t.Fail()
	}
}

func TestSumFloatVec(t *testing.T) {
	a := []float64{1, 2, 3, 4, 5.5}
	if !(gophy.SumFloatVec(a) == float64(15.5)) {
		fmt.Println(gophy.SumFloatVec(a))
		t.Fail()
	}
}

func TestLog1exp(t *testing.T) {
	var a float64 = 36.0
	var b float64 = -11.0
	var c float64 = 2.0
	if !(gophy.Log1exp(a) == 36.0) {
		fmt.Println(gophy.Log1exp(a))
		t.Fail()
	}
	if !(gophy.Log1exp(b) == math.Exp(-11.0)) {
		fmt.Println(gophy.Log1exp(b))
		fmt.Println(math.Exp(-11.0))
		t.Fail()
	}
	if !(gophy.Log1exp(c) == 2.1269280110429727) {
		fmt.Println(gophy.Log1exp(c))
		t.Fail()
	}
}

func TestSumLogExp(t *testing.T) {
	a := 5.0
	b := 2.0
	if !(gophy.SumLogExp(a, b) == 5.0485873515737421) {
		fmt.Println(gophy.SumLogExp(a, b))
		t.Fail()
	}
}

func TestCalcAIC(t *testing.T) {
	ln := -4570.7964
	k := float64(20)
	if !(gophy.CalcAIC(ln, k) == 9181.5928) {
		//fmt.Printf("%.20f", ln)
		fmt.Println(gophy.CalcAIC(ln, k))
		t.Fail()
	}
}

func TestCalcAICC(t *testing.T) {
	ln := -4570.7964
	k := float64(20)
	n := 500
	if !(gophy.CalcAICC(ln, k, n) == 9183.3464534446775) {
		fmt.Println(gophy.CalcAICC(ln, k, n))
		t.Fail()
	}
}

func TestCalcBIC(t *testing.T) {
	ln := -4570.7964
	k := float64(20)
	n := 500
	if !(gophy.CalcBIC(ln, k, n) == 9265.8849619684443) {
		fmt.Println(gophy.CalcBIC(ln, k, n))
		t.Fail()
	}
}

func TestNodeSlicePosition(t *testing.T) {
	tree := gophy.ReadTreeFromFile("test_files/10tips.tre")
	n := tree.Rt
	a := tree.Pre
	b := tree.Post
	if !(gophy.NodeSlicePosition(a, n) == 0) || !(gophy.NodeSlicePosition(b, n) == 18) {
		fmt.Println(gophy.NodeSlicePosition(a, n))
		t.Fail()
	}
}

func TestRound(t *testing.T) {
	a := 1.234567
	if !(gophy.Round(a, 0.5, 5) == 1.23457) {
		fmt.Println(gophy.Round(a, 0.5, 5))
		t.Fail()
	}
}

func TestNewSortedIdxSliceD(t *testing.T) {
	s := gophy.NewSortedIdxSliceD(1, 25, 3, 5, 4)
	sort.Sort(s)
	if !(Equal(s.IntSlice, []int{1, 3, 4, 5, 25})) {
		fmt.Println(s.IntSlice)
		t.Fail()
	}
	if !(Equal(s.Idx, []int{0, 2, 4, 3, 1})) {
		fmt.Println(s.Idx)
		t.Fail()
	}
}

func TestNewSortedIdxSlice(t *testing.T) {
	s := gophy.NewSortedIdxSlice([]int{1, 25, 3, 5, 4})
	sort.Sort(s)
	if !(Equal(s.IntSlice, []int{1, 3, 4, 5, 25})) {
		fmt.Println(s.IntSlice)
		t.Fail()
	}
	if !(Equal(s.Idx, []int{0, 2, 4, 3, 1})) {
		fmt.Println(s.Idx)
		t.Fail()
	}
}

func TestLogFact(t *testing.T) {
	if !(gophy.LogFact(10) == 15.096081587310044) {
		fmt.Println(gophy.LogFact(10))
		t.Fail()
	}
}
func TestCalcConfIntLgLike(t *testing.T) {
	x := math.Round(gophy.CalcConfIntLgLike(2, 0.95)*1000) / 1000
	y := math.Round(math.Log(0.95)*1000) / 1000
	if !(x == y) {
		fmt.Println(x, y)
		t.Fail()
	}
}

func TestCalcNormPDF(t *testing.T) {
	x := gophy.CalcNormPDF(1.000, 0.1, 0.1)
	fmt.Println(x)
	t.Fail()
}

func TestCalcLogNormPDF(t *testing.T) {
	x := gophy.CalcLogNormPDF(1.0, 0.25)
	fmt.Println(x)
	t.Fail()
}

func TestCalcLogNormPDFLog(t *testing.T) {
	x := gophy.CalcLogNormPDFLog(1.900, 0.0, 0.25)
	fmt.Println(math.Exp(x), gophy.CalcLogNormPDF(1.900, 0.25))
	t.Fail()
}

func TestCalcLogNormLogScalePDF(t *testing.T) {
	x := gophy.CalcLogNormLocScalePDF(1.0, 0.25, 0.3, 1.)
	fmt.Println(x)
	t.Fail()
}
