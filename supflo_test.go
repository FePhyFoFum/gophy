package gophy_test

import (
	"fmt"
	"math"
	"testing"

	"github.com/FePhyFoFum/gophy"
)

func TestGetLn(t *testing.T) {
	x := gophy.NewSupFlo(2.2, 1)
	if x.GetLn().Float64() != math.Log(22) {
		t.Fail()
	}
	fmt.Println(x.GetLn())
	fmt.Println(math.Log(22))
}

func TestFloat64(t *testing.T) {
	x := gophy.NewSupFlo(2.2, 1)
	if x.Float64() != 22 {
		t.Fail()
	}
	fmt.Println(x, x.Float64())
}

func TestMulEqFloat(t *testing.T) {
	x := gophy.NewSupFlo(2.2, 1)
	x.MulEqFloat(0.2)
	if x.Float64() != (22 * 0.2) {
		t.Fail()
	}
	fmt.Println(x, x.Float64(), 22*0.2)
}

func TestAdjust(t *testing.T) {
	x := gophy.NewSupFlo(2.2*10e100, 0)
	x.MulEqFloat(0.2)
	if x.Float64() != (4.4e+100) {
		t.Fail()
	}
	fmt.Println(x, x.Float64(), 4.4e+100)
}
