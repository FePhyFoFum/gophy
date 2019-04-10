package gophy

import (
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/mathext"
	"gonum.org/v1/gonum/stat/distuv"
)

func sumMatrix(m *mat.Dense) (s float64) {
	s = 0
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			s += m.At(i, j)
		}
	}
	return
}
func sumRow(m *mat.Dense, i int) (s float64) {
	s = 0
	for j := 0; j < 4; j++ {
		s += m.At(i, j)
	}
	return
}

//LikeResult likelihood value and site
type LikeResult struct {
	value float64
	site  int
}

func pointGamma(p float64, a float64, b float64) float64 {
	c := distuv.ChiSquared{K: 2 * (a)}
	return ((c.Quantile(p)) / (2 * b))
}

func GetGammaCats(alpha float64, cats int, median bool) []float64 {
	K := float64(cats)
	a, b := alpha, alpha
	factor := a / b * K
	rK := make([]float64, cats)
	if median == true {
		gap05 := 1.0 / (2.0 * K)
		for i := 0; i < cats; i++ {
			rK[i] = pointGamma((float64(i)*2.0+1.)*gap05, a, b)
		}
		t := 0.
		for i := 0; i < cats; i++ {
			t += rK[i]
		}
		for i := 0; i < cats; i++ {
			rK[i] *= factor / t
		}
		return rK
	} else {
		freqK := make([]float64, cats)
		for i := 0; i < cats-1; i++ {
			freqK[i] = pointGamma((float64(i)+1.0)/K, a, b)
		}
		for i := 0; i < cats-1; i++ {
			freqK[i] = mathext.GammaIncReg(a+1, freqK[i]*b)
		}
		rK[0] = freqK[0] * factor
		rK[cats-1] = (1 - freqK[cats-2]) * factor
		for i := 1; i < cats-1; i++ {
			rK[i] = (freqK[i] - freqK[i-1]) * factor
		}
		return rK
	}
}
