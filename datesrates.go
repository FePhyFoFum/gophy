package gophy

import "math"

//PLCalcLogLike penalized likelihood calculator
func PLCalcLogLike(rates []float64, durations []float64, chardurations []float64,
	lgFacCharDurations []float64) (ll float64) {
	for i := range rates {
		x := rates[i] * durations[i]
		l := -(chardurations[i]*math.Log(x) - x - lgFacCharDurations[i])
		ll += l
	}
	return
}
