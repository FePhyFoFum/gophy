package gophy

import (
	"fmt"
	"math"
	"os"
	"strconv"

	"gonum.org/v1/gonum/mat"
)

// MultStateModel multistate model struct
type MultStateModel struct {
	BF        []float64 // base frequencies
	EBF       []float64 //empirical base freqs
	R         *mat.Dense
	Q         *mat.Dense // common use
	CharMap   map[string][]int
	NumStates int
	//sync.RWMutex
	Ps map[float64]*mat.Dense
	X  *mat.Dense
	P  mat.Dense
	//for decomposing
	QS         *mat.Dense
	EigenVals  []float64  // to be exponentiated
	EigenVecs  *mat.Dense //
	EigenVecsI *mat.Dense
	X1         *mat.Dense
	X2         *mat.Dense
}

// NewMultStateModel get new MULTModel pointer
func NewMultStateModel() *MultStateModel {
	return &MultStateModel{}
}

//SetEqualBF set all equal
func (d *MultStateModel) SetEqualBF() {
	d.BF = make([]float64, d.NumStates)
	for i := 0; i < d.NumStates; i++ {
		d.BF[i] = 1. / float64(d.NumStates)
	}
}

//SetEmpiricalBF set all to empirical
func (d *MultStateModel) SetEmpiricalBF() {
	d.BF = d.EBF
}

//SetRateMatrix length is (((numstates*numstates)-numstates)/2) - 1
// or (numstates * numstates) - numstates
// this is for scaled branch lengths and matrices
// CHANGE THIS TO UNSCALED
func (d *MultStateModel) SetRateMatrix(params []float64) {
	d.R = mat.NewDense(d.NumStates, d.NumStates, nil)
	if len(params) == (((d.NumStates*d.NumStates)-d.NumStates)/2)-1 {
		//symm
		pcount := 0
		for i := 0; i < d.NumStates; i++ {
			for j := 0; j < d.NumStates; j++ {
				if j > i {
					continue
				}
				if i == j {
					d.R.Set(i, j, 0.0)
				} else if i == d.NumStates-1 && j == d.NumStates-2 {
					d.R.Set(i, j, 1.0)
					d.R.Set(j, i, 1.0)
				} else {
					d.R.Set(i, j, params[pcount])
					d.R.Set(j, i, params[pcount])
					pcount++
				}
			}
		}
	} else if len(params) == ((d.NumStates*d.NumStates)-d.NumStates)-1 {
		//nonsymm
		pcount := 0
		for i := 0; i < d.NumStates; i++ {
			for j := 0; j < d.NumStates; j++ {
				if i == j {
					d.R.Set(i, j, 0.0)
				} else if i == d.NumStates-1 && j == d.NumStates-2 {
					d.R.Set(i, j, 1.0)
				} else {
					d.R.Set(i, j, params[pcount])
					pcount++
				}
			}
		}
	} else {
		fmt.Println("WRONG MATRIX SIZE")
		os.Exit(1)
	}
}

// SetBaseFreqs needs to be done before doing other things
func (d *MultStateModel) SetBaseFreqs(basefreq []float64) {
	d.BF = basefreq
}

// SetupQJC setup Q matrix
//    This is scaled so that change is reflected in the branch lengths
//    You don't need to use the SetScaledRateMatrix
func (d *MultStateModel) SetupQJC() {
	d.Ps = make(map[float64]*mat.Dense)
	d.Q = mat.NewDense(d.NumStates, d.NumStates, nil)

	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			if i != j {
				d.Q.Set(i, j, 1./float64(d.NumStates-1))
			} else {
				d.Q.Set(i, j, -1.)
			}
		}
	}

}

// SetupQJC1Rate setup Q matrix with one rate, probably for anc multi state
//     These are unscaled so the branch lengths are going to be time or something else
//     and not relative to these rates
//     Will take BF from something else
func (d *MultStateModel) SetupQJC1Rate(rt float64) {
	d.Ps = make(map[float64]*mat.Dense)
	d.Q = mat.NewDense(d.NumStates, d.NumStates, nil)

	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			if i != j {
				d.Q.Set(i, j, rt)
			} else {
				d.Q.Set(i, j, -(float64(d.NumStates-1) * rt))
			}
		}
	}
}

// SetupQMk setup Q matrix
//    This is unscaled (so the branch lengths are going to be proportion to some other change
//    and not to these branch lengths)
//    Will take the BF from something else
func (d *MultStateModel) SetupQMk(rt []float64, sym bool) {
	d.Ps = make(map[float64]*mat.Dense)
	d.Q = mat.NewDense(d.NumStates, d.NumStates, nil)
	cc := 0
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			if i != j {
				if sym && j > i {
					d.Q.Set(i, j, rt[cc])
					d.Q.Set(j, i, rt[cc])
					cc++
				} else if sym == false {
					d.Q.Set(i, j, rt[cc])
					cc++
				}
			} else {
				d.Q.Set(i, j, 0.0)
			}
		}
	}
	for i := 0; i < d.NumStates; i++ {
		sumrow := 0.
		for j := 0; j < d.NumStates; j++ {
			sumrow += d.Q.At(i, j)
		}
		d.Q.Set(i, i, -sumrow)
	}
}

// SetupQGTR setup Q matrix
//    This is scaled (so the branch lengths are going to be proportional to these changes)
//    Use SetScaledRateMatrix and then do this
func (d *MultStateModel) SetupQGTR() {
	fl := make([]float64, d.NumStates*d.NumStates)
	for i := 0; i < d.NumStates*d.NumStates; i++ {
		fl[i] = 1.0
	}
	bigpi := mat.NewDense(d.NumStates, d.NumStates, fl)
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			if i != j {
				bigpi.Set(i, j, d.BF[i]*d.BF[j])
			} else {
				bigpi.Set(i, j, d.BF[i])
			}
		}
	}
	dQ := mat.NewDense(d.NumStates, d.NumStates, nil)
	dQ.MulElem(d.R, bigpi)
	for i := 0; i < d.NumStates; i++ {
		dQ.Set(i, i, 0.0)
	}
	s := sumMatrix(dQ)
	dQ.Scale(1/s, dQ)
	for i := 0; i < d.NumStates; i++ {
		dQ.Set(i, i, 0-sumRow(dQ, i))
	}
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			dQ.Set(i, j, dQ.At(i, j)/d.BF[j])
		}
	}
	m := dQ.T()
	d.Q = mat.NewDense(d.NumStates, d.NumStates, nil)
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			d.Q.Set(i, j, m.At(i, j))
		}
	}
}

// DecomposeQ this is just for NR optimization for branch lengths
func (d *MultStateModel) DecomposeQ() {
	d.QS = mat.NewDense(d.NumStates, d.NumStates, nil)
	d.EigenVecsI = mat.NewDense(d.NumStates, d.NumStates, nil)
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			d.QS.Set(i, j, d.Q.At(i, j))
		}
	}
	//decompose, each time you change the model
	var ES mat.Eigen
	d.EigenVecs = mat.NewDense(d.NumStates, d.NumStates, nil)
	ES.Factorize(d.QS, mat.EigenBoth) //true, true)
	TC := ES.VectorsTo(nil)
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			d.EigenVecs.Set(i, j, real(TC.At(i, j)))
		}
	}
	d.EigenVecsI.Inverse(d.EigenVecs)
	d.EigenVals = make([]float64, d.NumStates)
	for i := 0; i < d.NumStates; i++ {
		d.EigenVals[i] = 1.
	}
	TV := ES.Values(nil)
	for i := 0; i < d.NumStates; i++ {
		d.EigenVals[i] = real(TV[i])
	}

	d.X = mat.NewDense(d.NumStates, d.NumStates, nil)  // P
	d.X1 = mat.NewDense(d.NumStates, d.NumStates, nil) // first der
	d.X2 = mat.NewDense(d.NumStates, d.NumStates, nil) // second der
}

//SetScaledRateMatrix needs to be done before doing SetupQGTR
// just send along the rates and this will make them the whole matrix
// the scaled is that this is assuming that the last rate is 1
//THIS IS THE SAME? AS SETRATEMATRIX?
func (d *MultStateModel) SetScaledRateMatrix(params []float64, sym bool) {
	d.R = mat.NewDense(d.NumStates, d.NumStates, nil)
	cc := 0
	lasti := 0
	lastj := 0
	if sym {
		lasti = d.NumStates - 2
		lastj = d.NumStates - 1
	} else {
		lasti = d.NumStates - 1
		lastj = d.NumStates - 2
	}
	for i := 0; i < d.NumStates; i++ {
		d.R.Set(i, i, 0)
		for j := 0; j < d.NumStates; j++ {
			if i == j {
				continue
			}
			if sym && j > i {
				if j == lastj && i == lasti {
					d.R.Set(i, j, 1.0)
					d.R.Set(j, i, 1.0)
				} else {
					d.R.Set(i, j, params[cc])
					d.R.Set(j, i, params[cc])
					cc++
				}
			} else if sym == false {
				if j == lastj && i == lasti {
					d.R.Set(i, j, 1.0)
				} else {
					d.R.Set(i, j, params[cc])
					cc++
				}
			}
		}
	}
}

// ExpValue used for the matrix exponential
func (d *MultStateModel) ExpValue(iv []float64, blen float64) {
	for i, j := range iv {
		d.X.Set(i, i, math.Exp(j*blen))
	}
	return
}

// ExpValueFirstD get the first derivaite for NR
func (d *MultStateModel) ExpValueFirstD(blen float64) (x *mat.Dense) {
	x = mat.NewDense(d.NumStates, d.NumStates, nil)
	for i := 0; i < d.NumStates; i++ {
		d.X1.Set(i, i, d.EigenVals[i]*math.Exp(d.EigenVals[i]*blen))
	}
	x.Mul(d.EigenVecs, d.X1)
	x.Mul(x, d.EigenVecsI)
	return
}

// ExpValueSecondD get the second derivaite for NR
func (d *MultStateModel) ExpValueSecondD(blen float64) (x *mat.Dense) {
	x = mat.NewDense(d.NumStates, d.NumStates, nil)
	for i := 0; i < d.NumStates; i++ {
		d.X1.Set(i, i, (d.EigenVals[i]*d.EigenVals[i])*math.Exp(d.EigenVals[i]*blen))
	}
	x.Mul(d.EigenVecs, d.X1)
	x.Mul(x, d.EigenVecsI)
	return
}

// SetP use the standard spectral decom
func (d *MultStateModel) SetP(blen float64) {
	P := mat.NewDense(d.NumStates, d.NumStates, nil)
	d.ExpValue(d.EigenVals, blen)
	P.Mul(d.EigenVecs, d.X)
	P.Mul(P, d.EigenVecsI)
	// easier
	//d.Lock()
	d.Ps[blen] = P
	//d.Unlock()
}

// SetPSimple use the gonum matrixexp (seems faster)
func (d *MultStateModel) SetPSimple(blen float64) {
	P := mat.NewDense(d.NumStates, d.NumStates, nil)
	P.Scale(blen, d.Q)
	P.Exp(P)
	//d.Lock()
	d.Ps[blen] = P
	//d.Unlock()
}

// EmptyPDict save memory perhaps?
func (d *MultStateModel) EmptyPDict() {
	d.Ps = nil
	d.Ps = make(map[float64]*mat.Dense)
}

// GetPMap get the Ps from the dictionary
func (d *MultStateModel) GetPMap(blen float64) *mat.Dense {
	//d.RLock()
	if _, ok := d.Ps[blen]; ok {
		//d.RLock()
		return d.Ps[blen]
		//X := d.Ps[blen]
		//	d.RUnlock()
		//return X
	}
	//d.RUnlock()
	d.SetPSimple(blen)
	//d.RLock()
	//X := d.Ps[blen]
	//d.RUnlock()
	return d.Ps[blen]
}

// GetPCalc calculate P matrix
func (d *MultStateModel) GetPCalc(blen float64) *mat.Dense {
	var P mat.Dense
	P.Scale(blen, d.Q)
	P.Exp(&P)
	//d.ExpValue(d.EigenVals, blen)
	//P.Mul(&d.EigenVecs, d.X)
	//P.Mul(&P, d.EigenVecsT)
	return &P
}

//SetMap for getting the position in the array
func (d *MultStateModel) SetMap() {
	d.CharMap = make(map[string][]int)
	d.CharMap["-"] = make([]int, d.NumStates)
	d.CharMap["N"] = make([]int, d.NumStates)
	for i := 0; i < d.NumStates; i++ {
		d.CharMap[strconv.Itoa(i)] = []int{i}
		d.CharMap["-"][i] = i
		d.CharMap["N"][i] = i
	}
}

func (d *MultStateModel) GetCharMap() map[string][]int {
	return d.CharMap
}

func (d *MultStateModel) GetNumStates() int {
	return d.NumStates
}

func (d *MultStateModel) GetBF() []float64 {
	return d.BF
}

func (d *MultStateModel) GetStochMapMatrices(dur float64, from int, to int) (summed *mat.Dense, summedR *mat.Dense) {
	nstates := d.NumStates
	d.DecomposeQ()
	Ql := mat.NewDense(nstates, nstates, nil)
	Ql.Zero()
	Ql.Set(from, to, d.Q.At(from, to))
	W := mat.NewDense(nstates, nstates, nil)
	W.Zero()
	W.Set(from, from, 1.)
	summed = mat.NewDense(nstates, nstates, nil)
	summed.Zero()
	summedR = mat.NewDense(nstates, nstates, nil)
	summedR.Zero()
	for i := 0; i < nstates; i++ {
		Ei := mat.NewDense(nstates, nstates, nil)
		Ei.Zero()
		Ei.Set(i, i, 1.)
		Si := mat.NewDense(nstates, nstates, nil)
		Si.Mul(d.EigenVecs, Ei)
		Si.Mul(Si, d.EigenVecsI)
		for j := 0; j < nstates; j++ {
			dij := (d.EigenVals[i] - d.EigenVals[j]) * dur
			Ej := mat.NewDense(nstates, nstates, nil)
			Ej.Zero()
			Ej.Set(j, j, 1.)
			Sj := mat.NewDense(nstates, nstates, nil)
			Sj.Mul(d.EigenVecs, Ej)
			Sj.Mul(Sj, d.EigenVecsI)
			Iijt := 0.
			if math.Abs(dij) > 10 {
				Iijt = (math.Exp(d.EigenVals[i]*dur) - math.Exp(d.EigenVals[j]*dur)) / (d.EigenVals[i] - d.EigenVals[j])
			} else if math.Abs(dij) < 10e-20 {
				Iijt = dur * math.Exp(d.EigenVals[j]*dur) * (1. + dij/2. + math.Pow(dij, 2.)/6. + math.Pow(dij, 3.)/24.)
			} else {
				if d.EigenVals[i] == d.EigenVals[j] {
					//					if isImag {
					//						Iijt = dur * exp(d.EigenVals[j]*dur) * (exp(dij) - 1.) / dij
					//					} else {
					Iijt = dur * math.Exp(d.EigenVals[j]*dur) * (math.Expm1((dij))) / dij //(expm1(real(dij))) / dij
					//					}
				} else {
					//					if isImag {
					//						Iijt = -dur * exp(d.EigenVals[i]*dur) * (exp(-dij) - 1.) / dij
					//					} else {
					Iijt = -dur * math.Exp(d.EigenVals[i]*dur) * (math.Expm1((-dij))) / dij //(expm1(real(-dij))) / dij
					//					}
				}
			}
			temp := mat.NewDense(nstates, nstates, nil)
			temp.Mul(Si, Ql)
			temp.Mul(temp, Sj)
			temp.Scale(Iijt, temp)
			summed.Add(summed, temp)
			temp = mat.NewDense(nstates, nstates, nil)
			temp.Mul(Si, W)
			temp.Mul(temp, Sj)
			temp.Scale(Iijt, temp)
			summedR.Add(summedR, temp)
		}
	}
	return
}

func PrintMatrix(d *mat.Dense, diag bool) {
	r, c := d.Dims()
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			fmt.Print("  ")
			if diag && i == j {
				fmt.Print("- ")
			} else {
				fmt.Print(d.At(i, j), " ")
			}
		}
		fmt.Print("\n")
	}
}
