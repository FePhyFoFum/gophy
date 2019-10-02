package gophy

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

// DNAModel standard DNA struct
type DNAModel struct {
	BF      []float64 // base frequencies
	R       *mat.Dense
	Q       *mat.Dense // common use
	CharMap map[string][]int
	//sync.RWMutex
	Ps  map[float64]*mat.Dense
	PsL map[float64]*mat.Dense
	X   *mat.Dense
	P   mat.Dense
	//for decomposing
	QS         *mat.Dense
	EigenVals  []float64  // to be exponentiated
	EigenVecs  *mat.Dense //
	EigenVecsI *mat.Dense
	X1         *mat.Dense
	X2         *mat.Dense
	NumStates  int
}

// NewDNAModel get new DNAModel pointer
func NewDNAModel() *DNAModel {
	d := &DNAModel{}
	d.NumStates = 4
	return d
}

// DeepCopyDNAModel ...
func (d *DNAModel) DeepCopyDNAModel() *DNAModel {
	outm := NewDNAModel()
	outm.BF = []float64{0.25, 0.25, 0.25, 0.25}
	copy(outm.BF, d.BF)
	outm.Q = mat.NewDense(4, 4, nil)
	outm.R = mat.NewDense(4, 4, nil)
	outm.R.Copy(d.R)
	outm.Q.Copy(d.Q)
	outm.SetMap()
	return outm
}

// SetupQJC setup Q matrix
//    This is scaled so that change is reflected in the branch lengths
//    You don't need to use the SetScaledRateMatrix
func (d *DNAModel) SetupQJC() {
	d.BF = []float64{0.25, 0.25, 0.25, 0.25}
	d.Ps = make(map[float64]*mat.Dense)
	d.Q = mat.NewDense(4, 4, nil)

	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if i != j {
				d.Q.Set(i, j, 0.333333333333)
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
func (d *DNAModel) SetupQJC1Rate(rt float64) {
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

// SetupQGTR setup Q matrix
func (d *DNAModel) SetupQGTR() {
	bigpi := mat.NewDense(4, 4, []float64{1, 1, 1, 1,
		1, 1, 1, 1,
		1, 1, 1, 1,
		1, 1, 1, 1})
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if i != j {
				bigpi.Set(i, j, d.BF[i]*d.BF[j])
			} else {
				bigpi.Set(i, j, d.BF[i])
			}
		}
	}
	dQ := mat.NewDense(4, 4, nil)
	dQ.MulElem(d.R, bigpi)
	dQ.Set(0, 0, 0.0)
	dQ.Set(1, 1, 0.0)
	dQ.Set(2, 2, 0.0)
	dQ.Set(3, 3, 0.0)
	s := sumMatrix(dQ)
	dQ.Scale(1/s, dQ)
	for i := 0; i < 4; i++ {
		dQ.Set(i, i, 0-sumRow(dQ, i))
	}
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			dQ.Set(i, j, dQ.At(i, j)/d.BF[j])
		}
	}
	m := dQ.T()
	d.Q = mat.NewDense(4, 4, nil)
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			d.Q.Set(i, j, m.At(i, j))
		}
	}
}

// SetupQMk setup Q matrix
//    This is unscaled (so the branch lengths are going to be proportion to some other change
//    and not to these branch lengths)
//    Will take the BF from something else
func (d *DNAModel) SetupQMk(rt []float64, sym bool) {
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

//SetScaledRateMatrix needs to be done before doing SetupQGTR
// just send along the rates and this will make them the whole matrix
// the scaled is that this is assuming that the last rate is 1
func (d *DNAModel) SetScaledRateMatrix(params []float64, sym bool) {
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

// DecomposeQ this is just for NR optimization for branch lengths
func (d *DNAModel) DecomposeQ() {
	d.QS = mat.NewDense(4, 4, nil)
	d.EigenVecsI = mat.NewDense(4, 4, nil)
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			d.QS.Set(i, j, d.Q.At(i, j))
		}
	}
	//decompose, each time you change the model
	var ES mat.Eigen
	d.EigenVecs = mat.NewDense(4, 4, nil)
	ES.Factorize(d.QS, mat.EigenBoth) //true, true)
	TC := mat.NewCDense(0, 0, nil)
	//	TC := ES.VectorsTo(nil)
	ES.VectorsTo(TC)
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			d.EigenVecs.Set(i, j, real(TC.At(i, j)))
		}
	}
	d.EigenVecsI.Inverse(d.EigenVecs)
	d.EigenVals = []float64{1., 1., 1., 1.}
	TV := ES.Values(nil)
	for i := 0; i < 4; i++ {
		d.EigenVals[i] = real(TV[i])
	}

	d.X = mat.NewDense(4, 4, nil)  // P
	d.X1 = mat.NewDense(4, 4, nil) // first der
	d.X2 = mat.NewDense(4, 4, nil) // second der
}

//SetRateMatrix needs to be done before doing SetupQGTR
// just send along the 5 rates and this will make them the whole matrix
func (d *DNAModel) SetRateMatrix(params []float64) {
	d.R = mat.NewDense(4, 4, nil)
	d.R.Set(0, 0, 0)
	d.R.Set(1, 1, 0)
	d.R.Set(2, 2, 0)
	d.R.Set(3, 3, 0)
	d.R.Set(0, 1, params[0])
	d.R.Set(1, 0, params[0])
	d.R.Set(0, 2, params[1])
	d.R.Set(2, 0, params[1])
	d.R.Set(0, 3, params[2])
	d.R.Set(3, 0, params[2])
	d.R.Set(1, 2, params[3])
	d.R.Set(2, 1, params[3])
	d.R.Set(1, 3, params[4])
	d.R.Set(3, 1, params[4])
	d.R.Set(2, 3, 1.)
	d.R.Set(3, 2, 1.)
}

// SetBaseFreqs needs to be done before doing SetupQGTR
func (d *DNAModel) SetBaseFreqs(basefreq []float64) {
	d.BF = basefreq
}

// ExpValue used for the matrix exponential
func (d *DNAModel) ExpValue(iv []float64, blen float64) {
	for i, j := range iv {
		d.X.Set(i, i, math.Exp(j*blen))
	}
	return
}

// ExpValueFirstD get the first derivaite for NR
func (d *DNAModel) ExpValueFirstD(blen float64) (x *mat.Dense) {
	x = mat.NewDense(4, 4, nil)
	for i := 0; i < 4; i++ {
		d.X1.Set(i, i, d.EigenVals[i]*math.Exp(d.EigenVals[i]*blen))
	}
	x.Mul(d.EigenVecs, d.X1)
	x.Mul(x, d.EigenVecsI)
	return
}

// ExpValueSecondD get the second derivaite for NR
func (d *DNAModel) ExpValueSecondD(blen float64) (x *mat.Dense) {
	x = mat.NewDense(4, 4, nil)
	for i := 0; i < 4; i++ {
		d.X1.Set(i, i, (d.EigenVals[i]*d.EigenVals[i])*math.Exp(d.EigenVals[i]*blen))
	}
	x.Mul(d.EigenVecs, d.X1)
	x.Mul(x, d.EigenVecsI)
	return
}

// SetP use the standard spectral decom
func (d *DNAModel) SetP(blen float64) {
	P := mat.NewDense(4, 4, nil)
	d.ExpValue(d.EigenVals, blen)
	P.Mul(d.EigenVecs, d.X)
	P.Mul(P, d.EigenVecsI)
	// easier
	//d.Lock()
	d.Ps[blen] = P
	//d.Unlock()
}

// SetPSimple use the gonum matrixexp (seems faster)
func (d *DNAModel) SetPSimple(blen float64) {
	P := mat.NewDense(4, 4, nil)
	P.Scale(blen, d.Q)
	P.Exp(P)
	//d.Lock()
	d.Ps[blen] = P
	//d.Unlock()
}

// EmptyPDict save memory perhaps?
func (d *DNAModel) EmptyPDict() {
	d.Ps = nil
	d.Ps = make(map[float64]*mat.Dense)
}

// EmptyPLDict the logged one
func (d *DNAModel) EmptyPLDict() {
	d.PsL = nil
	d.PsL = make(map[float64]*mat.Dense)
}

// GetPMap get the Ps from the dictionary
func (d *DNAModel) GetPMap(blen float64) *mat.Dense {
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

// GetPMapLogged get the Ps from the dictionary
func (d *DNAModel) GetPMapLogged(blen float64) *mat.Dense {
	if _, ok := d.PsL[blen]; ok {
		return d.PsL[blen]
	}
	P := mat.NewDense(4, 4, nil)
	P.Scale(blen, d.Q)
	P.Exp(P)
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			P.Set(i, j, math.Log(P.At(i, j)))
		}
	}
	d.PsL[blen] = P
	return d.PsL[blen]
}

// GetPCalc calculate P matrix
func (d *DNAModel) GetPCalc(blen float64) *mat.Dense {
	var P mat.Dense
	P.Scale(blen, d.Q)
	P.Exp(&P)
	//d.ExpValue(d.EigenVals, blen)
	//P.Mul(&d.EigenVecs, d.X)
	//P.Mul(&P, d.EigenVecsT)
	return &P
}

/* TODO: get this in there
	R,r ==> {AG}
   Y,y ==> {CT}
   M,m ==> {AC}
   K,k ==> {GT}
   S,s ==> {CG}
   W,w ==> {AT}
   H,h ==> {ACT}
   B,b ==> {CGT}
   V,v ==> {ACG}
   D,d ==> {AGT}
   N,n ==> {ACGT}
*/

//SetMap for getting the position in the array
func (d *DNAModel) SetMap() {
	d.CharMap = make(map[string][]int)
	d.CharMap["A"] = []int{0}
	d.CharMap["C"] = []int{1}
	d.CharMap["G"] = []int{2}
	d.CharMap["T"] = []int{3}
	d.CharMap["-"] = []int{0, 1, 2, 3}
	d.CharMap["N"] = []int{0, 1, 2, 3}
	d.CharMap["R"] = []int{0, 2}
	d.CharMap["Y"] = []int{1, 3}
	d.CharMap["M"] = []int{0, 1}
	d.CharMap["K"] = []int{2, 3}
	d.CharMap["S"] = []int{1, 2}
	d.CharMap["W"] = []int{0, 3}
	d.CharMap["H"] = []int{0, 1, 3}
	d.CharMap["B"] = []int{1, 2, 3}
	d.CharMap["V"] = []int{0, 1, 2}
	d.CharMap["D"] = []int{0, 2, 3}
}

func (d *DNAModel) GetCharMap() map[string][]int {
	return d.CharMap
}

func (d *DNAModel) GetNumStates() int {
	return d.NumStates
}

func (d *DNAModel) GetBF() []float64 {
	return d.BF
}

func (d *DNAModel) GetStochMapMatrices(dur float64, from int, to int) (summed *mat.Dense, summedR *mat.Dense) {
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
