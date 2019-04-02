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

// NewDNAModel get new DNAModel pointer
func NewDNAModel() *DNAModel {
	return &DNAModel{}
}

// SetupQJC setup Q matrix
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

// this is just for NR optimization for branch lengths
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
	ES.Factorize(d.QS, true, true)
	TC := ES.VectorsTo(nil)
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

//SetNucMap for getting the position in the array
func (d *DNAModel) SetNucMap() {
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
