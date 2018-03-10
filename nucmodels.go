package gophy

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

// DNAModel standard DNA struct
type DNAModel struct {
	BF      []float64 // base frequencies
	R       *mat.Dense
	Q       *mat.Dense
	QS      *mat.SymDense
	CharMap map[string][]int
	//sync.RWMutex
	Ps         map[float64]*mat.Dense
	EigenVals  []float64 //to be exponentiated
	EigenVecs  mat.Dense
	EigenVecsT mat.Matrix
	X          *mat.SymDense
	P          mat.Dense
}

// NewDNAModel get new DNAModel pointer
func NewDNAModel() *DNAModel {
	return &DNAModel{}
}

// SetupQJC setup Q matrix
func (d *DNAModel) SetupQJC() {
	d.BF = []float64{0.25, 0.25, 0.25, 0.25}
	//just JC for now
	d.Ps = make(map[float64]*mat.Dense)
	d.Q = mat.NewDense(4, 4, nil)
	//d.QS = mat.NewSymDense(4, nil)
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if i != j {
				d.Q.Set(i, j, 0.333333333333)
				//d.QS.SetSym(i, j, 0.333333333333)
			} else {
				d.Q.Set(i, j, -1.)
				//d.QS.SetSym(i, j, -1.)
			}
		}
	}
	//decompose, each time you change the model
	//var ES mat.EigenSym
	//ES.Factorize(d.QS, true)
	//d.EigenVecs.EigenvectorsSym(&ES)
	//d.EigenVecsT = d.EigenVecs.T()
	//d.EigenVals = ES.Values(nil)
	//d.X = mat.NewSymDense(4, nil)
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

// ExpValue used for the matrix exponential
func (d *DNAModel) ExpValue(iv []float64, blen float64) {
	for i, j := range iv {
		d.X.SetSym(i, i, math.Exp(j*blen))
	}
	return
}

// ExpValueFirstD get the first derivaite for NR
func ExpValueFirstD(d []float64, blen float64) (x *mat.SymDense) {
	x = mat.NewSymDense(4, nil)
	for i, j := range d {
		x.SetSym(i, i, j*math.Exp(j*blen))
	}
	return
}

// ExpValueSecondD get the second derivative for NR
func ExpValueSecondD(d []float64, blen float64) (x *mat.SymDense) {
	x = mat.NewSymDense(4, nil)
	for i, j := range d {
		x.SetSym(i, i, (j*j)*math.Exp(j*blen))
	}
	return
}

// SetP use the standard spectral decom
func (d *DNAModel) SetP(blen float64) {
	P := mat.NewDense(4, 4, nil)
	d.ExpValue(d.EigenVals, blen)
	P.Mul(&d.EigenVecs, d.X)
	P.Mul(P, d.EigenVecsT)
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
// SetNucMap for getting the position in the array
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
