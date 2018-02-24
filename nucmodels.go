package gophy

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

// DNAModel standard DNA struct
type DNAModel struct {
	BF      mat.Vector // base frequencies
	R       *mat.Dense
	Q       *mat.Dense
	QS      *mat.SymDense
	CharMap map[string]int
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

// SetupQ setup Q matrix
func (d *DNAModel) SetupQ() {
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

// SetNucMap for getting the position in the array
func (d *DNAModel) SetNucMap() {
	d.CharMap = map[string]int{"A": 0, "C": 1, "G": 2, "T": 3}
}
