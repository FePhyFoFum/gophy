package gophy

import (
	"gonum.org/v1/gonum/mat"
)

type DNAModel struct {
	BF      mat.Vector // base frequencies
	R       *mat.Dense
	Q       *mat.Dense
	CharMap map[string]int
	Ps      map[float64]mat.Dense
}

func (d *DNAModel) SetupQ() {
	//just JC for now
	d.Ps = make(map[float64]mat.Dense)
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

func (d *DNAModel) SetP(blen float64) {
	X := mat.NewDense(4, 4, nil)
	X.Scale(blen, d.Q)
	var P mat.Dense
	P.Exp(X)
	d.Ps[blen] = P
}

func (d *DNAModel) GetPMap(blen float64) mat.Dense {
	if _, ok := d.Ps[blen]; ok {
		return d.Ps[blen]
	} else {
		d.SetP(blen)
	}
	return d.Ps[blen]
}

func (d *DNAModel) GetPCalc(blen float64) mat.Dense {
	X := mat.NewDense(4, 4, nil)
	X.Scale(blen, d.Q)
	var P mat.Dense
	P.Exp(X)
	return P
}

func (d *DNAModel) SetNucMap() {
	d.CharMap = map[string]int{"A": 0, "C": 1, "G": 2, "T": 3}
}
