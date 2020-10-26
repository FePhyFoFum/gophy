package gophy

import "gonum.org/v1/gonum/mat"

type DNAModel struct {
	M DiscreteModel
}

func NewDNAModel() *DNAModel {
	CharMap := make(map[string][]int)
	CharMap["A"] = []int{0}
	CharMap["C"] = []int{1}
	CharMap["G"] = []int{2}
	CharMap["T"] = []int{3}
	CharMap["-"] = []int{0, 1, 2, 3}
	CharMap["N"] = []int{0, 1, 2, 3}
	CharMap["R"] = []int{0, 2}
	CharMap["Y"] = []int{1, 3}
	CharMap["M"] = []int{0, 1}
	CharMap["K"] = []int{2, 3}
	CharMap["S"] = []int{1, 2}
	CharMap["W"] = []int{0, 3}
	CharMap["H"] = []int{0, 1, 3}
	CharMap["B"] = []int{1, 2, 3}
	CharMap["V"] = []int{0, 1, 2}
	CharMap["D"] = []int{0, 2, 3}
	dnam := DNAModel{}
	dnam.M.Alph = Nucleotide
	dnam.M.NumStates = 4
	dnam.M.CharMap = CharMap
	return &dnam
}

/* this should be set by the SetRateMatrix in discrete_model
func (d *DNAModelNew) SetRateMatrix(params []float64) {
	R := mat.NewDense(4, 4, nil)
	R.Set(0, 0, 0)
	R.Set(1, 1, 0)
	R.Set(2, 2, 0)
	R.Set(3, 3, 0)
	R.Set(0, 1, params[0])
	R.Set(1, 0, params[0])
	R.Set(0, 2, params[1])
	R.Set(2, 0, params[1])
	R.Set(0, 3, params[2])
	R.Set(3, 0, params[2])
	R.Set(1, 2, params[3])
	R.Set(2, 1, params[3])
	R.Set(1, 3, params[4])
	R.Set(3, 1, params[4])
	R.Set(2, 3, 1.)
	R.Set(3, 2, 1.)
	d.DiscreteModel.R = R
}
*/

// DeepCopyDNAModel ...
func (d *DNAModel) DeepCopyDNAModel() *DNAModel {
	outm := NewDNAModel()
	outm.M.BF = []float64{0.25, 0.25, 0.25, 0.25}
	copy(outm.M.BF, d.M.BF)
	outm.M.Q = mat.NewDense(4, 4, nil)
	outm.M.R = mat.NewDense(4, 4, nil)
	outm.M.R.Copy(d.M.R)
	outm.M.Q.Copy(d.M.Q)
	return outm
}
