package gophy

import "gonum.org/v1/gonum/mat"

//interfact for model (DNAModel, MULTModel)

// StateModel ...
type StateModel interface {
	GetPCalc(float64) *mat.Dense
	GetPMap(float64) *mat.Dense
	SetMap()
	GetNumStates() int
	GetBF() []float64
	EmptyPDict()
	GetStochMapMatrices(float64, int, int) (*mat.Dense, *mat.Dense)
	SetupQJC1Rate(float64)
	SetupQMk([]float64, bool)            // bool = false is AsyMK
	SetScaledRateMatrix([]float64, bool) // before setupQGTR
	SetupQGTR()                          //
	DecomposeQ()                         //
	ExpValueFirstD(float64) *mat.Dense   //BL
	ExpValueSecondD(float64) *mat.Dense  //BL
}
