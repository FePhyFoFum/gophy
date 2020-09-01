package gophy

import "gonum.org/v1/gonum/mat"

//interfact for model (DNAModel, MultStateModel, AAModel)

// StateModel ...
type StateModel interface {
	GetPCalc(float64) *mat.Dense
	GetPMap(float64) *mat.Dense
	SetMapDNA()
	SetMapProt()
	GetNumStates() int
	GetBF() []float64
	EmptyPDict()
	GetStochMapMatrices(float64, int, int) (*mat.Dense, *mat.Dense)
	SetupQJC1Rate(float64)
	SetupQMk([]float64, bool) // bool = false is AsyMK
	SetRateMatrix([]float64)
	SetScaledRateMatrix([]float64, bool) // before setupQGTR
	SetRateMatrixDNA([]float64)
	SetRateMatrixJTT()
	SetRateMatrixWAG()
	SetRateMatrixLG()
	SetupQGTR()                         //
	DecomposeQ()                        //
	ExpValueFirstD(float64) *mat.Dense  //BL
	ExpValueSecondD(float64) *mat.Dense //BL
	GetCharMap() map[string][]int
	SetBaseFreqs([]float64)
}
