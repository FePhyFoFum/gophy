package gophy

import (
	"strconv"
)

// MultStateModel multistate model struct
type MultStateModel struct {
	DiscreteModel
}

// NewMultStateModel get new MULTModel pointer
func NewMultStateModel(numstates int) *MultStateModel {
	CharMap := make(map[string][]int)
	CharMap["-"] = make([]int, numstates)
	CharMap["N"] = make([]int, numstates)
	for i := 0; i < numstates; i++ {
		CharMap[strconv.Itoa(i)] = []int{i}
		CharMap["-"][i] = i
		CharMap["N"][i] = i
	}
	dnam := MultStateModel{}
	dnam.DiscreteModel.Alph = MultiState
	dnam.DiscreteModel.NumStates = numstates
	dnam.DiscreteModel.CharMap = CharMap
	return &dnam
}
