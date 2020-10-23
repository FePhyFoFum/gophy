package gophy

import (
	"strconv"
)

// MultStateModel multistate model struct
type MultStateModelNew struct {
	DiscreteModelNew
}

// NewMultStateModel get new MULTModel pointer
func NewMultStateModelNew(numstates int) *MultStateModelNew {
	CharMap := make(map[string][]int)
	CharMap["-"] = make([]int, numstates)
	CharMap["N"] = make([]int, numstates)
	for i := 0; i < numstates; i++ {
		CharMap[strconv.Itoa(i)] = []int{i}
		CharMap["-"][i] = i
		CharMap["N"][i] = i
	}
	dnam := MultStateModelNew{}
	dnam.DiscreteModelNew.Alph = MultiState
	dnam.DiscreteModelNew.NumStates = numstates
	dnam.DiscreteModelNew.CharMap = CharMap
	return &dnam
}
