package gophy_test

import (
	"fmt"
	"math"
	"testing"

	"github.com/FePhyFoFum/gophy"
)

func TestPCalcLikePatternsMS(t *testing.T) {
	tfn := "test_files/10tips.nuc.fa.treefile"
	tr := gophy.ReadTreeFromFile(tfn)
	afn := "test_files/10tips.fa"
	seqs, patternsint, _, bf, numstates := gophy.ReadPatternsMSeqsFromFile(afn)
	patternval, _ := gophy.PreparePatternVecsMS(tr, patternsint, seqs, gophy.GetMap(numstates), numstates)
	modelparams := make([]float64, 5)
	for i := range modelparams {
		modelparams[i] = 1.0
	}
	x := gophy.NewModel()
	x.Alph = "mult"
	x.NumStates = numstates
	x.SetScaledRateMatrix(modelparams, true)
	x.SetMapMult()
	x.SetBaseFreqs(bf)
	x.SetupQGTR()
	lnl := gophy.PCalcLogLikePatternsMS(tr, x, patternval, 2)
	if math.Round(lnl*1000)/1000 != -4570.796 {
		fmt.Println(lnl)
		t.Fail()
	}
}

func TestPCalcLogLikePatternsMS(t *testing.T) {
	tfn := "test_files/10tips.nuc.fa.treefile"
	tr := gophy.ReadTreeFromFile(tfn)
	afn := "test_files/10tips.fa"
	seqs, patternsint, _, bf, numstates := gophy.ReadPatternsMSeqsFromFile(afn)
	patternval, _ := gophy.PreparePatternVecsMS(tr, patternsint, seqs, gophy.GetMap(numstates), numstates)
	modelparams := make([]float64, 5)
	for i := range modelparams {
		modelparams[i] = 1.0
	}
	x := gophy.NewModel()
	x.Alph = "mult"
	x.NumStates = numstates
	x.SetScaledRateMatrix(modelparams, true)
	x.SetMapMult()
	x.SetBaseFreqs(bf)
	x.EBF = bf
	x.SetupQGTR()
	lnl := gophy.PCalcLikePatternsMS(tr, x, patternval, 2)
	if math.Round(lnl*1000)/1000 != -4570.796 {
		t.Fail()
	}
}
