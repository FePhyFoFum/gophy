package gophy_test

import (
	"fmt"
	"math"
	"testing"

	"github.com/FePhyFoFum/gophy"
)

func TestPCalcLikePatterns(t *testing.T) {
	tfn := "test_files/10tips.nuc.fa.treefile"
	tr := gophy.ReadTreeFromFile(tfn)
	afn := "test_files/10tips.nuc.fa"
	seqs, patternsint, _, bf := gophy.ReadPatternsSeqsFromFile(afn)
	patternval, _ := gophy.PreparePatternVecs(tr, patternsint, seqs)
	x := gophy.NewModel()
	x.Alph = "nuc"
	x.NumStates = 4
	x.SetMapDNA()
	x.SetBaseFreqs(bf)
	modelparams := make([]float64, 5)
	for i := range modelparams {
		modelparams[i] = 1.0
	}
	x.SetRateMatrix(modelparams)
	x.SetupQGTR()
	lnl := gophy.PCalcLogLikePatterns(tr, x, patternval, 2)
	if math.Round(lnl*1000)/1000 != -4570.796 {
		fmt.Println(lnl)
		t.Fail()
	}
}

func TestPCalcLogLikePatterns(t *testing.T) {
	tfn := "test_files/10tips.nuc.fa.treefile"
	tr := gophy.ReadTreeFromFile(tfn)
	afn := "test_files/10tips.nuc.fa"
	seqs, patternsint, _, bf := gophy.ReadPatternsSeqsFromFile(afn)
	patternval, _ := gophy.PreparePatternVecs(tr, patternsint, seqs)
	x := gophy.NewModel()
	x.Alph = "nuc"
	x.NumStates = 4
	x.SetMapDNA()
	x.SetBaseFreqs(bf)
	modelparams := make([]float64, 5)
	for i := range modelparams {
		modelparams[i] = 1.0
	}
	x.SetRateMatrix(modelparams)
	x.SetupQGTR()
	lnl := gophy.PCalcLikePatterns(tr, x, patternval, 2)
	fmt.Println(lnl)
	if math.Round(lnl*1000)/1000 != -4570.796 {
		fmt.Println(lnl)
		t.Fail()
	}
}
