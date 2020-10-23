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
	seqs, patternsint, _, bf := gophy.ReadPatternsSeqsFromFile(afn, true)
	patternval, _ := gophy.PreparePatternVecs(tr, patternsint, seqs)
	x := gophy.NewDNAModel()
	x.M.SetBaseFreqs(bf)
	modelparams := make([]float64, 5)
	for i := range modelparams {
		modelparams[i] = 1.0
	}
	x.M.SetRateMatrix(modelparams)
	x.M.SetupQGTR()
	lnl := gophy.PCalcLogLikePatterns(tr, &x.M, patternval, 2)
	if math.Round(lnl*1000)/1000 != -4570.796 {
		fmt.Println(lnl)
		t.Fail()
	}
}

func TestPCalcLogLikePatterns(t *testing.T) {
	tfn := "test_files/10tips.nuc.fa.treefile"
	tr := gophy.ReadTreeFromFile(tfn)
	afn := "test_files/10tips.nuc.fa"
	seqs, patternsint, _, bf := gophy.ReadPatternsSeqsFromFile(afn, true)
	patternval, _ := gophy.PreparePatternVecs(tr, patternsint, seqs)
	x := gophy.NewDNAModel()
	x.M.SetBaseFreqs(bf)
	modelparams := make([]float64, 5)
	for i := range modelparams {
		modelparams[i] = 1.0
	}
	x.M.SetRateMatrix(modelparams)
	x.M.SetupQGTR()
	lnl := gophy.PCalcLikePatterns(tr, &x.M, patternval, 2)
	fmt.Println(lnl)
	if math.Round(lnl*1000)/1000 != -4570.796 {
		fmt.Println(lnl)
		t.Fail()
	}
}

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
	x := gophy.NewMultStateModel(numstates)
	x.M.SetScaledRateMatrix(modelparams, true)
	x.M.SetBaseFreqs(bf)
	x.M.SetupQGTR()
	lnl := gophy.PCalcLogLikePatterns(tr, &x.M, patternval, 2)
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
	x := gophy.NewMultStateModel(numstates)
	x.M.SetScaledRateMatrix(modelparams, true)
	x.M.SetBaseFreqs(bf)
	x.M.EBF = bf
	x.M.SetupQGTR()
	lnl := gophy.PCalcLikePatterns(tr, &x.M, patternval, 2)
	if math.Round(lnl*1000)/1000 != -4570.796 {
		t.Fail()
	}
}
