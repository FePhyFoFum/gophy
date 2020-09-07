package gophy

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/mat"
)

// DataType type for alphabet
type DataType string

// datatype constants
const (
	Nucleotide DataType = "nuc"
	AminoAcid           = "aa"
	MultiState          = "mult"
)

// Model overall model struct
type Model struct {
	Alph      DataType  // nuc, prot or multstate model
	BF        []float64 // base frequencies, order is A,C,G,T or A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V
	MBF       []float64 // model base frequencies
	EBF       []float64 // empirical base freqs
	R         *mat.Dense
	Q         *mat.Dense // common use
	Ex        string     // PAML-formatted exchangeabilities for Prot models
	CharMap   map[string][]int
	NumStates int
	//sync.RWMutex
	Ps  map[float64]*mat.Dense
	PsL map[float64]*mat.Dense
	X   *mat.Dense
	P   mat.Dense
	//for decomposing
	QS         *mat.Dense
	EigenVals  []float64  // to be exponentiated
	EigenVecs  *mat.Dense //
	EigenVecsI *mat.Dense
	X1         *mat.Dense
	X2         *mat.Dense
}

// DeepCopyDNAModel ...
func (d *Model) DeepCopyDNAModel() *Model {
	outm := NewModel()
	outm.Alph = "nuc"
	outm.BF = []float64{0.25, 0.25, 0.25, 0.25}
	copy(outm.BF, d.BF)
	outm.Q = mat.NewDense(4, 4, nil)
	outm.R = mat.NewDense(4, 4, nil)
	outm.R.Copy(d.R)
	outm.Q.Copy(d.Q)
	outm.SetMapDNA()
	return outm
}

// NewModel get new model pointer
func NewModel() *Model {
	return &Model{}
}

// SetEqualBF set all state frequencies equal
func (d *Model) SetEqualBF() {
	for i := range d.BF {
		d.BF[i] = 1. / float64(d.NumStates)
	}
}

//SetEmpiricalBF set all to empirical
func (d *Model) SetEmpiricalBF() {
	d.BF = d.EBF
}

// SetModelBF set all to frequencies from empirical model
func (d *Model) SetModelBF() {
	d.BF = d.MBF
}

// SetBaseFreqs needs to be done before doing SetupQGTR
func (d *Model) SetBaseFreqs(basefreq []float64) {
	d.BF = basefreq
}

// SetupQJC setup Q matrix
//    This is scaled so that change is reflected in the branch lengths
//    You don't need to use the SetScaledRateMatrix
func (d *Model) SetupQJC() {
	bf := []float64{}
	for i := 0; i < d.NumStates; i++ {
		bf[i] = 1. / float64(d.NumStates)
	}
	d.BF = bf
	d.Ps = make(map[float64]*mat.Dense)
	d.Q = mat.NewDense(d.NumStates, d.NumStates, nil)

	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			if i != j {
				d.Q.Set(i, j, 1./float64(d.NumStates-1))
			} else {
				d.Q.Set(i, j, -1.)
			}
		}
	}
}

// SetupQJC1Rate setup Q matrix with one rate, probably for anc multi state
//     These are unscaled so the branch lengths are going to be time or something else
//     and not relative to these rates
//     Will take BF from something else
func (d *Model) SetupQJC1Rate(rt float64) {
	d.Ps = make(map[float64]*mat.Dense)
	d.Q = mat.NewDense(d.NumStates, d.NumStates, nil)

	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			if i != j {
				d.Q.Set(i, j, rt)
			} else {
				d.Q.Set(i, j, -(float64(d.NumStates-1) * rt))
			}
		}
	}
}

// SetupQGTR sets up scaled GTR for (mainly) nucleotide models
// Could also be used to estimate PROTGTR but that's a bad idea
func (d *Model) SetupQGTR() {
	// NWH code, Foster-style scaling
	bigpi := mat.NewDiagDense(d.NumStates, d.BF)
	dQ := mat.NewDense(d.NumStates, d.NumStates, nil)
	dQ.Mul(d.R, bigpi)
	var offdSum float64
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			if i != j {
				offdSum += dQ.At(i, j)
			}
		}
	}
	var dSum float64
	d.Q = mat.NewDense(d.NumStates, d.NumStates, nil)
	for i := 0; i < d.NumStates; i++ {
		d.Q.Set(i, i, 0-sumRow(dQ, i))
		dSum += d.Q.At(i, i) * d.BF[i]
	}
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			if i == j {
				d.Q.Set(i, i, d.Q.At(i, i)/-dSum)
			} else {
				d.Q.Set(i, j, dQ.At(i, j)/-dSum)
			}
		}
	}
}

// SetupQMk setup Q matrix
//    This is unscaled (so the branch lengths are going to be proportion to some other change
//    and not to these branch lengths)
//    Will take the BF from something else
func (d *Model) SetupQMk(rt []float64, sym bool) {
	d.Ps = make(map[float64]*mat.Dense)
	d.Q = mat.NewDense(d.NumStates, d.NumStates, nil)
	cc := 0
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			if i != j {
				if sym && j > i {
					d.Q.Set(i, j, rt[cc])
					d.Q.Set(j, i, rt[cc])
					cc++
				} else if sym == false {
					d.Q.Set(i, j, rt[cc])
					cc++
				}
			} else {
				d.Q.Set(i, j, 0.0)
			}
		}
	}
	for i := 0; i < d.NumStates; i++ {
		sumrow := 0.
		for j := 0; j < d.NumStates; j++ {
			sumrow += d.Q.At(i, j)
		}
		d.Q.Set(i, i, -sumrow)
	}
}

// DecomposeQ this is just for NR optimization for branch lengths
func (d *Model) DecomposeQ() {
	d.QS = mat.NewDense(d.NumStates, d.NumStates, nil)
	d.EigenVecsI = mat.NewDense(d.NumStates, d.NumStates, nil)
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			d.QS.Set(i, j, d.Q.At(i, j))
		}
	}
	//decompose, each time you change the model
	var ES mat.Eigen
	d.EigenVecs = mat.NewDense(d.NumStates, d.NumStates, nil)
	ES.Factorize(d.QS, mat.EigenBoth) //true, true)
	TC := mat.NewCDense(d.NumStates, d.NumStates, nil)
	//	TC := ES.VectorsTo(nil)
	ES.VectorsTo(TC)
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			d.EigenVecs.Set(i, j, real(TC.At(i, j)))
		}
	}
	d.EigenVecsI.Inverse(d.EigenVecs)
	d.EigenVals = make([]float64, d.NumStates)
	for i := 0; i < d.NumStates; i++ {
		d.EigenVals[i] = 1.
	}
	TV := ES.Values(nil)
	for i := 0; i < d.NumStates; i++ {
		d.EigenVals[i] = real(TV[i])
	}

	d.X = mat.NewDense(d.NumStates, d.NumStates, nil)  // P
	d.X1 = mat.NewDense(d.NumStates, d.NumStates, nil) // first der
	d.X2 = mat.NewDense(d.NumStates, d.NumStates, nil) // second der
}

//SetRateMatrixDNA needs to be done before doing SetupQGTR
// just send along the 5 rates and this will make them the whole matrix
// NWH TO DO version of this for AA rates
func (d *Model) SetRateMatrixDNA(params []float64) {
	d.R = mat.NewDense(4, 4, nil)
	d.R.Set(0, 0, 0)
	d.R.Set(1, 1, 0)
	d.R.Set(2, 2, 0)
	d.R.Set(3, 3, 0)
	d.R.Set(0, 1, params[0])
	d.R.Set(1, 0, params[0])
	d.R.Set(0, 2, params[1])
	d.R.Set(2, 0, params[1])
	d.R.Set(0, 3, params[2])
	d.R.Set(3, 0, params[2])
	d.R.Set(1, 2, params[3])
	d.R.Set(2, 1, params[3])
	d.R.Set(1, 3, params[4])
	d.R.Set(3, 1, params[4])
	d.R.Set(2, 3, 1.)
	d.R.Set(3, 2, 1.)
}

// SetRateMatrixJTT set up JTT exchangeabilities
func (d *Model) SetRateMatrixJTT() {
	d.Ex = `58
	54  45
	81  16 528
	56 113  34  10
	57 310  86  49   9
	105  29  58 767   5 323
	179 137  81 130  59  26 119
	27 328 391 112  69 597  26  23
	36  22  47  11  17   9  12   6  16
	30  38  12   7  23  72   9   6  56 229
	35 646 263  26   7 292 181  27  45  21  14
	54  44  30  15  31  43  18  14  33 479 388  65
	15   5  10   4  78   4   5   5  40  89 248   4  43
	194  74  15  15  14 164  18  24 115  10 102  21  16  17
	378 101 503  59 223  53  30 201  73  40  59  47  29  92 285
	475  64 232  38  42  51  32  33  46 245  25 103 226  12 118 477
	9 126   8   4 115  18  10  55   8   9  52  10  24  53   6  35  12
	11  20  70  46 209  24   7   8 573  32  24   8  18 536  10  63  21  71
	298  17  16  31  62  20  45  47  11 961 180  14 323  62  23  38 112  25  16
	
	0.076748 0.051691 0.042645 0.051544 0.019803 0.040752 0.061830 0.073152 0.022944 0.053761 0.091904 0.058676 0.023826 0.040126 0.050901 0.068765 0.058565 0.014261 0.032102 0.066005`
	// from paml jones.dat

	scanner := bufio.NewScanner(strings.NewReader(d.Ex))
	res := [][]float64{}
	for scanner.Scan() {
		tmp := []float64{}
		line := strings.Fields(scanner.Text())
		for _, j := range line {
			f, err := strconv.ParseFloat(j, 64)
			if err != nil {
				log.Fatal(err)
			}
			tmp = append(tmp, f)
		}
		res = append(res, tmp)
	}
	d.R = mat.NewDense(20, 20, nil) // exchangeabilities - let diagonals be 0
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			if i == j {
				d.R.Set(i, j, 0.0)
			} else if j > i {
				d.R.Set(i, j, res[j-1][i])
			} else {
				d.R.Set(i, j, res[i-1][j])
			}
		}
	}
	d.MBF = res[20]
}

// SetRateMatrixWAG set up WAG exchangeabilities
func (d *Model) SetRateMatrixWAG() {
	d.Ex = `0.551571
	0.509848  0.635346
	0.738998  0.147304  5.429420
	1.027040  0.528191  0.265256  0.0302949
	0.908598  3.035500  1.543640  0.616783  0.0988179
	1.582850  0.439157  0.947198  6.174160  0.021352  5.469470
	1.416720  0.584665  1.125560  0.865584  0.306674  0.330052  0.567717
	0.316954  2.137150  3.956290  0.930676  0.248972  4.294110  0.570025  0.249410
	0.193335  0.186979  0.554236  0.039437  0.170135  0.113917  0.127395  0.0304501 0.138190
	0.397915  0.497671  0.131528  0.0848047 0.384287  0.869489  0.154263  0.0613037 0.499462  3.170970
	0.906265  5.351420  3.012010  0.479855  0.0740339 3.894900  2.584430  0.373558  0.890432  0.323832  0.257555
	0.893496  0.683162  0.198221  0.103754  0.390482  1.545260  0.315124  0.174100  0.404141  4.257460  4.854020  0.934276
	0.210494  0.102711  0.0961621 0.0467304 0.398020  0.0999208 0.0811339 0.049931  0.679371  1.059470  2.115170  0.088836  1.190630
	1.438550  0.679489  0.195081  0.423984  0.109404  0.933372  0.682355  0.243570  0.696198  0.0999288 0.415844  0.556896  0.171329  0.161444
	3.370790  1.224190  3.974230  1.071760  1.407660  1.028870  0.704939  1.341820  0.740169  0.319440  0.344739  0.967130  0.493905  0.545931  1.613280
	2.121110  0.554413  2.030060  0.374866  0.512984  0.857928  0.822765  0.225833  0.473307  1.458160  0.326622  1.386980  1.516120  0.171903  0.795384  4.378020
	0.113133  1.163920  0.0719167 0.129767  0.717070  0.215737  0.156557  0.336983  0.262569  0.212483  0.665309  0.137505  0.515706  1.529640  0.139405  0.523742  0.110864
	0.240735  0.381533  1.086000  0.325711  0.543833  0.227710  0.196303  0.103604  3.873440  0.420170  0.398618  0.133264  0.428437  6.454280  0.216046  0.786993  0.291148  2.485390
	2.006010  0.251849  0.196246  0.152335  1.002140  0.301281  0.588731  0.187247  0.118358  7.821300  1.800340  0.305434  2.058450  0.649892  0.314887  0.232739  1.388230  0.365369  0.314730
	
	0.0866279 0.043972  0.0390894 0.0570451 0.0193078 0.0367281 0.0580589 0.0832518 0.0244313 0.048466  0.086209  0.0620286 0.0195027 0.0384319 0.0457631 0.0695179 0.0610127 0.0143859 0.0352742 0.0708956`
	// from PAML wag.dat

	scanner := bufio.NewScanner(strings.NewReader(d.Ex))
	res := [][]float64{}
	for scanner.Scan() {
		tmp := []float64{}
		line := strings.Fields(scanner.Text())
		for _, j := range line {
			f, err := strconv.ParseFloat(j, 64)
			if err != nil {
				log.Fatal(err)
			}
			tmp = append(tmp, f)
		}
		res = append(res, tmp)
	}
	d.R = mat.NewDense(20, 20, nil) // exchangeabilities - let diagonals be 0
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			if i == j {
				d.R.Set(i, j, 0.0)
			} else if j > i {
				d.R.Set(i, j, res[j-1][i])
			} else {
				d.R.Set(i, j, res[i-1][j])
			}
		}
	}
	d.MBF = res[20]
}

// SetRateMatrixLG set up LG exchangeabilities
func (d *Model) SetRateMatrixLG() {
	d.Ex = `0.425093
	0.276818 0.751878
	0.395144 0.123954 5.076149
	2.489084 0.534551 0.528768 0.062556
	0.969894 2.807908 1.695752 0.523386 0.084808
	1.038545 0.363970 0.541712 5.243870 0.003499 4.128591
	2.066040 0.390192 1.437645 0.844926 0.569265 0.267959 0.348847
	0.358858 2.426601 4.509238 0.927114 0.640543 4.813505 0.423881 0.311484
	0.149830 0.126991 0.191503 0.010690 0.320627 0.072854 0.044265 0.008705 0.108882
	0.395337 0.301848 0.068427 0.015076 0.594007 0.582457 0.069673 0.044261 0.366317 4.145067
	0.536518 6.326067 2.145078 0.282959 0.013266 3.234294 1.807177 0.296636 0.697264 0.159069 0.137500
	1.124035 0.484133 0.371004 0.025548 0.893680 1.672569 0.173735 0.139538 0.442472 4.273607 6.312358 0.656604
	0.253701 0.052722 0.089525 0.017416 1.105251 0.035855 0.018811 0.089586 0.682139 1.112727 2.592692 0.023918 1.798853
	1.177651 0.332533 0.161787 0.394456 0.075382 0.624294 0.419409 0.196961 0.508851 0.078281 0.249060 0.390322 0.099849 0.094464
	4.727182 0.858151 4.008358 1.240275 2.784478 1.223828 0.611973 1.739990 0.990012 0.064105 0.182287 0.748683 0.346960 0.361819 1.338132
	2.139501 0.578987 2.000679 0.425860 1.143480 1.080136 0.604545 0.129836 0.584262 1.033739 0.302936 1.136863 2.020366 0.165001 0.571468 6.472279
	0.180717 0.593607 0.045376 0.029890 0.670128 0.236199 0.077852 0.268491 0.597054 0.111660 0.619632 0.049906 0.696175 2.457121 0.095131 0.248862 0.140825
	0.218959 0.314440 0.612025 0.135107 1.165532 0.257336 0.120037 0.054679 5.306834 0.232523 0.299648 0.131932 0.481306 7.803902 0.089613 0.400547 0.245841 3.151815
	2.547870 0.170887 0.083688 0.037967 1.959291 0.210332 0.245034 0.076701 0.119013 10.649107 1.702745 0.185202 1.898718 0.654683 0.296501 0.098369 2.188158 0.189510 0.249313
	
	0.079066 0.055941 0.041977 0.053052 0.012937 0.040767 0.071586 0.057337 0.022355 0.062157 0.099081 0.064600 0.022951 0.042302 0.044040 0.061197 0.053287 0.012066 0.034155 0.069147`
	// lg.dat from paml

	scanner := bufio.NewScanner(strings.NewReader(d.Ex))
	res := [][]float64{}
	for scanner.Scan() {
		tmp := []float64{}
		line := strings.Fields(scanner.Text())
		for _, j := range line {
			f, err := strconv.ParseFloat(j, 64)
			if err != nil {
				log.Fatal(err)
			}
			tmp = append(tmp, f)
		}
		res = append(res, tmp)
	}
	d.R = mat.NewDense(20, 20, nil) // exchangeabilities - let diagonals be 0
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			if i == j {
				d.R.Set(i, j, 0.0)
			} else if j > i {
				d.R.Set(i, j, res[j-1][i])
			} else {
				d.R.Set(i, j, res[i-1][j])
			}
		}
	}
	d.MBF = res[20]
}

//SetRateMatrix length is (((numstates*numstates)-numstates)/2) - 1
// or (numstates * numstates) - numstates
// this is for scaled branch lengths and matrices
// CHANGE THIS TO UNSCALED
func (d *Model) SetRateMatrix(params []float64) {
	d.R = mat.NewDense(d.NumStates, d.NumStates, nil)
	if len(params) == (((d.NumStates*d.NumStates)-d.NumStates)/2)-1 {
		//symm
		pcount := 0
		for i := 0; i < d.NumStates; i++ {
			for j := 0; j < d.NumStates; j++ {
				if j > i {
					continue
				}
				if i == j {
					d.R.Set(i, j, 0.0)
				} else if i == d.NumStates-1 && j == d.NumStates-2 {
					d.R.Set(i, j, 1.0)
					d.R.Set(j, i, 1.0)
				} else {
					d.R.Set(i, j, params[pcount])
					d.R.Set(j, i, params[pcount])
					pcount++
				}
			}
		}
	} else if len(params) == ((d.NumStates*d.NumStates)-d.NumStates)-1 {
		//nonsymm
		pcount := 0
		for i := 0; i < d.NumStates; i++ {
			for j := 0; j < d.NumStates; j++ {
				if i == j {
					d.R.Set(i, j, 0.0)
				} else if i == d.NumStates-1 && j == d.NumStates-2 {
					d.R.Set(i, j, 1.0)
				} else {
					d.R.Set(i, j, params[pcount])
					pcount++
				}
			}
		}
	} else {
		fmt.Println("WRONG MATRIX SIZE")
		os.Exit(1)
	}
}

//SetScaledRateMatrix needs to be done before doing SetupQGTR
// just send along the rates and this will make them the whole matrix
// the scaled is that this is assuming that the last rate is 1
//THIS IS THE SAME? AS SETRATEMATRIX?
func (d *Model) SetScaledRateMatrix(params []float64, sym bool) {
	d.R = mat.NewDense(d.NumStates, d.NumStates, nil)
	cc := 0
	lasti := 0
	lastj := 0
	if sym {
		lasti = d.NumStates - 2
		lastj = d.NumStates - 1
	} else {
		lasti = d.NumStates - 1
		lastj = d.NumStates - 2
	}
	for i := 0; i < d.NumStates; i++ {
		d.R.Set(i, i, 0)
		for j := 0; j < d.NumStates; j++ {
			if i == j {
				continue
			}
			if sym && j > i {
				if j == lastj && i == lasti {
					d.R.Set(i, j, 1.0)
					d.R.Set(j, i, 1.0)
				} else {
					d.R.Set(i, j, params[cc])
					d.R.Set(j, i, params[cc])
					cc++
				}
			} else if sym == false {
				if j == lastj && i == lasti {
					d.R.Set(i, j, 1.0)
				} else {
					d.R.Set(i, j, params[cc])
					cc++
				}
			}
		}
	}
}

// ExpValue used for the matrix exponential
func (d *Model) ExpValue(iv []float64, blen float64) {
	for i, j := range iv {
		d.X.Set(i, i, math.Exp(j*blen))
	}
	return
}

// ExpValueFirstD get the first derivaite for NR
func (d *Model) ExpValueFirstD(blen float64) (x *mat.Dense) {
	x = mat.NewDense(d.NumStates, d.NumStates, nil)
	for i := 0; i < d.NumStates; i++ {
		d.X1.Set(i, i, d.EigenVals[i]*math.Exp(d.EigenVals[i]*blen))
	}
	x.Mul(d.EigenVecs, d.X1)
	x.Mul(x, d.EigenVecsI)
	return
}

// ExpValueSecondD get the second derivaite for NR
func (d *Model) ExpValueSecondD(blen float64) (x *mat.Dense) {
	x = mat.NewDense(d.NumStates, d.NumStates, nil)
	for i := 0; i < d.NumStates; i++ {
		d.X1.Set(i, i, (d.EigenVals[i]*d.EigenVals[i])*math.Exp(d.EigenVals[i]*blen))
	}
	x.Mul(d.EigenVecs, d.X1)
	x.Mul(x, d.EigenVecsI)
	return
}

// SetP use the standard spectral decom
func (d *Model) SetP(blen float64) {
	P := mat.NewDense(d.NumStates, d.NumStates, nil)
	d.ExpValue(d.EigenVals, blen)
	P.Mul(d.EigenVecs, d.X)
	P.Mul(P, d.EigenVecsI)
	// easier
	//d.Lock()
	d.Ps[blen] = P
	//d.Unlock()
}

// SetPSimple use the gonum matrixexp (seems faster)
func (d *Model) SetPSimple(blen float64) {
	P := mat.NewDense(d.NumStates, d.NumStates, nil)
	P.Scale(blen, d.Q)
	P.Exp(P)
	//d.Lock()
	d.Ps[blen] = P
	//d.Unlock()
}

// EmptyPDict save memory perhaps?
func (d *Model) EmptyPDict() {
	d.Ps = nil
	d.Ps = make(map[float64]*mat.Dense)
}

// EmptyPLDict the logged one
func (d *Model) EmptyPLDict() {
	d.PsL = nil
	d.PsL = make(map[float64]*mat.Dense)
}

// GetPMap get the Ps from the dictionary
func (d *Model) GetPMap(blen float64) *mat.Dense {
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

// GetPMapLogged get the Ps from the dictionary
func (d *Model) GetPMapLogged(blen float64) *mat.Dense {
	if _, ok := d.PsL[blen]; ok {
		return d.PsL[blen]
	}
	P := mat.NewDense(d.NumStates, d.NumStates, nil)
	P.Scale(blen, d.Q)
	P.Exp(P)
	for i := 0; i < d.NumStates; i++ {
		for j := 0; j < d.NumStates; j++ {
			P.Set(i, j, math.Log(P.At(i, j)))
		}
	}
	d.PsL[blen] = P
	return d.PsL[blen]
}

// GetPCalc calculate P matrix
func (d *Model) GetPCalc(blen float64) *mat.Dense {
	var P mat.Dense
	P.Scale(blen, d.Q)
	P.Exp(&P)
	//d.ExpValue(d.EigenVals, blen)
	//P.Mul(&d.EigenVecs, d.X)
	//P.Mul(&P, d.EigenVecsT)
	return &P
}

//SetMapDNA for getting the position in the array
func (d *Model) SetMapDNA() {
	d.CharMap = make(map[string][]int)
	d.CharMap["A"] = []int{0}
	d.CharMap["C"] = []int{1}
	d.CharMap["G"] = []int{2}
	d.CharMap["T"] = []int{3}
	d.CharMap["-"] = []int{0, 1, 2, 3}
	d.CharMap["N"] = []int{0, 1, 2, 3}
	d.CharMap["R"] = []int{0, 2}
	d.CharMap["Y"] = []int{1, 3}
	d.CharMap["M"] = []int{0, 1}
	d.CharMap["K"] = []int{2, 3}
	d.CharMap["S"] = []int{1, 2}
	d.CharMap["W"] = []int{0, 3}
	d.CharMap["H"] = []int{0, 1, 3}
	d.CharMap["B"] = []int{1, 2, 3}
	d.CharMap["V"] = []int{0, 1, 2}
	d.CharMap["D"] = []int{0, 2, 3}
}

// SetMapProt for getting the position in the array
func (d *Model) SetMapProt() {
	d.CharMap = make(map[string][]int)
	d.CharMap["A"] = []int{0}
	d.CharMap["R"] = []int{1}
	d.CharMap["N"] = []int{2}
	d.CharMap["D"] = []int{3}
	d.CharMap["C"] = []int{4}
	d.CharMap["Q"] = []int{5}
	d.CharMap["E"] = []int{6}
	d.CharMap["G"] = []int{7}
	d.CharMap["H"] = []int{8}
	d.CharMap["I"] = []int{9}
	d.CharMap["L"] = []int{10}
	d.CharMap["K"] = []int{11}
	d.CharMap["M"] = []int{12}
	d.CharMap["F"] = []int{13}
	d.CharMap["P"] = []int{14}
	d.CharMap["S"] = []int{15}
	d.CharMap["T"] = []int{16}
	d.CharMap["W"] = []int{17}
	d.CharMap["Y"] = []int{18}
	d.CharMap["V"] = []int{19}
	d.CharMap["-"] = []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}
	d.CharMap["X"] = []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}
	d.CharMap["B"] = []int{2, 3}
	d.CharMap["Z"] = []int{5, 6}
}

//SetMapMult for getting the position in the array
func (d *Model) SetMapMult() {
	d.CharMap = make(map[string][]int)
	d.CharMap["-"] = make([]int, d.NumStates)
	d.CharMap["N"] = make([]int, d.NumStates)
	for i := 0; i < d.NumStates; i++ {
		d.CharMap[strconv.Itoa(i)] = []int{i}
		d.CharMap["-"][i] = i
		d.CharMap["N"][i] = i
	}
}

// GetCharMap get the int map for states with ambiguities
func (d *Model) GetCharMap() map[string][]int {
	return d.CharMap
}

// GetNumStates return the number of states
func (d *Model) GetNumStates() int {
	return d.NumStates
}

// GetBF return the base frequencies
func (d *Model) GetBF() []float64 {
	return d.BF
}

// GetStochMapMatrices return matrices for stochastic mapping
func (d *Model) GetStochMapMatrices(dur float64, from int, to int) (summed *mat.Dense, summedR *mat.Dense) {
	nstates := d.NumStates
	d.DecomposeQ()
	Ql := mat.NewDense(nstates, nstates, nil)
	Ql.Zero()
	Ql.Set(from, to, d.Q.At(from, to))
	W := mat.NewDense(nstates, nstates, nil)
	W.Zero()
	W.Set(from, from, 1.)
	summed = mat.NewDense(nstates, nstates, nil)
	summed.Zero()
	summedR = mat.NewDense(nstates, nstates, nil)
	summedR.Zero()
	for i := 0; i < nstates; i++ {
		Ei := mat.NewDense(nstates, nstates, nil)
		Ei.Zero()
		Ei.Set(i, i, 1.)
		Si := mat.NewDense(nstates, nstates, nil)
		Si.Mul(d.EigenVecs, Ei)
		Si.Mul(Si, d.EigenVecsI)
		for j := 0; j < nstates; j++ {
			dij := (d.EigenVals[i] - d.EigenVals[j]) * dur
			Ej := mat.NewDense(nstates, nstates, nil)
			Ej.Zero()
			Ej.Set(j, j, 1.)
			Sj := mat.NewDense(nstates, nstates, nil)
			Sj.Mul(d.EigenVecs, Ej)
			Sj.Mul(Sj, d.EigenVecsI)
			Iijt := 0.
			if math.Abs(dij) > 10 {
				Iijt = (math.Exp(d.EigenVals[i]*dur) - math.Exp(d.EigenVals[j]*dur)) / (d.EigenVals[i] - d.EigenVals[j])
			} else if math.Abs(dij) < 10e-20 {
				Iijt = dur * math.Exp(d.EigenVals[j]*dur) * (1. + dij/2. + math.Pow(dij, 2.)/6. + math.Pow(dij, 3.)/24.)
			} else {
				if d.EigenVals[i] == d.EigenVals[j] {
					//					if isImag {
					//						Iijt = dur * exp(d.EigenVals[j]*dur) * (exp(dij) - 1.) / dij
					//					} else {
					Iijt = dur * math.Exp(d.EigenVals[j]*dur) * (math.Expm1((dij))) / dij //(expm1(real(dij))) / dij
					//					}
				} else {
					//					if isImag {
					//						Iijt = -dur * exp(d.EigenVals[i]*dur) * (exp(-dij) - 1.) / dij
					//					} else {
					Iijt = -dur * math.Exp(d.EigenVals[i]*dur) * (math.Expm1((-dij))) / dij //(expm1(real(-dij))) / dij
					//					}
				}
			}
			temp := mat.NewDense(nstates, nstates, nil)
			temp.Mul(Si, Ql)
			temp.Mul(temp, Sj)
			temp.Scale(Iijt, temp)
			summed.Add(summed, temp)
			temp = mat.NewDense(nstates, nstates, nil)
			temp.Mul(Si, W)
			temp.Mul(temp, Sj)
			temp.Scale(Iijt, temp)
			summedR.Add(summedR, temp)
		}
	}
	return
}
