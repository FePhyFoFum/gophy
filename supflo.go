package gophy

import (
	"math"
)

//SupFlo super float
type SupFlo struct {
	mant        float64
	exp         int
	stilldouble bool
	uplimit     float64
	lowlimit    float64
}

//NewSupFlo get a new one
func NewSupFlo(m float64, e int) *SupFlo {
	sf := SupFlo{mant: m, exp: e, uplimit: 1e+175, lowlimit: 1e-175}
	if e == 0 {
		sf.stilldouble = true
	} else {
		sf.adjustDecimal()
	}
	return &sf
}

//SetFloat64 set it to a float64 value
func (s *SupFlo) SetFloat64(m float64) {
	s.exp = 0
	s.stilldouble = true
	s.mant = m
}

func isnormal(n float64) bool {
	ret := true
	if math.IsInf(n, 1) {
		return false
	} else if math.IsInf(n, -1) {
		return false
	} else if math.IsNaN(n) {
		return false
	}
	return ret
}

func (s *SupFlo) adjustDecimal() {
	if s.exp != 0 {
		s.stilldouble = false
	}
	if math.Abs(s.mant) < s.uplimit &&
		math.Abs(s.mant) > s.lowlimit {
		return
	}
	s.stilldouble = false
	if isnormal(s.mant) == false {
		s.exp = 0
		s.stilldouble = true
	} else if s.mant == 0 {
		s.exp = 0
		s.stilldouble = true
	} else {
		exa := int(math.Round(math.Log10(s.mant) - 0.5))
		s.mant = s.mant / math.Pow10(exa)
		s.exp += exa
		/*
			for math.Abs(s.mant) >= 10 {
				s.mant *= 0.1
				s.exp++
			}
			for math.Abs(s.mant) < 1 {
				s.mant *= 10.0
				s.exp--
			}*/
	}
}

// Mul s * x return result
func (s *SupFlo) Mul(x *SupFlo) *SupFlo {
	result := NewSupFlo(s.mant*x.mant, s.exp+x.exp)
	result.adjustDecimal()
	/*
		if result.stilldouble == true {
			if math.Abs(result.mant) > s.uplimit ||
				math.Abs(result.mant) < s.lowlimit {
				result.adjustDecimal()
			}
		} else {
			result.adjustDecimal()
		}*/
	return result
}

// MulFloat s * x return result and x is a float64
func (s *SupFlo) MulFloat(x float64) *SupFlo {
	result := NewSupFlo(s.mant*x, s.exp)
	result.adjustDecimal()
	return result
}

// Div s/x return sup
func (s *SupFlo) Div(x *SupFlo) *SupFlo {
	result := NewSupFlo(s.mant/x.mant, s.exp-x.exp)
	result.adjustDecimal()
	return result
}

// Add s+x return result
func (s *SupFlo) Add(x *SupFlo) *SupFlo {
	//only tricky thing is converting them to same exponent
	if s.mant != 0.0 {
		exponentdif := x.exp - s.exp
		result := NewSupFlo(s.mant+(x.mant*(math.Pow10(exponentdif))), s.exp)
		result.adjustDecimal()
		return result
	}
	result := NewSupFlo(x.mant, x.exp)
	result.adjustDecimal()
	return result
}

//Sub s-x return result
func (s *SupFlo) Sub(x *SupFlo) *SupFlo {
	//only tricky thing is converting them to same exponent
	if s.mant != 0 {
		exponentdif := x.exp - s.exp
		result := NewSupFlo(s.mant-(x.mant*(math.Pow10(exponentdif))), s.exp)
		result.adjustDecimal()
		return result
	}
	result := NewSupFlo(-1.0*x.mant, x.exp)
	result.adjustDecimal()
	return result
}

//PlusPlus ++
func (s *SupFlo) PlusPlus() {
	s.mant++
	s.adjustDecimal()
}

//MinMin --
func (s *SupFlo) MinMin() {
	s.mant--
	s.adjustDecimal()
}

//MulEq *=
func (s *SupFlo) MulEq(x *SupFlo) {
	s.mant *= x.mant
	s.exp += x.exp
	s.adjustDecimal()
}

//MulEqFloat *=
func (s *SupFlo) MulEqFloat(x float64) {
	s.mant *= x
	if s.stilldouble == true {
		if math.Abs(s.mant) > s.uplimit ||
			math.Abs(s.mant) < s.lowlimit {
			s.adjustDecimal()
		}
	} else {
		s.adjustDecimal()
	}
}

//DivEq /=
func (s *SupFlo) DivEq(x *SupFlo) {
	s.mant /= x.mant
	s.exp -= x.exp
	s.adjustDecimal()
}

//AddEq +=
func (s *SupFlo) AddEq(x *SupFlo) {
	//only tricky thing is converting them to same exponent
	if s.mant != 0 {
		if s.stilldouble == true && x.stilldouble == true {
			s.mant += x.mant
			if math.Abs(s.mant) > s.uplimit || math.Abs(s.mant) < s.lowlimit {
				s.adjustDecimal()
			}
		} else {
			s.mant += (x.mant * (math.Pow10(x.exp - s.exp)))
			s.adjustDecimal()
		}
	} else {
		s.mant = x.mant
		s.exp = x.exp
		if s.stilldouble == false || x.stilldouble == false {
			s.adjustDecimal()
		}
	}
}

//SubEq -=
func (s *SupFlo) SubEq(x SupFlo) {
	//only tricky thing is converting them to same exponent
	if s.mant != 0 {
		exponentdif := x.exp - s.exp
		s.mant = s.mant - (x.mant * (math.Pow10(exponentdif)))
		s.adjustDecimal()
	} else {
		s.mant = -1.0 * x.mant
		s.exp = x.exp
		s.adjustDecimal()
	}
}

//Less <
func (s SupFlo) Less(x SupFlo) bool {
	if s.exp < x.exp {
		return true
	} else if s.exp == x.exp && s.mant < x.mant {
		return true
	} else {
		return false
	}
}

//Greater >
func (s SupFlo) Greater(x SupFlo) bool {
	if s.exp > x.exp {
		return true
	} else if s.exp == x.exp && s.mant > x.mant {
		return true
	} else {
		return false
	}
}

//LessEq <=
func (s SupFlo) LessEq(x SupFlo) bool {
	if s.exp < x.exp {
		return true
	} else if s.exp == x.exp && s.mant <= x.mant {
		return true
	} else {
		return false
	}
}

//GreaterEq >=
func (s SupFlo) GreaterEq(x SupFlo) bool {
	if s.exp > x.exp {
		return true
	} else if s.exp == x.exp && s.mant >= x.mant {
		return true
	} else {
		return false
	}
}

//GetExp just get it
func (s *SupFlo) GetExp() int {
	return s.exp
}

//GetMant just get it
func (s *SupFlo) GetMant() float64 {
	return s.mant
}

//Float64 just get it
func (s *SupFlo) Float64() float64 {
	return s.mant * math.Pow10(s.exp)
}

//GetLn ln
func (s *SupFlo) GetLn() *SupFlo {
	//ln(a * 10^b) = ln(a) + ln(10^b) = ln(a) + log10 (10^b) / log10 (e^1) = ln(a) + b/log10(e^1)
	result := NewSupFlo(math.Log(s.mant)+(1.0*float64(s.exp))/math.Log10(math.Exp(1)), 0)
	result.adjustDecimal()
	return result
}

//Abs absolute value
func (s *SupFlo) Abs() *SupFlo {
	if s.mant < 0 {
		result := NewSupFlo(-s.mant, s.exp)
		return result
	} else {
		result := NewSupFlo(s.mant, s.exp)
		return result
	}
}

//SwitchSign change the sign
func (s *SupFlo) SwitchSign() {
	s.mant = -1 * s.mant
}
