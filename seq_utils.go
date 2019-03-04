package gophy

//GetSitePatterns return site pattens when the datatype for the alignment is a map[string]string
func GetSitePatterns(seqs map[string]string, nsites int, seqnames []string) (patterns map[string][]int,
	patternsint map[int]float64, gapsites []int, constant []int, uninformative []int) {
	patterns = make(map[string][]int)
	for k := 0; k < nsites; k++ {
		tp := ""
		As, Cs, Gs, Ts, gapcount := 0, 0, 0, 0, 0
		for _, j := range seqnames {
			tp += string(seqs[j][k])
			switch c := string(seqs[j][k]); c {
			case "A":
				As++
			case "C":
				Cs++
			case "G":
				Gs++
			case "T":
				Ts++
			case "-":
				gapcount++
			case "N":
				gapcount++
			default:
				//fmt.Println(c)
			}
		}
		efc := len(seqs) - gapcount
		if gapcount == len(seqs) {
			gapsites = append(gapsites, k)
			continue
		} else if As >= efc || Cs >= efc || Gs >= efc || Ts >= efc {
			constant = append(constant, k)
		}
		twocount := 0
		if As >= 2 {
			twocount++
		}
		if Cs >= 2 {
			twocount++
		}
		if Gs >= 2 {
			twocount++
		}
		if Ts >= 2 {
			twocount++
		}
		if twocount < 2 {
			uninformative = append(uninformative, k)
		}
		if _, ok := patterns[tp]; !ok {
			patterns[tp] = make([]int, 0)
		}
		patterns[tp] = append(patterns[tp], k)
	}
	patternsint = make(map[int]float64) // key is first site, value is the number of that one
	for _, j := range patterns {
		patternsint[j[0]] = float64(len(j))
	}
	return
}
