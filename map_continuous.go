package gophy

import (
	"fmt"
	"os"

	//"io/ioutil"
	"strconv"
	"strings"
)

//readContinuous will read in a tab-separated file of continuous traits and return a map
func readContinuous(traitfl string) map[string][]float64 {
	traitlines := ReadLine(traitfl)
	tm := make(map[string][]float64)
	var err error
	var curtr float64
	ntraits := 0
	for i, line := range traitlines {
		if line == "" {
			continue
		}
		ss := strings.Split(line, "\t")
		if i == 0 {
			ntraits = len(ss[1:])
		}
		curtax := ss[0]
		var curtraits []float64
		if len(ss[1:]) != ntraits {
			fmt.Println("not all of the taxa in this file share the same number of traits. problem caught at line ", i+1)
		}
		for _, v := range ss[1:] {
			if v != "?" {
				curtr, err = strconv.ParseFloat(v, 64)
				if err != nil {
					fmt.Println("couldn't parse a trait to a float starting at line ", i+1)
					fmt.Println(curtax, v)
				}
			} else {
				curtr = -1000000.0 //this just initializes a missing value in an obvious way (want to keep a float for data imputation later)
			}
			curtraits = append(curtraits, curtr)
		}
		tm[curtax] = curtraits
	}
	return tm
}

//makeMissingDataSlice will intialize the Mis attribute of node t by identifying traits with LARGE (ie. missing) values and updating Mis accordingly.
func makeMissingDataSlice(t *Node) {
	for _, tr := range t.ContData { //create slice marking missing data
		if tr == -1000000.0 {
			t.Mis = append(t.Mis, true)
		} else {
			t.Mis = append(t.Mis, false)
		}
	}
}

//MapContinuous maps the traits from a file to the tips of a tree and initializes slices of the same length for the internal nodes
func MapContinuous(t *Tree, traitfl string) {
	traits := readContinuous(traitfl)
	var z float64
	z = 0.0
	ntraits := 0
	for _, n := range t.Post {
		if len(n.Chs) == 0 {
			if _, ok := traits[n.Nam]; !ok {
				fmt.Println(n.Nam, " is in the tree but missing from the trait matrix\n(L.87 in map_continuous.go)")
				os.Exit(0)
			}
			n.ContData = traits[n.Nam]
			if ntraits == 0 {
				ntraits = len(n.ContData)
			}
			makeMissingDataSlice(n)
		} else {
			for count := 0; count < ntraits; count++ {
				n.ContData = append(n.ContData, z)
				n.Mis = append(n.Mis, false)
			}
			n.ConPruneLen = make([]float64, len(n.ContData))
		}
		n.ConPruneLen = make([]float64, len(n.ContData))
	}
}
