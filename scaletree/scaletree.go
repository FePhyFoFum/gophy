package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/FePhyFoFum/gophy"
)

//this is for scaling a tree to be ultrametric
// just splitting the difference and all that

//mrca filename should be
//name1,name2 date

func readTreeFile(treefilename string) (gophy.Tree, *gophy.Node) {
	f, err := os.Open(treefilename)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	defer f.Close()
	scanner := bufio.NewReader(f)
	fmt.Fprintln(os.Stderr, "reading trees")
	var t gophy.Tree
	var rt *gophy.Node
	for {
		ln, err := scanner.ReadString('\n')
		if len(ln) > 0 {
			rt = gophy.ReadNewickString(ln)
			t.Instantiate(rt)
			break
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			break
		}
	}
	return t, rt
}

func extractDates(treefilename string, mrcafilename string) {
	f, err := os.Create(mrcafilename)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	defer f.Close()
	t, _ := readTreeFile(treefilename)
	gophy.SetHeights(&t)
	for _, i := range t.Pre {
		if len(i.Chs) > 0 {
			//fmt.Println(i.Newick(true))
			lf := i.Chs[0].GetTipNames()[0]
			rt := i.Chs[1].GetTipNames()[0]
			f.WriteString(lf + "," + rt + " " + fmt.Sprintf("%f", i.Height) + "\n")
		}
	}
}

func main() {
	tfn := flag.String("t", "", "tree filename")
	mfn := flag.String("m", "", "mrca filename")
	ext := flag.String("e", "", "extract dates from tree")
	flag.Parse()
	if len(*tfn) == 0 {
		os.Exit(0)
	}
	if len(*mfn) == 0 && len(*ext) == 0 {
		os.Exit(0)
	}
	fmt.Fprintln(os.Stderr, "treefile:", *tfn)
	if len(*mfn) > 0 {
		fmt.Fprintln(os.Stderr, "mrcafile:", *mfn)
	} else {
		fmt.Fprintln(os.Stderr, "extract dates from:", *ext)
		*mfn = "extracted_dates.txt"
		extractDates(*ext, *mfn)
	}
	// read tree file
	t, rt := readTreeFile(*tfn)
	nmsnds := make(map[string]*gophy.Node)
	for _, i := range t.Post {
		if len(i.Nam) > 0 {
			nmsnds[i.Nam] = i
		}
		i.Len = 0.0
	}
	//end tree file reading
	//read mrca file
	f, err := os.Open(*mfn)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	defer f.Close()
	mrcas := make(map[*gophy.Node]float64) //map is node and float is date
	scanner := bufio.NewReader(f)
	fmt.Fprintln(os.Stderr, "reading mrcas")
	for {
		ln, err := scanner.ReadString('\n')
		if len(ln) > 0 {
			//fmt.Fprintln(os.Stderr, strings.Trim(ln, "\n"))
			spls1 := strings.Split(strings.Trim(ln, "\n"), " ") //mrca = 0, date = 1
			spls2 := strings.Split(spls1[0], ",")
			nds := make([]*gophy.Node, 2)
			nds[0] = nmsnds[spls2[0]]
			nds[1] = nmsnds[spls2[1]]
			ff, ferr := strconv.ParseFloat(spls1[1], 64)
			if ferr != nil {
				fmt.Fprintln(os.Stderr, "problem parsing", spls1[1], "as float64")
			}
			nd := gophy.GetMrca(nds, rt)
			if nd.Par != nil {
				for len(nd.Par.Chs) == 1 {
					nd = nd.Par
					if nd.Par == nil {
						break
					}
				}
			}
			mrcas[nd] = ff
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			break
		}
	}
	setFeasibleTimes(mrcas, rt, t)
	fmt.Println(t.Rt.Newick(true) + ";")
}

func setFeasibleTimes(mrcas map[*gophy.Node]float64,
	rt *gophy.Node, t gophy.Tree) {
	for _, i := range t.Pre {
		if _, ok := mrcas[i]; ok {
			//fmt.Println("WORKINGWITH", mrcas[i], i)
			if i.FData["max"] < mrcas[i] && i.FData["max"] != 0.0 {
				fmt.Println("PROBLEM CALIBRATION", mrcas[i], i.FData["max"], mrcas)
				os.Exit(0)
			}
			i.FData["date"] = mrcas[i]
			i.FData["max"] = mrcas[i]
			i.FData["min"] = mrcas[i]
			for _, k := range i.PreorderArray() {
				if k.FData["max"] == 0 || k.FData["max"] > mrcas[i] {
					k.FData["max"] = mrcas[i]
				}
			}
			cur := i
			for {
				if cur.FData["max"] == 0 || cur.FData["min"] < mrcas[i] {
					cur.FData["min"] = mrcas[i]
				}
				if cur.Par == nil {
					break
				} else {
					cur = cur.Par
				}
			}
			i.SData["fixed"] = "fixed"
		} else {
			i.SData["fixed"] = "free"
		}
	}
	for _, i := range t.Pre {
		//fmt.Println(i, i.FData)
		if i.SData["fixed"] == "free" {
			if len(i.Chs) == 0 {
				i.FData["date"] = 0
				i.FData["max"] = 0
				i.FData["min"] = 0
				i.SData["fixed"] = "fixed"
			} else {
				//fmt.Println(" GOOD DATE", i, i.FData["max"], i.FData["min"], i.FData["date"])
				i.FData["date"] = getEqual(i, i.FData["min"], i.FData["max"])
				for _, j := range i.Chs {
					if j.SData["fixed"] == "free" {
						if j.FData["max"] > i.FData["date"] {
							j.FData["max"] = i.FData["date"]
						}
					}
				}
			}
		}
	}
	fixDates(t)
}

//just splitting the difference
func getOrder(node *gophy.Node, minv float64, maxv float64) float64 {
	num := (maxv - minv) / 2.
	return maxv - num
}

func getEqual(node *gophy.Node, minv float64, maxv float64) float64 {
	num := (maxv - minv) / (getMaxFreeNodes(node) + 1)
	return maxv - num
}

func fixDates(t gophy.Tree) {
	for _, i := range t.Post {
		if i.Par != nil {
			i.Len = i.Par.FData["date"] - i.FData["date"]
			if i.Len < 0 {
				fmt.Println(i)
				fmt.Println(i.Len, i.Par.FData["date"], i.FData["date"])
				os.Exit(0)
			}
		}
	}
}

func getMaxFreeNodes(node *gophy.Node) float64 {
	v := 0.
	if len(node.Chs) == 0 {
		return 0.
	}
	for _, i := range node.GetTips() {
		vl := 1.
		going := true
		cur := i
		for going {
			if cur.SData["fixed"] == "fixed" {
				vl = 1.
			}
			par := cur.Par
			if par == node {
				going = false
				break
			} else {
				vl += 1.
				cur = par
			}
		}
		if vl > v {
			v = vl
		}
	}
	return v
}
