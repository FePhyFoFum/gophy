package gophy

import (
	"strconv"
	"strings"
)

type Cluster struct {
	Sites         []int // this stores all of the sites with a preference for this cluster
	BranchLengths []float64
	LogLike       float64
	SiteWeights   map[int]float64 // this will store the probability that each site in the MATRIX belongs here.
}

func (c *Cluster) CalcLL(tree *Node) {
	c.LogLike = SubUnrootedLogLikeParallel(tree, c.Sites, 6)
}

func (c *Cluster) WriteClusterPhylip(nodes []*Node) string {
	seqs := make(map[string][]string)
	for _, n := range nodes {
		if n.Nam != "" && len(n.Chs) == 0 {
			var s []string
			for _, tr := range c.Sites {
				var con string
				if n.Mis[tr] != true {
					con = strconv.FormatFloat(n.ContData[tr], 'f', 6, 64)
				} else {
					con = "?"
				}
				s = append(s, con)
			}
			seqs[n.Nam] = s
		}
	}
	filestr := strconv.Itoa(len(seqs)) + "\t" + strconv.Itoa(len(c.Sites)) + "\n"
	for k, v := range seqs {
		filestr += k + "\t" + strings.Join(v, "\t") + "\n"
	}
	return filestr
}
