package clustering

import (
	"math"
)

type SiteConfiguration struct {
	Sites         map[int]map[int]bool
	AIC           float64
	ClusterTrees  map[int]string
	ClusterSizes  map[int]int
	ClusterString string
}

func (c *SiteConfiguration) CalcClusterSizes() {
	var sizes = map[int]int{}
	for _, v := range c.Sites {
		size := len(v)
		if _, ok := sizes[size]; !ok {
			sizes[size] = 1
		} else {
			sizes[size]++
		}
	}
	c.ClusterSizes = sizes
}

func (c *SiteConfiguration) Equals(check *SiteConfiguration) (equal bool) {
	equal = false
	diff := math.Abs(c.AIC - check.AIC)
	if diff > 0.5 {
		return
	} else if diff < 0.01 {
		equal = true
		return
	}
	if len(c.ClusterTrees) == 1 && len(check.ClusterTrees) == 1 {
		equal = true
		return
	}
	if len(c.Sites) != len(check.Sites) {
		return
	}
	if len(c.ClusterSizes) == 0 {
		c.CalcClusterSizes()
	}
	if len(check.ClusterSizes) == 0 {
		check.CalcClusterSizes()
	}
	for k, v := range c.ClusterSizes {
		if _, ok := check.ClusterSizes[k]; !ok {
			equal = false
			return
		}
		if check.ClusterSizes[k] != v {
			equal = false
			return
		}
	}
	allmatch := true
	for _, v := range check.Sites {
		match := false
		for _, v2 := range c.Sites {
			if len(v) != len(v2) {
				continue
			}
			elementnotfound := false
			for m := range v {
				if _, ok := v2[m]; !ok {
					elementnotfound = true
					break
				}
			}
			if elementnotfound == false {
				match = true
				break
			}
		}
		if match == false {
			//v != v2
			allmatch = false
			break
		}
	}
	equal = allmatch
	return
}
