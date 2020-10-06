package clustering

import (
	"fmt"
	"math"
	"math/rand"

	"github.com/FePhyFoFum/gophy"
)

func (s *HCSearch) RunEM() {
	for i := 0; i < s.Gen; i++ {
		if i%s.PrintFreq == 0 {
			fmt.Println("ITERATION", i)
		}
		s.updateClusters()
		s.updateMixtureBranchLengths()
	}
	fmt.Println(len(s.Clusters), s.ClusterString())
}

func (s *HCSearch) SplitEM() {
	clen := len(s.Clusters)
	if clen > 1 {
		for i := 0; i < s.SplitGen; i++ {
			s.updateClusters()
			s.updateMixtureBranchLengths()
			//s.updateClassificationBranchLengths()
			//fmt.Println(i)
			//fmt.Println(s.ClusterString())
		}
		//fmt.Println(clen, s.ClusterString())
	}
}

func (s *HCSearch) updateMixtureBranchLengths() {
	for _, v := range s.Clusters {
		if len(v.Sites) == 0 {
			continue
		}
		assignClusterLengths(s.PreorderNodes, v)
		IterateLengthsWeighted(s.Tree, v, 20)
	}
}

func (s *HCSearch) updateClusters() {
	var weights map[int]float64
	for k, v := range s.SiteAssignments {
		weights = s.siteClusterUpdate(k, v)
		for l, c := range s.Clusters {
			c.SiteWeights[k] = weights[l]
			//fmt.Println(k, l, weights[l])
		}
	}
}

func (s *HCSearch) clusterLL() {
	for _, c := range s.Clusters {
		curll := 0.0
		for site := range c.Sites {
			curll += SingleSiteLL(s.Tree.Rt, site)
		}
		c.LogLike = curll
	}
}

func (s *HCSearch) siteClusterUpdate(site int, siteClusterLab int) (weights map[int]float64) {
	siteCluster := s.Clusters[siteClusterLab]
	bestLL := -100000000000.
	var bestClustLab int
	var bestClust *Cluster
	llsum := 0.0
	llmap := make(map[int]float64)
	weights = map[int]float64{}
	for k, v := range s.Clusters {
		for i, n := range s.PreorderNodes { //assign current cluster's branch lengths
			n.BMLen = v.BranchLengths[i]
		}
		var constrain float64
		constrain = (math.Log(float64(len(v.Sites))+(s.Alpha/float64(len(s.Clusters)))) / (s.NumPoints + s.Alpha - 1.))
		//constrain = 0
		curll := SingleSiteLL(s.Tree.Rt, site) + constrain
		llmap[k] = curll
		llsum += curll

		if curll > bestLL {
			bestLL = curll
			bestClustLab = k
			bestClust = v
		}
	}
	for k, v := range llmap {
		weights[k] = v / llsum
	}
	if bestClustLab != siteClusterLab { //move the site to its new cluster if the best cluster has changed
		var newSlice []int
		for _, s := range siteCluster.Sites {
			if s != site {
				newSlice = append(newSlice, s)
			}
		}
		siteCluster.Sites = newSlice
		bestClust.Sites = append(bestClust.Sites, site)
		s.SiteAssignments[site] = bestClustLab
	}
	return
}

func InitEMSearch(tree *gophy.Tree, gen int, k int, pr int, alpha float64) *HCSearch {
	s := new(HCSearch)
	s.Tree = tree
	s.PreorderNodes = tree.Pre
	s.Gen = gen
	s.K = k
	//s.startingClustersEMOnly()
	s.singleStartingCluster()
	s.perturbAndUpdate(3)
	s.PrintFreq = pr
	s.Alpha = alpha
	s.ExpandPenalty = math.Log(s.Alpha / (s.Alpha + s.NumPoints))
	return s
}

/*
this gives starting clusters when K is unknown (for basically a prior-free DPP-style mixture model)
func (search *HCSearch) startingClusters() {
	clus := make(map[int]*Cluster)
	lab := 0
	siteClust := make(map[int]int)
	for k := range search.Tree.ContData {
		cur := new(Cluster)
		cur.Sites = append(cur.Sites, k)
		ClusterMissingTraitsEM(search.Tree, cur, 10)
		clus[lab] = cur
		siteClust[k] = lab
		lab++
	}
	search.Clusters = clus
	search.SiteAssignments = siteClust
}
*/

func (search *HCSearch) startingClustersEMOnly() {
	clus := make(map[int]*Cluster)
	siteClust := make(map[int]int)
	var clustLabs []int
	for i := 0; i < search.K; i++ { //create clusters
		cur := new(Cluster)
		clus[i] = cur
		clustLabs = append(clustLabs, i)
		cur.SiteWeights = map[int]float64{}
	}
	for k := range search.Tree.Rt.ContData {
		lab := rand.Intn(search.K)
		/*
			var lab int
			if k < 50 {
				lab = 0
			} else {
				lab = 1
			}
		*/
		cur := clus[lab]
		cur.Sites = append(cur.Sites, k)
		siteClust[k] = lab
	}
	stWt := 1.0 / float64(search.K)
	search.Clusters = clus
	search.SiteAssignments = siteClust
	for _, cur := range search.Clusters {
		if len(cur.Sites) == 0 {
			for range search.PreorderNodes {
				cur.BranchLengths = append(cur.BranchLengths, rand.Float64())
			}
			continue
		}

		for i := range search.SiteAssignments {
			cur.SiteWeights[i] = stWt
		}
		//ClusterMissingTraitsEM(search.Tree, cur, 10)
		IterateLengthsWeighted(search.Tree, cur, 40)
	}

}
