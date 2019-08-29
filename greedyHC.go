package gophy

import (
	"bufio"
	"bytes"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"strconv"
)

type HCSearch struct {
	Tree            *Tree
	PreorderNodes   []*Node
	Clusters        map[int]*Cluster
	SiteAssignments map[int]int
	Gen             int
	Threads         int
	Workers         int
	RunName         string
	LogOutFile      string
	K               int
	PrintFreq       int
	CurrentAIC      float64
	NumTraits       float64
	Criterion       int
	SavedConfig     []*SiteConfiguration
	CurBestAIC      float64
	JoinLikes       map[int]map[int]float64
	SplitGen        int
	Alpha           float64
	NumPoints       float64
	ExpandPenalty   float64
	MinK            int
}

func (s *HCSearch) NewSiteConfig() *SiteConfiguration {
	config := new(SiteConfiguration)
	var sitemap = map[int]map[int]bool{}
	var treemap = map[int]string{}
	count := 0
	for _, v := range s.Clusters {
		sitemap[count] = map[int]bool{}
		for _, site := range v.Sites {
			sitemap[count][site] = true
		}
		for i, node := range s.PreorderNodes {
			node.BMLen = v.BranchLengths[i]
		}
		treemap[count] = s.Tree.Rt.BMPhylogram()
		count++
	}
	config.AIC = s.CurrentAIC
	//fmt.Println(config.AIC)
	config.Sites = sitemap
	config.ClusterTrees = treemap
	config.ClusterString = s.ClusterString()
	return config
}

func (s *HCSearch) RunSingleHC() {
	f, err := os.Create(s.RunName)
	if err != nil {
		log.Fatal(err)
	}
	w := bufio.NewWriter(f)

	count := 0
	for {
		clcount := len(s.Clusters)
		fmt.Println(clcount)
		if clcount == 1 {
			break
		}
		quit := s.bestClusterJoin()
		if quit == true {
			break
		}
		count++
		if count > 1000000 {
			fmt.Println("couldn't find clusters to join")
			os.Exit(0)
		}
	}
	fmt.Println(s.ClusterString())
	for _, c := range s.Clusters {
		for i, n := range s.PreorderNodes {
			n.BMLen = c.BranchLengths[i]
		}
		fmt.Fprint(w, s.Tree.Rt.BMPhylogram()+";"+"\n")
	}
	err = w.Flush()
	if err != nil {
		log.Fatal(err)
	}
	f.Close()
}

func (s *HCSearch) CalcRelLikes() (denom float64) {
	denom = 0.
	var deltaAICi float64
	for _, c := range s.SavedConfig {
		if c.AIC == s.CurBestAIC {
			denom += 1.0
			continue
		}
		deltaAICi = c.AIC - s.CurBestAIC
		denom += math.Exp(-0.5 * deltaAICi)
	}
	return
}

func (s *HCSearch) CheckCluster(checkConfig *SiteConfiguration) (keep bool) {
	keep = true
	denom := s.CalcRelLikes()
	relLike := math.Exp(-0.5*checkConfig.AIC - s.CurBestAIC)
	denom += relLike
	weight := relLike / denom
	if weight < .01 {
		keep = false
		return
	}
	for _, c := range s.SavedConfig {
		seen := checkConfig.Equals(c)
		if seen == true {
			keep = false
			return
		}

	}
	return
}

func (s *HCSearch) PerturbedRun() {
	count := 0
	var bestClust string
	bestAIC := 1000000000.
	var bestK int
	fmt.Println("iteration\tcurrent score\tbest score\tbest K")
	for {
		clcount := len(s.Clusters)
		//fmt.Println(clcount, s.CurrentAIC)
		var quit bool
		if clcount == s.MinK {
			quit = true
		} else {
			//fmt.Println(len(s.Clusters), s.CurrentAIC)
			quit = s.bestClusterJoin()
		}
		if quit == true {
			count++
			if clcount >= s.MinK {
				config := s.NewSiteConfig()
				var keep bool
				if s.CurrentAIC < bestAIC {
					bestAIC = s.CurrentAIC
					bestK = clcount
					s.CurBestAIC = bestAIC
					bestClust = s.ClusterString()
					keep = true
				} else {
					keep = s.CheckCluster(config)
				}
				if keep == true {
					s.SavedConfig = append(s.SavedConfig, config)
				}
			}
			if count > s.Gen {
				break
			}
			fmt.Println(count, s.CurrentAIC, bestAIC, bestK)
			s.perturbClusters()
			/*for _, c := range s.Clusters {
				assignClusterLengths(s.Tree.Pre, c)
				fmt.Println(s.Tree.Rt.BMPhylogram())
			}*/
			s.JoinLikes = make(map[int]map[int]float64)
		}
		if count > 10000000000 {
			fmt.Println("couldn't find clusters to join")
			os.Exit(0)
		}
	}
	s.RefineSavedClusterings()
	//fmt.Println(len(s.SavedConfig))
	s.WriteBestClusters()
	s.WriteClusterTrees()
	fmt.Println(bestAIC, bestClust)
}

func (s *HCSearch) classLogLike() (err bool) {
	for _, c := range s.Clusters {
		if len(c.BranchLengths) != len(s.PreorderNodes) {
			err = true
			return
		}
		assignClusterLengths(s.PreorderNodes, c)
		c.CalcLL(s.Tree.Rt)
	}
	err = false
	return
}

func (s *HCSearch) WriteBestClusters() {
	f, err := os.Create(s.RunName + "_bestClusters")
	if err != nil {
		log.Fatal(err)
	}
	w := bufio.NewWriter(f)
	for _, c := range s.SavedConfig {
		fmt.Fprint(w, strconv.FormatFloat(c.AIC, 'f', 6, 64)+"\t"+c.ClusterString+"\n")
	}
	err = w.Flush()
	if err != nil {
		log.Fatal(err)
	}
	f.Close()
}

func (s *HCSearch) WriteClusterTrees() {
	for i, c := range s.SavedConfig {
		stringCount := strconv.Itoa(i)
		f, err := os.Create(s.RunName + "_config" + stringCount + "TREES")
		if err != nil {
			log.Fatal(err)
		}
		w := bufio.NewWriter(f)
		for lab, sitels := range c.Sites {
			var buffer bytes.Buffer
			buffer.WriteString("| ")
			//fmt.Println(sitels)
			for site := range sitels {
				cur := strconv.Itoa(site)
				buffer.WriteString(cur)
				buffer.WriteString(" ")
			}
			buffer.WriteString("| \n\n")
			fmt.Fprint(w, buffer.String())
			fmt.Fprint(w, c.ClusterTrees[lab])
			fmt.Fprint(w, ";\n\n")
		}

		err = w.Flush()
		if err != nil {
			log.Fatal(err)
		}
		f.Close()
	}
}

func (s *HCSearch) RefineSavedClusterings() {
	var kept []*SiteConfiguration
	for _, c := range s.SavedConfig {
		if c.AIC == s.CurBestAIC {
			kept = append(kept, c)
			continue
		}
		denom := s.CalcRelLikes()
		relLike := math.Exp(-0.5*c.AIC - s.CurBestAIC)
		weight := relLike / denom
		if weight > 0.01 {
			kept = append(kept, c)
		}
	}
	s.SavedConfig = kept
}

func (s *HCSearch) checkAndAddK() {
	if len(s.Clusters) < s.K {
		biggest := MaxClustLab(s.Clusters)
		for {
			biggest++
			newC := new(Cluster)
			s.Clusters[biggest] = newC
			for range s.PreorderNodes {
				newC.BranchLengths = append(newC.BranchLengths, rand.Float64())
			}
			newC.SiteWeights = map[int]float64{}
			if len(s.Clusters) == s.K {
				break
			}
		}
	}
}

func assignClusterLengths(ns []*Node, c *Cluster) {
	for i, n := range ns {
		n.BMLen = c.BranchLengths[i]
	}
}

func (s *HCSearch) calcClusterSiteWeights() {
	ll := 0.
	llsum := 0.
	cll := make(map[int]float64)
	for site := range s.SiteAssignments {
		llsum = 0.
		for lab, c := range s.Clusters {
			assignClusterLengths(s.PreorderNodes, c)
			var constrain float64
			if len(c.Sites) > 0 {
				//constrain = math.Log(float64(len(v.Sites)) / (float64(len(s.Tree.ContData)))+0.5)
				constrain = math.Log(float64(len(c.Sites)) / (s.NumPoints + s.Alpha))
			} else {
				//constrain = 0
				constrain = s.ExpandPenalty
			}
			ll = SingleSiteLL(s.Tree.Rt, site) + constrain
			llsum += ll
			cll[lab] = ll
			//fmt.Println(c.SiteWeights)
		}
		for k, v := range cll {
			s.Clusters[k].SiteWeights[site] = v / llsum
		}
	}
}

func (s *HCSearch) updateClassificationBranchLengths() {
	for _, v := range s.Clusters {
		if len(v.Sites) == 0 {
			continue
		}
		for i, n := range s.PreorderNodes { //assign current cluster's branch lengths
			n.BMLen = v.BranchLengths[i]
		}
		ClusterMissingTraitsEM(s.Tree, v, 20)
	}
}

func (s *HCSearch) perturbAndUpdate(it int) {
	for i := 0; i < it; i++ {
		s.perturbSites()
		s.updateClassificationBranchLengths()
	}
	s.calcClusterSiteWeights()
}

func (s *HCSearch) perturbClusters() {
	s.perturbAndUpdate(3)
	s.checkAndAddK()
	//s.calcClusterSiteWeights()
	s.SplitEM()
	s.removeEmptyK()
	//for _, c := range s.Clusters {
	//fmt.Println(c.Sites)
	//}
	//s.updateClassificationBranchLengths()
	err := s.classLogLike()
	if err == true {
		fmt.Println("There was a problem dealing with the branch length clusters encountered in perturbClusters()")
		os.Exit(1)
	}
	s.CurrentAIC = s.calcAIC()
}

func (s *HCSearch) removeEmptyK() {
	for k, v := range s.Clusters {
		if len(v.Sites) == 0 {
			delete(s.Clusters, k)
		}
	}
}

func (s *HCSearch) perturbSites() {
	/*for {
		nclust := len(s.Clusters)
		if nclust == s.K {
			break
		}
		s.searchExpand()
	}*/
	var weights map[int]float64
	for k, v := range s.SiteAssignments {
		weights = s.reseatSites(k, v)
		for l, c := range s.Clusters {
			c.SiteWeights[k] = weights[l]
		}
	}
}

func (s *HCSearch) randomizeBranchLengths() {
	for _, n := range s.PreorderNodes[1:] {
		r := rand.Float64()
		n.BMLen = r
	}
}

func (s *HCSearch) searchExpand() {
	bestAlone := -10000000000.
	var aloneBL []float64
	var expSite int
	for site, curlab := range s.SiteAssignments {
		if len(s.Clusters[curlab].Sites) == 1 {
			continue
		}
		bestLL := -1000000000000.
		for _, v := range s.Clusters {
			assignClusterLengths(s.PreorderNodes, v)
			curll := SingleSiteLL(s.Tree.Rt, site)
			if curll > bestLL {
				bestLL = curll
			}
		}
		selfLL, selfBL := s.checkSelfLL(site, bestLL)
		if selfLL > bestAlone {
			bestAlone = selfLL
			expSite = site
			aloneBL = selfBL
		}
	}
	s.expandSingleSite(expSite, bestAlone, aloneBL)
}

func (s *HCSearch) expandSingleSite(site int, selfLL float64, aloneBL []float64) {
	var oldsites []int
	oldclust := s.Clusters[s.SiteAssignments[site]]
	for _, i := range oldclust.Sites {
		if i != site {
			oldsites = append(oldsites, i)
		}
	}
	oldclust.Sites = oldsites
	newLab := MaxClustLab(s.Clusters) + 1
	s.SiteAssignments[site] = newLab
	selfClust := new(Cluster)
	selfClust.SiteWeights = make(map[int]float64)
	selfClust.LogLike = selfLL //valid because this is a single site cluster
	selfClust.BranchLengths = aloneBL
	selfClust.Sites = append(selfClust.Sites, site)
	s.Clusters[newLab] = selfClust
}

func (s *HCSearch) checkSelfLL(site int, bestLL float64) (selfLL float64, brlens []float64) {
	var sendsites []int
	sendsites = append(sendsites, site)
	s.randomizeBranchLengths()
	GreedyIterateLengthsMissing(s.Tree, sendsites, 10)
	selfLL = SingleSiteLL(s.Tree.Rt, site)
	for _, n := range s.PreorderNodes {
		brlens = append(brlens, n.BMLen)
	}
	return
}

func (s *HCSearch) reseatSites(site int, siteClusterLab int) (weights map[int]float64) {
	siteCluster := s.Clusters[siteClusterLab]
	bestLL := -1000000000000.
	var bestClustLab int
	var bestClust *Cluster
	nclust := len(s.Clusters)
	llsum := 0.0
	llmap := make(map[int]float64)
	for k, v := range s.Clusters {
		assignClusterLengths(s.PreorderNodes, v)
		var constrain float64
		if len(v.Sites) > 0 {
			constrain = math.Log(float64(len(v.Sites)) / (s.NumPoints + s.Alpha))
		} else {
			constrain = s.ExpandPenalty
		}
		curll := SingleSiteLL(s.Tree.Rt, site) + constrain
		llmap[k] = curll
		llsum += curll
		if curll > bestLL {
			bestLL = curll
			bestClustLab = k
			bestClust = v
		}
	}
	//newself := false
	if len(siteCluster.Sites) != 1 && nclust < s.K {
		sendsites := []int{site}
		s.randomizeBranchLengths()
		GreedyIterateLengthsMissing(s.Tree, sendsites, 10)
		selfLL := SingleSiteLL(s.Tree.Rt, site) + s.ExpandPenalty
		if selfLL > bestLL {
			newLab := MaxClustLab(s.Clusters) + 1
			llmap[newLab] = selfLL
			llsum += selfLL
			bestLL = selfLL
			bestClustLab = newLab
			selfClust := new(Cluster)
			selfClust.SiteWeights = make(map[int]float64)
			selfClust.LogLike = selfLL //valid because this is a single site cluster
			for _, n := range s.PreorderNodes {
				selfClust.BranchLengths = append(selfClust.BranchLengths, n.BMLen)
			}
			s.Clusters[newLab] = selfClust
			//fmt.Println(s.Clusters[newLab])
			bestClust = selfClust
			//newself = true
		}
	}
	weights = map[int]float64{}
	for k, v := range llmap {
		weights[k] = v / llsum
	}
	if bestClustLab != siteClusterLab { //move the site to its new cluster if the best cluster has changed
		if len(siteCluster.Sites) == 1 { // delete the cluster containing the current site if it is a singleton
			delete(s.Clusters, siteClusterLab)
		} else { // delete the current site from its existing cluster if it shares the cluster with other sites
			var newSlice []int
			for _, s := range siteCluster.Sites {
				if s != site {
					newSlice = append(newSlice, s)
				}
			}
			siteCluster.Sites = newSlice
		}
		//fmt.Println("BEST", bestClust)
		bestClust.Sites = append(bestClust.Sites, site)
		s.SiteAssignments[site] = bestClustLab
		//if newself == false {
		//	assignClusterLengths(s.PreorderNodes, bestClust)
		//	bestClust.CalcLL(s.Tree)
		//}
	}
	return
}

func (s *HCSearch) calcAIC() (aic float64) {
	params := 0.
	blCount := float64(len(s.PreorderNodes)) - 1. //subtract 1 because you don't estimate the root node length
	ll := 0.
	for _, c := range s.Clusters {
		params += blCount
		ll += c.LogLike
	}
	if s.Criterion == 0 {
		aic = (2. * params) - (2. * ll)
	} else if s.Criterion == 1 {
		aic = (s.NumTraits * params) - (2. * ll)
	} else if s.Criterion == 2 {
		aic = (2. * params) - (2. * ll)
		aic -= ((2. * params) * (params + 1.)) / (s.NumTraits - params - 2.)
	}
	//aic += (((2. * (params * params)) + (2. * params)) / (s.NumTraits - params - 1.))
	return
}

func (s *HCSearch) unjoinedLikeParams(exclude1, exclude2 int) (params, ll float64) {
	params = 0.
	blCount := float64(len(s.PreorderNodes)) - 1. //subtract 1 because you don't estimate the root node length
	ll = 0.
	for i, c := range s.Clusters {
		if i != exclude1 && i != exclude2 {
			params += blCount
			ll += c.LogLike
		}
	}
	return
}

func (s *HCSearch) bestClusterJoin() (quit bool) {
	blCount := float64(len(s.PreorderNodes)) - 1. //subtract 1 because you don't estimate the root node length
	bestAIC := 1000000000000.
	//var best1, best2 *Cluster
	var bestll float64
	var bestBranchLengths []float64
	var bestClusterSites []int
	var deleteK int
	var newK int
	clen := len(s.Clusters)
	for i, c1 := range s.Clusters {
		for j, c2 := range s.Clusters {
			if i >= j {
				continue
			}
			var proposedSites []int
			for _, site := range c1.Sites {
				proposedSites = append(proposedSites, site)
			}
			for _, site := range c2.Sites {
				proposedSites = append(proposedSites, site)
			}
			params, ll := s.unjoinedLikeParams(i, j)
			//fmt.Println(params, len(s.Clusters), len(s.PreorderNodes)-1)
			precalc := true
			var clustll float64
			if _, ok := s.JoinLikes[i][j]; !ok {
				if _, ok := s.JoinLikes[i]; !ok {
					s.JoinLikes[i] = map[int]float64{}
				}
				if _, ok := s.JoinLikes[j][i]; !ok {
					s.JoinLikes[j] = map[int]float64{}
					precalc = false // pairwise calc needs to be done
				} else {
					clustll = s.JoinLikes[j][i]
				}
			} else {
				clustll = s.JoinLikes[i][j]
			}
			if precalc == false {
				for _, n := range s.PreorderNodes[1:] {
					r := rand.Float64()
					n.BMLen = r
				}
				if clen <= 15 {
					GreedyIterateLengthsMissing(s.Tree, proposedSites, 20)
				} else if clen <= 25 && clen > 15 {
					GreedyIterateLengthsMissing(s.Tree, proposedSites, 20)
				} else if clen > 25 {
					GreedyIterateLengthsMissing(s.Tree, proposedSites, 10)
				}
				clustll = SubUnrootedLogLikeParallel(s.Tree.Rt, proposedSites, 4)
				s.JoinLikes[i][j] = clustll
			}
			ll += clustll
			params += blCount
			aic := 0.
			if s.Criterion == 0 {
				aic = (2. * params) - (2. * ll)
			} else if s.Criterion == 1 {
				aic = (s.NumTraits * params) - (2. * ll)
				//fmt.Println(aic)
			} else if s.Criterion == 2 {
				aic = (2. * params) - (2. * ll)
				aic -= ((2. * params) * (params + 1.)) / (s.NumTraits - params - 2.)
			}
			if aic < bestAIC {
				bestAIC = aic
				deleteK = j
				newK = i
				bestll = clustll
				bestClusterSites = proposedSites
				var bl []float64
				for _, n := range s.PreorderNodes {
					bl = append(bl, n.BMLen)
				}
				bestBranchLengths = bl
			}
		}
	}
	if bestAIC < s.CurrentAIC { //+10. {
		newLab := MaxClustLab(s.Clusters) + 1
		addClust := new(Cluster)
		addClust.SiteWeights = make(map[int]float64)
		addClust.Sites = bestClusterSites
		for _, site := range bestClusterSites {
			s.SiteAssignments[site] = newLab
		}
		addClust.LogLike = bestll
		addClust.BranchLengths = bestBranchLengths
		s.CurrentAIC = bestAIC
		delete(s.Clusters, deleteK)
		delete(s.Clusters, newK)
		s.Clusters[newLab] = addClust
		quit = false
	} else {
		quit = true
	}
	return
}

func InitGreedyHC(tree *Tree, gen int, pr int, crit int, rstart bool, k int, runName string, splitgen int, alpha float64, minK int) *HCSearch {
	s := new(HCSearch)
	s.Tree = tree
	s.RunName = runName
	s.PreorderNodes = tree.Pre
	s.Gen = gen
	s.Criterion = crit
	s.K = k
	if rstart == false {
		s.startingClustersAllSeparate()
	} else if rstart == true {
		s.singleStartingCluster()
		//s.randomStartingClusters()
	}

	s.PrintFreq = pr
	s.JoinLikes = make(map[int]map[int]float64)
	s.SplitGen = splitgen
	s.NumPoints = float64(len(s.Tree.Rt.ContData))
	s.Alpha = alpha
	s.ExpandPenalty = math.Log(s.Alpha / (s.Alpha + s.NumPoints))
	s.MinK = minK
	return s
}

func TransferGreedyHC(tree *Tree, gen int, pr int, crit int, clus map[int]*Cluster, siteAssign map[int]int, runName string, splitgen int, alpha float64) *HCSearch {
	s := new(HCSearch)
	s.Tree = tree
	s.RunName = runName
	s.PreorderNodes = tree.Pre
	s.Gen = gen
	s.Criterion = crit
	s.PrintFreq = pr
	s.Clusters = clus
	for k, v := range s.Clusters {
		if len(v.Sites) == 0 {
			delete(s.Clusters, k) //delete empty clusters
		}
	}
	s.SiteAssignments = siteAssign
	s.NumTraits = math.Log(float64(len(clus))) //* float64(tipcount)
	s.CurrentAIC = s.calcAIC()
	s.JoinLikes = make(map[int]map[int]float64)
	s.SplitGen = splitgen
	s.NumPoints = float64(len(s.Tree.Rt.ContData))
	s.Alpha = alpha
	s.ExpandPenalty = math.Log(s.Alpha / (s.Alpha + s.NumPoints))
	return s
}

//this gives starting clusters when K is unknown
func (search *HCSearch) startingClustersAllSeparate() {
	clus := make(map[int]*Cluster)
	lab := 0
	siteClust := make(map[int]int)
	for k := range search.Tree.Rt.ContData {
		cur := new(Cluster)
		cur.Sites = append(cur.Sites, k)
		ClusterMissingTraitsEM(search.Tree, cur, 10)
		cur.LogLike = SingleSiteLL(search.Tree.Rt, k)
		clus[lab] = cur
		siteClust[k] = lab
		lab++

	}
	search.Clusters = clus
	search.SiteAssignments = siteClust
	tipcount := 0
	for _, n := range search.PreorderNodes {
		if len(n.Chs) == 0 {
			tipcount++
		}
	}
	search.NumTraits = math.Log(float64(len(clus))) * float64(tipcount)
	search.CurrentAIC = search.calcAIC()
}

//ClusterString will return a string of the current set of clusters
func (s *HCSearch) ClusterString() string {
	var buffer bytes.Buffer
	cSet := s.Clusters
	for _, like := range cSet {
		buffer.WriteString("(")
		for ind, site := range like.Sites {
			cur := strconv.Itoa(site)
			buffer.WriteString(cur)
			stop := len(like.Sites) - 1
			if ind != stop {
				buffer.WriteString(",")
			}
		}
		buffer.WriteString(")")
	}
	return buffer.String()
}

func (search *HCSearch) combineAndCalcAIC() {
	newCluster := new(Cluster)
	for k, c := range search.Clusters {
		for _, site := range c.Sites {
			newCluster.Sites = append(newCluster.Sites, site)
		}
		delete(search.Clusters, k)
	}
	BMOptimBLEM(search.Tree, 10)
	ll := CalcUnrootedLogLike(search.Tree.Rt, true)
	newCluster.LogLike = ll
	fmt.Println("Single cluster AIC:", search.calcAIC())
}

func (search *HCSearch) singleStartingCluster() {
	clus := make(map[int]*Cluster)
	siteClust := make(map[int]int)
	var clustLabs []int
	cur := new(Cluster)
	clus[0] = cur
	clustLabs = append(clustLabs, 0)
	cur.SiteWeights = make(map[int]float64)
	for k := range search.Tree.Rt.ContData {
		//cur := clus[0]
		cur.Sites = append(cur.Sites, k)
		siteClust[k] = 0
	}
	if len(cur.Sites) != 0 {
		ClusterMissingTraitsEM(search.Tree, cur, 10)
		clustll := 0.0
		for _, site := range cur.Sites {
			curll := SingleSiteLL(search.Tree.Rt, site)
			clustll += curll
		}
		cur.LogLike = clustll
	}
	search.Clusters = clus
	search.SiteAssignments = siteClust
	//search.calcClusterSiteWeights()
	tipcount := 0
	for _, n := range search.PreorderNodes {
		if len(n.Chs) == 0 {
			tipcount++
		}
	}
	search.NumTraits = math.Log(float64(len(clus))) * float64(tipcount)
	search.CurrentAIC = search.calcAIC()

}

func (search *HCSearch) randomStartingClusters() {
	clus := make(map[int]*Cluster)
	siteClust := make(map[int]int)
	var clustLabs []int
	for i := 0; i < 20; i++ { //create clusters
		cur := new(Cluster)
		clus[i] = cur
		clustLabs = append(clustLabs, i)
		cur.SiteWeights = make(map[int]float64)
	}
	for k := range search.Tree.Rt.ContData {
		lab := rand.Intn(20)
		cur := clus[lab]
		cur.Sites = append(cur.Sites, k)
		siteClust[k] = lab
	}
	for k, cur := range clus {
		if len(cur.Sites) != 0 {
			ClusterMissingTraitsEM(search.Tree, cur, 10)
			clustll := 0.0
			for _, site := range cur.Sites {
				curll := SingleSiteLL(search.Tree.Rt, site)
				clustll += curll
			}
			cur.LogLike = clustll
		} else {
			delete(clus, k) //delete empty clusters
		}
	}

	search.Clusters = clus
	search.SiteAssignments = siteClust
	//search.calcClusterSiteWeights()
	tipcount := 0
	for _, n := range search.PreorderNodes {
		if len(n.Chs) == 0 {
			tipcount++
		}
	}
	search.NumTraits = math.Log(float64(len(clus))) * float64(tipcount)
	search.CurrentAIC = search.calcAIC()
}
