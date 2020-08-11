package clustering

//MaxClustLab returns the maximum value in a map of ints used like a set
func MaxClustLab(l map[int]*Cluster) (biggest int) {
	biggest = -10000000
	for i := range l {
		if i > biggest {
			biggest = i
		}
	}
	return
}
