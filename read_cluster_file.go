package gophy

import (
	"fmt"
	"os"
	"strconv"
	"strings"
)

func ReadMCLoutput(clfl string) (clusters map[int]*Cluster) {
	lines := ReadLine(clfl)
	clusters = map[int]*Cluster{}
	for lab, line := range lines {
		if line == "" {
			continue
		}
		ss := strings.Split(line, "\t")
		cur := new(Cluster)
		for _, site := range ss {
			con, err := strconv.Atoi(site)
			if err != nil {
				fmt.Println(site)
				fmt.Println("Couldn't convert a label in line ", lab, " of the input file")
				os.Exit(1)
			}
			cur.Sites = append(cur.Sites, con)
		}
		clusters[lab] = cur
	}
	return
}
