package main

import (
	"flag"
	"fmt"
	"os"
	"time"
    "gophy"
)

func main() {
	if len(os.Args) < 2 {
		fmt.Println("seq -fn filename")
		os.Exit(1)
	}
	wks := flag.Int("wks", 1, "an int")
	var fn string
	flag.StringVar(&fn, "fn", "", "filename")

	flag.Parse()

	if len(fn) == 0 {
		fmt.Println("need a filename")
		os.Exit(1)
	}

	sqs := gophy.ReadSeqsFromFile(fn)
	jobs := make(chan []int, len(sqs)*len(sqs))
	results := make(chan float32, len(sqs)*len(sqs))

	start := time.Now()

	for w := 1; w <= *wks; w++ {
		go gophy.PNW(sqs, jobs, results)
	}

	for i, _ := range sqs {
		for j, _ := range sqs {
			if i < j {
				jobs <- []int{i, j}
			}
		}
	}
	close(jobs)

	for i, _ := range sqs {
		for j, _ := range sqs {
			if i < j {
				<-results
			}
		}
	}
	end := time.Now()
	fmt.Println(end.Sub(start))
}
