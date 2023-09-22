package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	"sync"
)

var (
	min_range int
	max_range int
)

func init() {
	flag.IntVar(&min_range, "min", 26, "The minimum read length (inclusive)")
	flag.IntVar(&max_range, "max", 30, "The maximum read length (inclusive)")
}

type Sequences map[string]string
type Frequencies map[string]uint16

func main() {
	flag.Parse()

	if min_range > max_range {
		panic("min read length cannot be greater than max read length")
	}

	// read in the sequence json file
	jsonFile, err := os.Open("data/sequence_dict.json")
	panic_on_err(err)
	defer jsonFile.Close()
	byteValue, _ := ioutil.ReadAll(jsonFile)
	var sequences Sequences
	json.Unmarshal(byteValue, &sequences)

	fmt.Printf("Generating Frequency Dicts for [%d, %d]...\n", min_range, max_range)
	freq_dicts_map := map[int][]Frequencies{}
	channels := make(map[int]chan Frequencies)

	for i := min_range; i <= max_range; i++ {
		channels[i] = make(chan Frequencies)
	}

	for i := min_range; i <= max_range; i++ {
		for _, sequence := range sequences {
			go func(seq string, read_length int) { channels[read_length] <- get_freq_dict(seq, read_length) }(sequence, i)
		}
	}

	for i := min_range; i <= max_range; i++ {
		for j := 0; j < len(sequences); j++ {
			freq_dicts_map[i] = append(freq_dicts_map[i], <-channels[i])
		}
	}

	fmt.Println("Merging Frequency Dicts...")

	var wg sync.WaitGroup
	wg.Add(max_range - min_range + 1)
	for i := min_range; i <= max_range; i++ {
		go func(read_length int) {
			defer wg.Done()
			merged := merge_freq_dict_list(freq_dicts_map[read_length])
			fmt.Printf("Read Length %d merged with length %d. Filtering...\n", read_length, len(merged))

			filtered := filter_freq_dict(merged)
			fmt.Printf("Read Length %d filtered with length %d. Saving...\n", read_length, len(filtered))

			// write to a json file
			jsonData, err := json.Marshal(filtered)
			panic_on_err(err)

			filename := fmt.Sprintf("data/dupe_seq/dupe_seq_%d.json", read_length)
			file, err := os.Create(filename)
			panic_on_err(err)
			defer file.Close()

			_, err = file.Write(jsonData)
			panic_on_err(err)

			fmt.Printf("Saved %s\n", filename)

		}(i)
	}

	wg.Wait()

}

// get a frequency dict for one gene, given its sequence and a read_length
func get_freq_dict(sequence string, read_length int) (freq_dict Frequencies) {
	freq_dict = make(Frequencies)
	for i := 0; i < len(sequence)-read_length+1; i++ {
		value, exists := freq_dict[sequence[i:i+read_length]]
		if exists {
			freq_dict[sequence[i:i+read_length]] = value + 1
		} else {
			freq_dict[sequence[i:i+read_length]] = 1
		}
	}
	return
}

// merge two freq dicts, d1 and d2. if the same key exists in both maps,
// the resulting map's value will take the sum of the value in d1 and d2
func merge_freq_dicts(d1, d2 Frequencies) (d3 Frequencies) {
	d3 = make(Frequencies)

	for sequence, count := range d1 {
		d3[sequence] = count
	}

	for sequence, count_d2 := range d2 {
		count_d1, exists_in_d1 := d1[sequence]
		if exists_in_d1 {
			d3[sequence] = count_d1 + count_d2
		} else {
			d3[sequence] = count_d2
		}
	}

	return
}

// merge a list of frequency maps, freq_dicts, into one frequency map
// if the same key exists in multiple maps, the value in the result will
// be the sum of the values across the maps
func merge_freq_dict_list(freq_dicts []Frequencies) (merged Frequencies) {

	n := len(freq_dicts)

	if n == 1 {
		return freq_dicts[0]
	}

	// else
	f1 := freq_dicts[0 : n/2]
	f2 := freq_dicts[n/2 : n]

	// fmt.Printf("Dividing %d into %d and %d\n", n, len(f1), len(f2))

	c := make(chan Frequencies)
	defer close(c)
	go func() { c <- merge_freq_dict_list(f1) }()
	go func() { c <- merge_freq_dict_list(f2) }()
	l1 := <-c
	l2 := <-c

	return merge_freq_dicts(l1, l2)
}

// take in a frequency map freq and return a new map of only keys with counts > 0
func filter_freq_dict(freq Frequencies) (filtered Frequencies) {
	filtered = make(Frequencies)
	for sequence, count := range freq {
		if count > 1 {
			filtered[sequence] = count
		}
	}
	return
}

func panic_on_err(err error) {
	if err != nil {
		panic(err)
	}
}
