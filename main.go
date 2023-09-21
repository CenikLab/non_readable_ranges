package main

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"os"
)

type Sequences map[string]string
type Frequencies map[string]int

func main() {
	jsonFile, err := os.Open("human_sequence_dict.json")
	if err != nil {
		panic(err)
	}
	defer jsonFile.Close()

	byteValue, _ := ioutil.ReadAll(jsonFile)
	var sequences Sequences
	json.Unmarshal(byteValue, &sequences)

	fmt.Println("Generating Frequency Dicts...")

	freq_dicts_map := map[int][]Frequencies{}
	channels := map[int]chan Frequency

	// min_range := 21
	// max_range := 21 //40 // this is inclusive

	// for i := min_range; i <= max_range; i++ {
	// 	freq_dicts_map[i] := {}
	// }

	// freq_dicts := []Frequencies{}

	c := make(chan Frequencies)
	for _, sequence := range sequences {
		go func(seq string) { c <- get_freq_dict(seq, 40) }(sequence)
	}

	for i := 0; i < len(sequences); i++ {
		freq := <-c
		freq_dicts_map[21] = append(freq_dicts_map[21], freq)
	}

	fmt.Println("Merging Frequency Dicts...")
	freq_dicts_to_merge := freq_dicts_map[21]
	merged := merge_freq_dict_list(freq_dicts_to_merge)

	fmt.Println(len(merged))

	// filter to > 1

	// save to json file

}

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
	go func() { c <- merge_freq_dict_list(f1) }()
	go func() { c <- merge_freq_dict_list(f2) }()
	l1 := <-c
	l2 := <-c

	return merge_freq_dicts(l1, l2)
}
