package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	// "sync"
)

var (
	min_range        int
	max_range        int
	filter_cds_range bool
)

type Sequences map[string]string
type Frequencies map[string]uint16
type CDS_Ranges map[string][2]int

type gene_dupe_info struct {
	gene         string
	read_length  int
	dupe_indices []int
}

type dupe_file_info struct {
	read_length int
	data        Frequencies
}

func init() {
	flag.IntVar(&min_range, "min", 26, "The minimum read length (inclusive)")
	flag.IntVar(&max_range, "max", 30, "The maximum read length (inclusive)")
	flag.BoolVar(&filter_cds_range, "filter_cds_range", false, "Whether to filter out indices outside the CDS")
}

func main() {
	flag.Parse()

	// read in the sequence json file
	jsonFile, err := os.Open("data/sequence_dict.json")
	panic_on_err(err)
	defer jsonFile.Close()
	byteValue, _ := ioutil.ReadAll(jsonFile)
	var sequences Sequences
	json.Unmarshal(byteValue, &sequences)

	jsonFile2, err2 := os.Open("data/cds_ranges.json")
	panic_on_err(err2)
	defer jsonFile2.Close()
	byteValue2, _ := ioutil.ReadAll(jsonFile2)
	var cds_ranges CDS_Ranges
	json.Unmarshal(byteValue2, &cds_ranges)

	fmt.Println("Reading in dupe_seq JSON files...")
	dupe_maps := make(map[int]Frequencies)

	file_read_channel := make(chan dupe_file_info)

	for i := min_range; i <= max_range; i++ {
		go func(read_length int) {
			filename := fmt.Sprintf("data/dupe_seq/dupe_seq_%d.json", read_length)

			// check if file exists
			if _, err := os.Stat(filename); err == nil {
				jsonFile, err := os.Open(filename)
				panic_on_err(err)
				byteValue, _ := ioutil.ReadAll(jsonFile)
				dupe_map := make(Frequencies)
				json.Unmarshal(byteValue, &dupe_map)
				jsonFile.Close()
				file_read_channel <- dupe_file_info{read_length, dupe_map}
			}
		}(i)
	}

	for i := min_range; i <= max_range; i++ {
		file_read_info := <-file_read_channel
		dupe_maps[file_read_info.read_length] = file_read_info.data
	}

	fmt.Println("Searching gene sequences for dupes...")
	dupe_idx_maps := make(map[int]map[string][]int)
	c := make(chan gene_dupe_info)

	for i := min_range; i <= max_range; i++ {
		dupe_idx_maps[i] = make(map[string][]int)
		for gene, sequence := range sequences {
			go func(gene string, sequence string, read_length int) {
				c <- gene_dupe_info{gene, read_length,
					get_dupe_idxs(sequence, read_length, dupe_maps[read_length], cds_ranges[gene])}
			}(gene, sequence, i)
		}
	}

	for i := min_range; i <= max_range; i++ {
		for j := 0; j < len(sequences); j++ {
			dupe_info := <-c
			dupe_idx_maps[dupe_info.read_length][dupe_info.gene] = dupe_info.dupe_indices
		}
	}

	// write to a json file
	jsonData, err := json.Marshal(dupe_idx_maps)
	panic_on_err(err)

	filename := "data/dupe_idx_maps.json"
	file, err := os.Create(filename)
	panic_on_err(err)
	defer file.Close()

	_, err = file.Write(jsonData)
	panic_on_err(err)

	fmt.Printf("Saved %s\n", filename)

}

func get_dupe_idxs(sequence string, read_length int, dupe_seq Frequencies, cds_range [2]int) (indices []int) {
	indices = []int{}
	begin := 0
	end := len(sequence) - read_length + 1
	if filter_cds_range {
		begin = cds_range[0]
		end = cds_range[1] - read_length + 1
	}

	for i := begin; i < end; i++ {
		_, exists := dupe_seq[sequence[i:i+read_length]]
		if exists {
			indices = append(indices, i)
		}
	}
	return
}

func panic_on_err(err error) {
	if err != nil {
		panic(err)
	}
}
