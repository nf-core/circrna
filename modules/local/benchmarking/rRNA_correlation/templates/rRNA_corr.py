#!/usr/bin/env python3
import numpy as np
from collections import defaultdict

def process_meta_and_paths(input_string):
    # 1. Strip the first and last character
    stripped_input = input_string.replace("[[", "[").replace("]]", "]")
    
    # 2. Split by "[", resulting in a list
    split_data = stripped_input.split("[")
    
    # 3. For every list member, strip all occurrences of "]"
    cleaned_data = [item.replace("]", "").strip() for item in split_data if item.strip()]
    
    # 4. For every list member, split by ","
    parsed_data = [item.split(",") for item in cleaned_data]
    
    parsed_data = [sublist[:-1] if sublist[-1] == '' else sublist for sublist in parsed_data]

    
    return parsed_data


def calculate_favcount_summary_ratio(file_path):
    with open(file_path.strip(), 'r') as file:
        lines = file.readlines()
    assigned = int(lines[1].split()[1])
    unassigned_multimapping = int(lines[9].split()[1])
    others_sum = sum(int(lines[i].split()[1]) for i in list(range(2, 9)) + list(range(10, 15)))
    ratio = (assigned + unassigned_multimapping)
    #others_sum if others_sum != 0 else 1
    return ratio



def getpairs (beds, favcounts):
    pairs = []
    for bed in beds:
        fav = None
        for f in favcounts:
            if f[0] == bed[0]:
                fav = f
                break
        pairs.append([bed[0],bed[4],bed[-1],fav[-1]])
    return pairs

def get_totals(pairs):
    for i in range(len(pairs)):
            bed = pairs[i][2].strip(" ")
            fav = pairs[i][3]
            with open(bed, 'r') as file:
                bed_value = sum(1 for line in file)
            pairs[i][2] = bed_value
            pairs[i][3] = calculate_favcount_summary_ratio(fav)
    return pairs


def sort_tools(value_pairs):
    grouped_data = defaultdict(list)
    for entry in value_pairs:
        grouped_data[entry[1]].append(entry)    
    return dict(grouped_data)


def compute_correlation(value_pairs):
    results = defaultdict(float)    
    for tool in value_pairs.keys():
        bed_list, fav_list = zip(*[sublist[2:] for sublist in value_pairs[tool]])
        correlation = np.corrcoef(bed_list, fav_list)[0,1]
        results[tool] = correlation 
    return dict(results)
        
        
    

# Input data from Nextflow variables
real_bed = "$bed_real_list"
bench_bed = "$bed_bench_list"
real_rRNA = "$rRNA_real_list"
bench_rRNA = "$rRNA_bench_list"

real_bed = process_meta_and_paths(real_bed)
bench_bed = process_meta_and_paths(bench_bed)
real_rRNA = process_meta_and_paths(real_rRNA)
bench_rRNA = process_meta_and_paths(bench_rRNA)

real_pairs = getpairs(real_bed,real_rRNA)
bench_pairs = getpairs(bench_bed,bench_rRNA)


real_value_pairs = get_totals(real_pairs)
bench_value_pairs = get_totals(bench_pairs)

real_value_pairs = sort_tools(real_value_pairs)
bench_value_pairs = sort_tools(bench_value_pairs)

real_corr = compute_correlation(real_value_pairs)
bench_corr = compute_correlation(bench_value_pairs)

with open("real_corr.txt", "w") as file:
    for key, value in real_corr.items():
        # Remove leading space in the key if required
        line = f"{key.strip()}:\\t{value}\\n"
        file.write(line)
        
with open("bench_corr.txt", "w") as file:
    for key, value in bench_corr.items():
        # Remove leading space in the key if required
        line = f"{key.strip()}:\\t{value}\\n"
        file.write(line)

#TODO: both relative or both total. both?
#TODO: both correllations the same can't be!!!

