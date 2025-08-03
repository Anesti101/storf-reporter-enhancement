# --------------------------------------------------------
# Import: OrderedDict from collections
# Purpose: Maintain the insertion order of ORFs throughout
# the filtering process, especially after sorting by
# length or type. This ensures predictable and stable output
# when iterating or exporting results.
# --------------------------------------------------------
from collections import OrderedDict 
import argparse
from utilss import sortORFs, sortORFs_by_strand  # or from your_module import sortORFs


# --------------------------------------------------------
# Function: storf_type_score
# Purpose: Assign a priority score based on StORF type
# Logic:
#   - If the StORF type is 'Con-StORF', return score 0 (higher priority)
#   - Otherwise, return score 1 (regular StORF)
# Input:
#   storf: tuple containing ORF information
# Returns:
#   int: score (0 or 1)
# --------------------------------------------------------

def storf_type_score(storf):
    return 0 if storf[1] == 'Con-StORF' else 1


# --------------------------------------------------------
# Function: tile_filtering
# Purpose: Filter overlapping ORFs based on length or StORF type priority
# Logic:
#   - Sort ORFs by chosen strategy (length or type priority)
#   - Remove overlapping or nested ORFs using user-specified overlap threshold
#   - Return filtered and re-sorted ORFs
# Input:
#   storfs: OrderedDict of ORFs with positions as keys
#   options: object with filtering and ordering preferences
# Returns:
#   final_filtered_storfs: OrderedDict of filtered ORFs
# --------------------------------------------------------

def tile_filtering(storfs, options):
    # check which sorting strategy to use, default being 'length'
    if hasattr(options, 'priority_strategy'):
        strategy = options.priority_strategy
    else:
        strategy = 'length'
    

    # appply sorting strategy 
    if strategy == 'length':
        # sort by ORF length (index 3), descending
        storfs = sorted(storfs.items(), key=lambda s: s[1][3], reverse=True)
    elif strategy == 'storf_type':
        # sort by storf type score (con-storf first), then lenght 
        storfs = sorted(storfs.items(), key=lambda s: (storf_type_score(s), -s[1][3]))
    else:
         # Unsupported strategy
        raise ValueError(f"Unsupported priority strategy: {strategy}")
    

    # Convert sorted list to list for mutability
    ordered_by_priority = list(storfs)
    num_storfs = len(ordered_by_priority)
    i = 0

    # Strat filtering process
    while i < num_storfs:
        pos_x, data_x = ordered_by_priority[i] # get the first ORF
        start_x = int(pos_x.split(',')[0]) # get start position 
        stop_x = int(pos_x.split(',')[-1]) # get start position 

        j = i + 1 # start checking the next ORF

        # Check for overlap with subsequent ORFs
        while j < num_storfs:
            pos_y, data_y = ordered_by_priority[j] # get the next ORF 
            start_y = int(pos_y.split(',')[0]) # get start position 
            stop_y = int(pos_y.split(',')[-1]) # get start position 

            # No overlap
            if start_y >= stop_x or stop_y <= start_x:
                j += 1 # continue to next ORF

            # Fully nested ORF 
            elif start_y >= start_x and stop_y <= stop_x: 
                ordered_by_priority.pop(j) # remove nested ORf 
                num_storfs -= 1 # update the number of ORFs
            else:
                # Calculate overlap length 
                overlap_start = max(start_x, start_y) # get the start of the overlap
                overlap_end = min(stop_x, stop_y) # get the end of the overlap
                overlap = max(0, overlap_end - overlap_start + 1) # +1 to include both ends 

                # Remove if overlap exceeds threshold
                if overlap >= options.overlap_nt: 
                    ordered_by_priority.pop(j) # remove the overlapping stORF
                    num_storfs -=1 # update the number of ORFs
                else:
                    j += 1 # continue to next ORF

        i += 1 # update the index of the first ORf
        
        # Convert filtered list to OrderedDict 
    filtered_storfs = OrderedDict(ordered_by_priority) # convert to OrderedDict 

        # final sorting based on options 
    if options.storf_order == 'start_pos': # sort by start position
        final_filtered_storfs = sortORFs(filtered_storfs) # sort by start position 
    elif options.storf_order == 'strand': # sort by strand 
        final_filtered_storfs = OrderedDict() # initialise empty OrderedDict 
        storf_nums = sorted([item[-1] for item in filtered_storfs.values()]) # get the sorted list of storf numbers
        for num in storf_nums: # iterate through the sorted list of storf numbers 
            for key, value in filtered_storfs.items(): # iterate through the filtered storfs 
                 if value[-1] == num: # check if the storf number matches 
                    final_filtered_storfs.update({key: value}) # update the final filtered storfs 
                    break
    else:
        final_filtered_storfs = filtered_storfs # keep the order as is 

    return final_filtered_storfs # return the filtered ORfs 

 
# -------------------------------------------------------
# MAIN FUNCTION
# Purpose: Run tile_filtering on example data using argparse
# Logic:
#   - Parse user-supplied options from command line
#   - Create mock ORFs as test input
#   - Run the filter
#   - Print the results
# --------------------------------------------------------

def main():
    # Create the argument parser for command-line usage
    parser = argparse.ArgumentParser(description='Standalone StORF Filtering Tool using tile_filtering().')

    # Required input FASTA-like dictionary in mock form (position: [info])
    parser.add_argument('-priority', dest='priority_strategy', default='length',
                        choices=['length', 'storf_type'],
                        help='Filtering strategy. "length" prioritises longer ORFs. '
                             '"storf_type" prioritises Con-StORFs.')

    parser.add_argument('-olap', dest='overlap_nt', default=50, type=int,
                        help='Default = 50: Max nt overlap allowed between ORFs.')

    parser.add_argument('-so', dest='storf_order', default='start_pos',
                        choices=['start_pos', 'strand'],
                        help='How to sort final ORFs: by start position or strand order.')

    # Parse arguments
    options = parser.parse_args()

    # ------------------------------------------
    # EXAMPLE TEST DATA — You can replace this with file loading
    # Keys = "start,stop" strings; values = [sequence, frame, strand, length, type, id]
    # ------------------------------------------
    test_storfs = OrderedDict({
        "100,300": ["ATG...", 1, '+', 201, 'StORF', 0],
        "150,320": ["ATG...", 1, '+', 171, 'StORF', 1],
        "310,500": ["ATG...", 1, '+', 191, 'Con-StORF', 2],
        "490,700": ["ATG...", 1, '+', 211, 'StORF', 3],
        "650,850": ["ATG...", 1, '+', 201, 'Con-StORF', 4],
        "900,1100": ["ATG...", 1, '+', 201, 'StORF', 5]
    })

    # ------------------------------------------
    # Run the tile filtering function
    # ------------------------------------------
    filtered_storfs = tile_filtering(test_storfs, options)

    # ------------------------------------------
    # Output the filtered results
    # ------------------------------------------
    print("\nFiltered StORFs:\n")
    for pos, data in filtered_storfs.items():
        print(f"{pos}: {data}")

if __name__ == "__main__":
    main()

            
        

    
    

    
     


