from bwt import *
import numpy as np
import argparse
import time
from sys import getsizeof
from scl.core.encoded_stream import EncodedBlockReader


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--decompressed", help="Search from decompresssed source", action='store_true')
parser.add_argument("-i", "--input", required=True, type=str,
                    help="Input file directory. \nShould be in the appropriate form,\
                          depending on whether the searching is performed on compressed or decompressed source.")
parser.add_argument("-q", "--query", type=str,
                    help="Query sequence")
parser.add_argument("-c", "--doublecheck", action='store_true',
                    help="Check if the result matches the naive implementation using the built-in Python function ")


def occurence_from_sequence_naive(T, query):
    indices = []
    string = T[:]

    while len(string)>0:
        idx = string.find(query)
        if idx<=0:
            break
        indices.append(idx)
        string = string[idx+1:]

    return len(indices)


if __name__ == "__main__":

    DELIMITER = chr(130)
    args = parser.parse_args()
    query = args.query

    if args.decompressed :

        assert args.input.endswith(".txt"), \
            "input format error : .txt file is required for searching over decompressed data"
        
        inp = args.input
        print(f"Counting occurence of \'{query}\' in '{inp}' using FM index from BWT ...")

        f = open(inp,"r")
        T = f.read()
        block = DataBlock(list(T))

        bwt_transform = BurrowsWheelerTransform()

        bwt_block = bwt_transform.forward(block)
        bwt = "".join(bwt_block.data_list)
        keys = sorted(set(bwt))

        start_pre = time.time()
        FM_index_from_BWT = FMINDEX(bwt, is_bwt=True, keys=keys)
        loc_list = FM_index_from_BWT.get_loc_list()
        end_pre = time.time()

        start = time.time()
        lo, hi = FM_index_from_BWT.count_occurence(query)
        end = time.time()

        n = 0 if lo==None else (hi-lo)

        if args.doublecheck:
            n_naive = occurence_from_sequence_naive(T, query)
            assert n==n_naive, \
                "Detected occurences does not match with result from built-in .find() function"

        print(f"\t> Total {n} occurences detected")
        print(f'\t> Elapsed time : {end-start:.5f}sec (+time for computing FM index tables : {end_pre-start_pre:.5f}sec)')
        print(f'\t> Memory usage : {getsizeof(FM_index_from_BWT.OCCtable[DELIMITER])}')
        
    else :
        assert args.input.endswith(".bwtz"),\
             "input format error : .bwtz file is required for searching over compressed data"
        inp = args.input
        print(f"Counting occurence of '{query}' in '{inp}' using FM index, directly from Compressed ...")

        with EncodedBlockReader(inp) as reader:
            data_block = reader.get_block()
            FM_index_from_compressed = FMINDEX(data_block, is_bwt=False, delimiter=DELIMITER)

        start = time.time()
        lo, hi = FM_index_from_compressed.count_occurence(query)
        end = time.time()

        n = 0 if lo==None else (hi-lo)

        if args.doublecheck:

            bwt = BurrowsWheelerTransform()
            mtf = MoveToFrontTransform()
            rle = RunLengthEncoding()
            decoder = HuffmanEmpiricalDecoder()

            decoded_block = bwt.inverse(mtf.inverse(rle.inverse(decoder.decode_block(data_block)[0])))
            T = "".join(decoded_block.data_list)

            n_naive = occurence_from_sequence_naive(T, query)
            assert n==n_naive, \
                "Detected occurences does not match with result from built-in .find() function"
        
        print(f"\t> Total {n} occurences detected")
        print(f'\t> Elapsed time : {end-start:.5f}sec')
        print(f'\t> Memory usage : {getsizeof(FM_index_from_compressed.OCCtable[DELIMITER])}')