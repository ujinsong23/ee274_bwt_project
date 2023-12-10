from typing import Tuple
import argparse
import numpy as np
import bitarray
from scl.core.data_encoder_decoder import DataDecoder, DataEncoder
from scl.core.data_block import DataBlock
from scl.core.data_stream import TextFileDataStream
from scl.compressors.huffman_coder import HuffmanEncoder,HuffmanDecoder
from scl.core.encoded_stream import EncodedBlockReader, EncodedBlockWriter
from scl.core.prob_dist import ProbabilityDist
from scl.utils.bitarray_utils import BitArray, uint_to_bitarray, bitarray_to_uint


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--decompress", help="Decompression mode if called.", action='store_true')
parser.add_argument("-i", "--input", required=True, type=str, help="Input file directory")
parser.add_argument("-o", "--output", type=str,
                    help="Output file directory, locate same as input if undefined")
parser.add_argument("--delimiter", default=130, type=int,
                    help="Decimal ascii code for delimiter. It should never occurs in the input sequence. 130 if undefined")


DELIMITER = chr(130)

def encode_prob_dist(prob_dist: ProbabilityDist) -> BitArray:
    """Encode a probability distribution as a bit array

    Args:
        prob_dist (ProbabilityDist): probability distribution over 0, 1, 2, ..., 255
            (note that some probabilities might be missing if they are 0).

    Returns:
        BitArray: encoded bit array
    """
    prob_dict = prob_dist.prob_dict # dictionary mapping symbols to probabilities

    #########################
    # ADD CODE HERE
    # Following functions from utils can be useful to implement this
    # bits = BitArray(), bits.frombytes(byte_array), uint_to_bitarray

    n = 12
    #n = np.ceil(np.log2(max(prob_dist.alphabet)))
    encoded_probdist_bitarray = BitArray()
    
    # prob_list_round=np.floor(np.array(list(prob_dict.values()))*N)
    # max_len = int(np.ceil(np.max(np.log2(prob_list_round))))
    max_len = int(np.ceil(np.log2(float('{0:.22f}'.format(np.max(list(prob_dict.values())))[2:]))))
    encoded_probdist_bitarray += uint_to_bitarray(max_len, 20)

    num_nonzero_symbol = len(prob_dict)
    encoded_probdist_bitarray += uint_to_bitarray(num_nonzero_symbol,n)

    for s in prob_dict:
        encoded_probdist_bitarray+=uint_to_bitarray(s,n)
        #encoded_probdist_bitarray+=uint_to_bitarray(int(np.floor(prob_dict[s],)*N),max_len)
        encoded_probdist_bitarray+=uint_to_bitarray(int(('{0:.22f}'.format(prob_dict[s]))[2:]),max_len)
    #########################

    return encoded_probdist_bitarray
 
def decode_prob_dist(bitarray: BitArray) -> Tuple[ProbabilityDist, int]:
    """Decode a probability distribution from a bit array

    Args:
        bitarray (BitArray): bitarray encoding probability dist followed by arbitrary data

    Returns:
        prob_dist (ProbabilityDist): the decoded probability distribution
        num_bits_read (int): the number of bits read from bitarray to decode probability distribution
    """
    #########################
    # ADD CODE HERE
    # Following functions from utils can be useful to implement this
    # bitarray.tobytes() and bitarray_to_uint()

    n = 12
    max_len = bitarray_to_uint(bitarray[:20])
    num_nonzero_symbol = bitarray_to_uint(bitarray[20:20+n])

    num_bits_read = 20+n+num_nonzero_symbol*(max_len+n)
    prob_dict = dict()

    for i in range(num_nonzero_symbol):
        idx = 20+n+i*(n+max_len)
        
        s = bitarray_to_uint(bitarray[idx:idx+n])

        prob = str(bitarray_to_uint(bitarray[idx+n:idx+n+max_len]))
        prob = float("0."+"0"*(22-len(prob))+prob)
        prob_dict[s]=prob

    #########################

    prob_dist = ProbabilityDist(prob_dict)
    return prob_dist, num_bits_read

class FMINDEX:
    def __init__(self,source,is_bwt,keys=[],delimiter=DELIMITER):
        
        self.delimiter = delimiter
        self.is_bwt = is_bwt
        self.OCCtable = {}
        if self.is_bwt:
            assert type(source)==str
            if keys:
                self.keys = keys
            else:
                self.keys = sorted(set(source))
            self.n = len(source)
            
            # initialize OCC and C with source(bwt)
            self.OCC_from_bwt(source)
        else :
            assert isinstance(source, bitarray.bitarray)
            if keys:
                self.keys = keys
            else:
                self.keys = [chr(c) for c in range(128)]
            
            self.Mapping = {}
            self.Mapping2 = {}
            # initialize OCC and C with source(bit array of compressed text)
            self.OCC_from_compressed(source)

    def OCC_from_bwt(self,bwt):
        
        self.OCCtable = {c:np.zeros(self.n+1,dtype=int) for c in self.keys}
        for k in range(self.n):
            self.OCCtable[bwt[k]][k+1:]=self.OCCtable[bwt[k]][k]+1
        self.C = self.C_from_Occ()

    def create_OCC_entry_compressed(self, encoded_bitarray: DataBlock, decoded_length, offset):

        num_bits_consumed = 0
        added_length = 1
        last_decoded_symbol = -1
        prvoffset = offset
        # continue decoding until we reach leaf node

        while num_bits_consumed < len(encoded_bitarray):

            curr_node = self.HuffmanTree.root_node
            prev_num_bits_consumed = num_bits_consumed
            while not curr_node.is_leaf_node:
                bit = encoded_bitarray[num_bits_consumed]
                if bit == 0:
                    curr_node = curr_node.left_child
                else:
                    curr_node = curr_node.right_child
                num_bits_consumed += 1
            # as we reach the leaf node, the decoded symbol is the id of the node
            decoded_symbol = curr_node.id

            if decoded_symbol < 150 and last_decoded_symbol != -1:

                num_bits_consumed = prev_num_bits_consumed ## Don't read the next symbol
                for chars in self.keys:
 
                    if len(self.OCCtable[chars])==0:
                        u = 0
                    else:
                        u = self.OCCtable[chars][-1]

                    if chars != self.MTFTable[0]:
                        self.OCCtable[chars].append(u)
                    else:
                        self.OCCtable[chars].append(u+added_length)
                        
                decoded_length = decoded_length + added_length
                self.Mapping[offset] = decoded_length-1
                self.Mapping2[offset] = self.MTFTable[0]
                offset += 1

                break

            elif decoded_symbol < 150:
                
                c = self.MTFTable[decoded_symbol]
                self.MTFTable.pop(decoded_symbol)
                self.MTFTable.insert(0, c)
                last_decoded_symbol = decoded_symbol

            else:
                if decoded_symbol == 150:
                    added_length = added_length * 2
                else:
                    added_length = added_length * 2 + 1

        if prvoffset == offset:
            for chars in self.keys:
 
                if len(self.OCCtable[chars])==0:
                    u = 0
                else:
                    u = self.OCCtable[chars][-1]

                if chars != self.MTFTable[0]:
                    self.OCCtable[chars].append(u)
                else:
                   self.OCCtable[chars].append(u+added_length)
                        
            decoded_length = decoded_length + added_length
            self.Mapping[offset] = decoded_length-1
            self.Mapping2[offset] = self.MTFTable[0]
            offset += 1

        return num_bits_consumed, added_length, offset



    def OCC_from_compressed(self, compressed: DataBlock): 
        ## compressed : DataBlock read from compressed text

        #from scl.HWs.HW1.hw1_p5 import decode_prob_dist
        prob_dist, num_bits_read_prob_dist_encoder = decode_prob_dist(compressed)
        huffman_decoder = HuffmanDecoder(prob_dist)
        self.HuffmanTree = huffman_decoder.tree
        self.keys.append(self.delimiter)
        self.MTFTable = self.keys[:]
        #self.MTFTable = [chr(c) for c in range(128)]
        # self.MTFTable.append(self.delimiter)

        for c in self.MTFTable:
            self.OCCtable[c] = []
        num_bits_consumed = num_bits_read_prob_dist_encoder
        length = 0
        offset = 0
        while num_bits_consumed < len(compressed):
            
            num_bits, added_length, newoffset = self.create_OCC_entry_compressed(compressed[num_bits_consumed:], length, offset)
            num_bits_consumed += num_bits
            length += added_length
            offset = newoffset

        self.n = length
        self.C = self.C_from_Occ()

    def OCC(self,c,i):

        if self.is_bwt:
            return self.OCCtable[c][i]
        
        elif i == 0:
            return 0

        elif i == 1:
            if self.OCCtable[c][i-1] !=0:
                return 1
            return 0
        
        elif i == -1:
            k = max(list(self.Mapping.keys()))
            return self.OCCtable[c][k]

        else:
            i = i - 1
            s = 0
            k = max(list(self.Mapping.keys()))
            e = k
            while s!=e:
                m = (s+e)//2+1
                if self.Mapping[m] <= i:
                    s = m
                else:
                    e = m - 1
            if e == k:
                return self.OCCtable[c][e]
            elif self.Mapping2[e+1] != c:
                return self.OCCtable[c][e]
            else: 
                return self.OCCtable[c][e] + i - self.Mapping[e]
                
    def C_from_Occ(self):
        alphabet = self.keys
        C = {alphabet[0]:0}
        for i in range(len(alphabet)-1):
            c = alphabet[i]
            c_next = alphabet[i+1]
            C[c_next] = int(C[c]+self.OCC(c,-1))
        return C
    
    def count_occurence(self,query):
       
        lo, hi = 0, int(self.n)
        for q in query[::-1]:

            if lo>=hi or not(q in self.C.keys()):
                print(f'query \'{query}\' cannot be found!')
                return None, None
            
            #print(lo, hi, q, self.C[q], self.OCC(q, lo), self.OCC(q, hi))
      
            lo = self.C[q]+self.OCC(q,lo)
            hi = self.C[q]+self.OCC(q,hi)

        return lo, hi

    def get_loc_list(self):
        i=0
        loc_list = np.zeros(self.n,dtype=int)
        for t in range(self.n,0,-1):
            for c in self.C.keys():
                if self.OCC(c,i)!=self.OCC(c,i+1):
                    break
            i = int(self.C[c]+self.OCC(c,i))
            loc_list[i] = t-2
        return loc_list

    def locate_occurence(self, query):
        loc_list = self.get_loc_list()
        lo, hi = self.count_occurence(query)
        if lo==None:
            return None
        indices = loc_list[lo:hi]
        return indices
    
class BurrowsWheelerTransform:
    def __init__(self, delimiter=DELIMITER):
        # NOTE: we assume delimiter not present in the input
        self.delimiter = delimiter

    def forward(self, data_block: DataBlock):
        """
        Generates the forward transform of BWT
        NOTE: for consistency all forward and inverse functions take in as input
        a DataBlock
        """

        # create a string using data_block
        s = "".join(data_block.data_list)
        s = s + self.delimiter

        ###############################################
        # ADD CODE HERE
        # to generate bwt_str (BWT transformed string)
        # Note: remember to add the delimiter to the string!
        
        n = len(s)
        M = n+1
        init_ = []
        rep = []
        rep2 = []
        lenrep2 = []
        see = []

        for i in range(n):
            init_.append(M*ord(s[i])+i)
            rep.append(0)
            rep2.append([])
            lenrep2.append(0)


        init_.sort()

        k = 0
        m = 0
        for (j,x) in enumerate(init_):
            rep[x%M] = k
            rep2[k].append(x%M)
            lenrep2[k] += 1
            m = m + 1
            if j+1 < n and init_[j+1]//M!=init_[j]//M:
                k = k + m
                m = 0

        for i in range(n):
            if lenrep2[i] > 1:
                see.append(i)

        rmvtme = 0
        stage = 1
        while stage < n:

            rmvsee = []
            apdsee = []
            if len(see) == 0:
                break

            for rank in see:
                tmp = []
                tmprep2 = []

                for k1 in rep2[rank]:
                    if stage + k1 < n:
                        tmp.append(M*rep[stage+k1]+k1)
                    else:
                        tmp.append(M*rep[stage]+k1)

                tmp.sort()
                lentmp = len(tmp)
                k = 0
                m = 0

                for (j,x) in enumerate(tmp):

                    rep[x%M] = k + rank
                    if(k != 0):

                        rep2[k + rank].append(x%M)
                        lenrep2[rank] -= 1
                        lenrep2[k + rank] += 1

                    else:
                        tmprep2.append(x%M)

                    m = m + 1
                    if j+1 < lentmp and tmp[j+1]//M!=tmp[j]//M:
                        if (k != 0 and lenrep2[k + rank] > 1):
                            apdsee.append(rank+k)
                        k = k + m
                        m = 0
                    elif j+1 == lentmp:
                        if (k != 0 and lenrep2[k+rank] > 1):
                            apdsee.append(rank+k)

                rep2[rank] = tmprep2

                if len(rep2[rank]) < 2:
                    rmvsee.append(rank)

            for j in rmvsee:
                see.remove(j)
            see.extend(apdsee)

            stage = stage << 1

        bwt = []
        bwt_str = ""
        for i in s:
            bwt.append(i)
        for (j,i) in enumerate(rep):
            if j==0:
                bwt[i] = self.delimiter
            else:
                bwt[i] = s[j-1]
        for i in bwt:
            bwt_str = bwt_str + i

        ###############################################

        data_bwt_block = DataBlock(list(bwt_str))
        return data_bwt_block

    def inverse(self, bwt_block: DataBlock):
        """
        Generates the inverse of the BWT.
        NOTE: for consistency all forward and inverse functions take in as input
        a DataBlock
        """
        s = "".join(bwt_block.data_list)
        N = len(s)

        ###############################################
        # ADD CODE HERE
        # to generate output_str
        # Note: remember to remove the delimiter from the string!

        inverse_bwt = []
        i = s.index(self.delimiter)

        FM = FMINDEX(s,is_bwt=True)
        C = FM.C

        for _ in range(N):
            c = s[i]
            inverse_bwt.insert(0,c)
            i = int(C[c]+FM.OCC(c,i))

        return DataBlock(inverse_bwt[:-1])
        ###############################################
    
class MoveToFrontTransform:
    def __init__(self):
        # NOTE: this table should work for this HW
        self.table = [chr(c) for c in range(128)]
        self.table.append(DELIMITER)

    def forward(self, data_block):
        table = self.table.copy()
        output_data_list = []
        for c in data_block.data_list:
            rank = table.index(c)  # Find the rank of the character in the dictionary [O(k)]
            output_data_list.append(rank)  # Update the encoded text

            # Update the table
            table.pop(rank)
            table.insert(0, c)
        return DataBlock(output_data_list)

    def inverse(self, data_block_mtf):
        decoded_data_list = []
        table = self.table.copy()

        for rank in data_block_mtf.data_list:
            c = table[rank]
            decoded_data_list.append(c)

            # Update the dictionary
            table.pop(rank)
            table.insert(0, c)

        return DataBlock(decoded_data_list)

class RunLengthEncoding:
    def __init__(self):
        pass

    def forward(self, data_block):
        
        output_data_list = []
        c = 0
        k = 0 ## Store consecutive # of 0s.
        sp0 = 150 ## Special char to note 0
        sp1 = 151 ## Special char to note 1
        for x in data_block.data_list:
            if x != 0:
                if k != 0:
                    rle = uint_to_bitarray(k+1)
                    for b in rle:
                        if b == 0:
                            output_data_list.append(sp0)
                        elif c == 1:
                            output_data_list.append(sp1)
                        else:
                            c = 1
                    k = 0
                    c = 0
                output_data_list.append(x)
            elif x == 0:
                k = k + 1
        if k != 0:
            rle = uint_to_bitarray(k+1)
            for b in rle:
                if b == 0:
                    output_data_list.append(sp0)
                elif c == 1:
                    output_data_list.append(sp1)
                else:
                    c = 1
            
        return DataBlock(output_data_list)
    
    def inverse(self, data_block):

        output_data_list = []
        sp0 = 150
        sp1 = 151 
        k = 1
        for x in data_block.data_list:
            if x != sp0 and x != sp1:
                if k > 1:
                    output_data_list.extend([0]*(k-1))
                    k = 1
                output_data_list.append(x)
            elif x == sp0:
                k = k * 2
            else:
                k = k * 2 + 1
        if k > 1:
            output_data_list.extend([0]*(k-1))
            k = 1
        return DataBlock(output_data_list)

class HuffmanEmpiricalEncoder(DataEncoder):
    def encode_block(self, data_block: DataBlock):

        # get the empirical distribution of the data block
        prob_dist = data_block.get_empirical_distribution()

        # uncomment below to print Huffman tree
        #print_huffman_tree(prob_dist)

        # create Huffman encoder for the empirical distribution
        huffman_encoder = HuffmanEncoder(prob_dist)
        # encode the data with Huffman code
        encoded_data = huffman_encoder.encode_block(data_block)
        # return the Huffman encoding prepended with the encoded probability distribution
        return encode_prob_dist(prob_dist) + encoded_data

    def encode_file(self, input_file_path: str, encoded_file_path: str, block_size: int):
        """utility wrapper around the encode function using Uint8FileDataStream
        Args:
            input_file_path (str): path of the input file
            encoded_file_path (str): path of the encoded binary file
            block_size (int): choose the block size to be used to call the encode function
        """
        # call the encode function and write to the binary file
        with Uint8FileDataStream(input_file_path, "rb") as fds:
            with EncodedBlockWriter(encoded_file_path) as writer:
                self.encode(fds, block_size=block_size, encode_writer=writer)

class HuffmanEmpiricalDecoder(DataDecoder):
    def decode_block(self, encoded_block: DataBlock):
        # first decode the probability distribution
        prob_dist, num_bits_read_prob_dist_encoder = decode_prob_dist(encoded_block)
        # now create Huffman decoder
        huffman_decoder = HuffmanDecoder(prob_dist)
        # now apply Huffman decoding
        decoded_data, num_bits_read_huffman = huffman_decoder.decode_block(
            encoded_block[num_bits_read_prob_dist_encoder:]
        )
        # verify we read all the bits provided
        assert num_bits_read_huffman + num_bits_read_prob_dist_encoder == len(encoded_block)
        return decoded_data, len(encoded_block)

    def decode_file(self, encoded_file_path: str, output_file_path: str):
        """utility wrapper around the decode function using Uint8FileDataStream
        Args:
            encoded_file_path (str): input binary file
            output_file_path (str): output (text) file to which decoded data is written
        """
        # read from a binary file and decode data and write to a binary file
        with EncodedBlockReader(encoded_file_path) as reader:
            with Uint8FileDataStream(output_file_path, "wb") as fds:
                self.decode(reader, fds)

def test_bwt_transform():
    bwt_transform = BurrowsWheelerTransform(delimiter='$')

    sample_inputs = ["BANANA", "abracadabraabracadabraabracadabra", "hakunamatata"]
    expected_bwt_outputs = ["ANNB$AA", "arrrdddaa$rrrcccaaaaaaaaaaaabbbbbb", "athntm$aauaak"]
    for sample_input, expected_bwt_str in zip(sample_inputs, expected_bwt_outputs):
        print("\n" + "-" * 20)
        print(f"Input string: {sample_input}")

        # Get the BWT transformed string
        block = DataBlock(list(sample_input))
        bwt_block = bwt_transform.forward(block)
        bwt_str = "".join(bwt_block.data_list)
        print(f"BWT transfomed string: {bwt_str}")
        assert bwt_str == expected_bwt_str

        # get the inverse BWT
        inv_bwt = bwt_transform.inverse(bwt_block)
        inv_bwt_str = "".join(inv_bwt.data_list)
        print(f"I-BWT: {inv_bwt_str}")
        assert sample_input == inv_bwt_str

def test_mtf_transform():
    mtf = MoveToFrontTransform()

    sample_inputs = [
        "BANANA",
        "BNN~AAA",
        "abracadabraabracadabraabracadabra",
        "rrdd~aadrrrcccraaaaaaaaaaaabbbbbba"
    ]
    for sample_input in sample_inputs:

        # create MTF forward transforms for the given strings
        block = DataBlock(list(sample_input))
        mtf_block = mtf.forward(block)
        print("\n" + "-" * 20)
        print(f"Input str: {sample_input}")
        print(f"MTF: {mtf_block.data_list}")

        inv_mtf = mtf.inverse(mtf_block)
        inv_mtf_str = "".join(inv_mtf.data_list)
        print(f"MTF inverse: {inv_mtf_str}")
        assert inv_mtf_str == sample_input

def test_bwt_mtf_rle():

    bwt = BurrowsWheelerTransform()
    mtf = MoveToFrontTransform()
    rle = RunLengthEncoding()

    sample_inputs = ["BANANA", "abracadabraabracadabraabracadabra", "hakunamatata"]
    for sample_input in sample_inputs:
        print("\n" + "-" * 20)
        print(f"Input string: {sample_input}")

        # Get the BWT transformed string
        block = DataBlock(list(sample_input))
        bwt_block = bwt.forward(block)
        bwt_str = "".join(bwt_block.data_list)
        print(f"BWT transfomed string: {bwt_str}")

        mtf_bwt_block = mtf.forward(bwt_block)
        print(f"MTF + BWT transfomed string: {mtf_bwt_block.data_list}")

        rle_mtf_bwt_block = rle.forward(mtf_bwt_block)
        print(f"RLE + MTF + BWT transfomed string: {rle_mtf_bwt_block.data_list}")
        rleinv_block = rle.inverse(rle_mtf_bwt_block)
        print(f"RLEINV + MTF + BWT transfomed string: {rleinv_block.data_list}")
        assert rleinv_block.data_list == mtf_bwt_block.data_list

def test_bwt_mtf_entropy():
    DATA_BLOCK_SIZE = 100000
    FILE_PATH = "sherlock_ascii.txt"

    # read in DATA_BLOCK_SIZE bytes
    with TextFileDataStream(FILE_PATH, "r") as fds:
        data_block = fds.get_block(block_size=DATA_BLOCK_SIZE)
    print()
    print(f"Input data: 0-order Empirical Entropy: {data_block.get_entropy():.4f}")

    bwt = BurrowsWheelerTransform()
    data_bwt_transformed = bwt.forward(data_block)
    # print(''.join(data_block.data_list))
    # print(''.join(data_bwt_transformed.data_list))
    print(f"Input data + BWT: 0-order Empirical Entropy: {data_bwt_transformed.get_entropy():.4f}")

    mtf = MoveToFrontTransform()
    data_mtf_transformed = mtf.forward(data_block)
    print(f"Input data + MTF: 0-order Empirical Entropy: {data_mtf_transformed.get_entropy():.4f}")

    bwt = BurrowsWheelerTransform()
    mtf = MoveToFrontTransform()
    data_bwt_transformed = bwt.forward(data_block)
    data_bwt_mtf_transformed = mtf.forward(data_bwt_transformed)
    # print(data_bwt_mtf_transformed.data_list)
    print(f"Input data + BWT + MTF: 0-order Empirical Entropy: {data_bwt_mtf_transformed.get_entropy():.4f}")
    print(len(data_bwt_mtf_transformed.data_list))

    rle = RunLengthEncoding()
    data_bwt_mtf_rle_transformed = rle.forward(data_bwt_mtf_transformed)
    print(f"Input data + BWT + MTF + RLE: 0-order Empirical Entropy: {data_bwt_mtf_rle_transformed.get_entropy():.4f}")
    print(len(data_bwt_mtf_rle_transformed.data_list))

def test_BWT_compression():
    
    DATA_BLOCK_SIZE = 100000
    FILE_PATH = "sherlock_ascii.txt"

    # read in DATA_BLOCK_SIZE bytes
    with TextFileDataStream(FILE_PATH, "r") as fds:
        data_block = fds.get_block(block_size=DATA_BLOCK_SIZE)
  
    bwt = BurrowsWheelerTransform()
    mtf = MoveToFrontTransform()
    rle = RunLengthEncoding()

    data_bwt_transformed = bwt.forward(data_block)
    data_bwt_mtf_transformed = mtf.forward(data_bwt_transformed)
    data_bwt_mtf_rle_transformed = rle.forward(data_bwt_mtf_transformed)

    encoder = HuffmanEmpiricalEncoder()   
    encoded_block = encoder.encode_block(data_bwt_mtf_rle_transformed)

    decoder = HuffmanEmpiricalDecoder()
    decoded_block = decoder.decode_block(encoded_block)[0]

    decoded_bwt_mtf = rle.inverse(decoded_block)
    decoded_bwt = mtf.inverse(decoded_bwt_mtf)
    decoded_result = bwt.inverse(decoded_bwt)

    assert decoded_result.data_list == data_block.data_list

if __name__ == "__main__":

    args = parser.parse_args()

    BLOCKSIZE = 10_000_000 
    DELIMITER = chr(args.delimiter)
    bwt = BurrowsWheelerTransform()
    mtf = MoveToFrontTransform()
    rle = RunLengthEncoding()
    
    if args.decompress:

        assert args.input.endswith(".bwtz"), "input format error : .bwtz file is required for BWT decompression"
        inp = args.input

        if args.output:
            assert args.input.endswith(".txt"), "output format error : .txt file is required for BWT decompression"
            out = args.output
        else :
            out = inp[:-5]+'_decompressed.txt'

        decoder = HuffmanEmpiricalDecoder()


        with EncodedBlockReader(inp) as reader:
            with TextFileDataStream(out, "w") as fds:
            
                while True:

                    data_block = reader.get_block()
                    if data_block == None:
                        break

                    decoded_block = bwt.inverse(mtf.inverse(rle.inverse(decoder.decode_block(data_block)[0])))
                    fds.write_block(decoded_block)
                    
    else:

        assert args.input.endswith(".txt"), "input format error : .txt file is required for BWT compression"
        inp = args.input
        
        if args.output:
            assert args.input.endswith(".bwtz"), "output format error : .bwtz file is required for BWT decompression"
            out = args.output
        else :
            out = inp[:-4]+'.bwtz'

        encoder = HuffmanEmpiricalEncoder()

        with EncodedBlockWriter(out) as writer:
            with TextFileDataStream(inp, "r") as fds:
                
                while True:

                    data_block = fds.get_block(block_size = BLOCKSIZE)
                    if data_block == None:
                        break   

                    encoded_block = encoder.encode_block(rle.forward(mtf.forward(bwt.forward(data_block))))
                    writer.write_block(encoded_block)
