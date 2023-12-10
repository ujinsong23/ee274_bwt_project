from bwt import *
import argparse
from scl.core.encoded_stream import EncodedBlockReader, EncodedBlockWriter


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--decompress", help="Decompression mode if called.", action='store_true')
parser.add_argument("-i", "--input", required=True, type=str, help="Input file directory")
parser.add_argument("-o", "--output", type=str,
                    help="Output file directory, locate same as input if undefined")
parser.add_argument("--delimiter", default=130, type=int,
                    help="Decimal ascii code for delimiter. It should never occurs in the input sequence. 130 if undefined")



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
