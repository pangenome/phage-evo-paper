#!/usr/bin/env python3

import sys
import os
import re

def main(separator=' !! '):
    try:
        for line in sys.stdin:
            line = line.rstrip()
            if line.startswith('>'):
                header_info = line.split(separator)
                idx = re.findall(r'(?<=##\s).+', header_info[0])[0]
                product = 'Unknown'
                if len(header_info) == 7:
                    product = header_info[3]
                line = f'>{idx}: {product}'
            print(line)
    except BrokenPipeError:
        pass

if __name__ == "__main__":
    main(*sys.argv[1:])
