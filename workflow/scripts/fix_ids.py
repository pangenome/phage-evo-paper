#!/usr/bin/env python3
import os
import sys
import re

def main(file: str, codes: str) -> None:
    read_table = lambda f: iter(map(lambda x: x.strip().split('\t'), f.readlines()))
    with open(codes, 'r') as f:
        codes = dict(map(lambda x: reversed(x[1:]), read_table(f)))
    with open(file, 'r') as f:
        for line in iter(map(str.strip, f.readlines())):
            to_replace = re.findall(r'^S\d+#|(?<=\t)S\d+#', line)
            c = len(to_replace)
            assert c in (0,2)
            if c == 2:
                for i in set(to_replace):
                    line = re.sub(i, codes[i.replace('#','')] + '#' , line)
            print(line)
    
            
                
                
            
        
        
        

if __name__ == '__main__':
    main(*sys.argv[1:])
