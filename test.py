#!/usr/bin/env python
#-*- coding: utf-8 -*-
import sys
from parserQE.parser import loadFile

def main(name):
    print(f'It is worked with: {name}')
    loadFile(name)







if __name__ == '__main__':
    if(len(sys.argv) == 2):
        main(sys.argv[1])
