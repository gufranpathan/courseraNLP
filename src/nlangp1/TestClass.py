'''
Created on Apr 20, 2014

@author: Gufran.Pathan
'''
import unittest
from tagger import *

class Test(unittest.TestCase):


    def testName(self):
        #Print the new word counts with "_RARE_" replacement for comparision with gene.counts.replaced
        for elem in wordtag_counts:
            print str(wordtag_counts[elem]) + " " + " WORDTAG " + elem[1] + " " + elem[0]


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()