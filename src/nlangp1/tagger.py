'''
Created on Apr 16, 2014

@author: Gufran.Pathan
'''
import sys
import operator
from nlangp1.hmm import viterbi

class Tagger:
    def __init__(self):
        self.words_count = {}  # Store counts of words; if word<5, will be replaced by the word "_RARE_"
        self.wordtag_counts = {}  # Stores the count of tags for each word
        self.igene_emission_params = {}  # Stores the emission parameter e(word|I-GENE)
        self.tagO_emission_params = {}  # Stores the emission parameter e(word|O)
        self.emission_probs = {"I-GENE":self.igene_emission_params, "O":self.tagO_emission_params}
        self.ngram_counts = {1:{}, 2:{}, 3:{}}  # format==>{1grams:{(u1,v1,w1): count,(u2,v2,w2):count},2grams:{},3grams:{}} 
                                                # where u,v,w are sequence of tags
        self.trigram_probs = {}  # Trigram probablities of trigram (u,v,w) stored here. Calculated by est_trigram_prob fn
        

    '''
                    <------------------------------------------------------------------------------------>
                                        Read words from gene.dev, club by sentence and yield list
                    <------------------------------------------------------------------------------------>
    '''

    def read_sentences(self, test_file_loc):
        "Lazily read sentences from a handle."
        test_file = open(test_file_loc, 'r')
        sentence = []
        for l in test_file:
            if l.strip():
                sentence.append(l.strip())
            else:
                yield sentence
                sentence = []


    '''
                    <------------------------------------------------------------------------------------>
                                                Output gene.p2.out
                    <------------------------------------------------------------------------------------>
    '''
    def printer(self, test_file_loc):
        for input in self.viterbi(test_file_loc):
            tags = input[0]
            sentence = input[1]
            for elem in zip(tags, sentence):  # Take i'th element of both lists and print
                print "%s %s" % (elem[1], elem[0])
            print
            
    '''
                    <------------------------------------------------------------------------------------>
                                        Tag the test data using the Viterbi Algorithm
                    <------------------------------------------------------------------------------------>
    '''
    def viterbi(self, test_file_loc):
        f = 1
        k_tags = [i[0] for i in self.ngram_counts.get(1).keys()]
        start_tags = {'*'}
        stop_tags = {'STOP'}
        
        # Retrieve Backpointers
        def ret_bp(n):
            
            # Add y^n,y^n-1 and 'STOP' to tag list
            y_n, y_np = max(n_pivalues.iteritems(), key=operator.itemgetter(1))[0]
            tags.extend([y_n, y_np, 'STOP'])
            
            # Implement: For k = (n - 2)...1, yk = bp(k + 2; yk+1; yk+2)
            for k in reversed(range(1, n - 1)):
                tags.insert(0, back_pointer.get((k + 2, tags[0], tags[1])))
            
        
        # Pi function: pi(k; u; v) = max for all w (pi(k-1;w;u) * q(v|w,u) * e(x|v))
        def pi(k, u, v):
            if k == 0:
                return 1
            else:
                if k == 1 or k == 2:
                    w_tag = start_tags
                else:
                    w_tag = k_tags
            max_these = {}
            for w in w_tag:
                a = pi_values.get((k - 1, w, u)) * self.trigram_probs.get((w, u, v), 0.0) * self.get_emission_param(elem[k - 1], v)
                max_these[w] = a
            back_pointer[k, u, v] = max(max_these.iteritems(), key=operator.itemgetter(1))[0]
            return max(max_these.values())

        # Call the pi function for each sentence and return max-arg
        for elem in self.read_sentences(test_file_loc):
            n = len(elem)
            pi_values = {}  # format==> {(k,u,v):'value'}
            pi_values[(0, '*', '*')] = 1
            n_pivalues = {}  # format==> {(u,v):'value'}
            back_pointer = {}  # format==>    {(k,u,v):'w'}
            tags = []
            for k in range(1, n + 1):
                if k == 1:
                    u_tag = start_tags
                    v_tag = k_tags
                else:
                    u_tag = v_tag = k_tags
                for u in u_tag:
                    for v in v_tag:
                        pi_values[(k, u, v)] = pi(k, u, v)
                        if k == n:
                            n_pivalues[(u, v)] = pi_values.get((k, u, v)) * self.trigram_probs.get((u, v, 'STOP'), 0.0)
                            
            ret_bp(n)  # Call retrieve backpointers function
            yield [tags, elem]
                    
        
    '''
                    <------------------------------------------------------------------------------------>
                                            Tag the test data using Only emission parameters e(word|tag)
                    <------------------------------------------------------------------------------------>
    '''
    
    def tag_test_data(self, test_file_loc):
        test_file = open(test_file_loc, 'r')
        for line in test_file:
            line = line.strip()
            if line:
                if self.igene_emission_params.has_key(line) or self.tagO_emission_params.has_key(line):
                    igene_prob = self.igene_emission_params.get(line, 0.0)
                    tag0_prob = self.tagO_emission_params.get(line, 0.0)
                else:
                    igene_prob = self.igene_emission_params['_RARE_']
                    tag0_prob = self.tagO_emission_params['_RARE_']
                if igene_prob > tag0_prob:
                    print line + " I-GENE"
                else:
                    print line + " O"
            else:
                print ""   
    '''
                        <------------------------------------------------------------------------------------>
                                            Count total number of "I-GENE"-Tags and "O"-Tags
                        <------------------------------------------------------------------------------------>
    '''
    
    def count_tag(self):
        igene_count = 0.0
        tag0_count = 0.0
        for elem in self.wordtag_counts:
            if elem[1] == "I-GENE":
                igene_count += self.wordtag_counts[elem]
            else:
                tag0_count += self.wordtag_counts[elem]
        return [igene_count, tag0_count]
    
    
    '''
                    <------------------------------------------------------------------------------------>
                                    Estimate Trigram Parameter q(w|u,v) = count(u,v,w)/count(v,w)
                                                    Takes li=[u,v,w] as an argument
                    <------------------------------------------------------------------------------------>
    '''
    def est_trigram_prob(self):
        for elem in self.ngram_counts.get(3):
            prob = self.ngram_counts.get(3).get(elem) / self.ngram_counts.get(2).get(tuple(elem[:2]))
            self.trigram_probs[elem] = prob
    
    '''
                <------------------------------------------------------------------------------------------------------->
                            Calculate emission parameter [e(x|y)=c(x,y)/c(y)] x-->word | y-->I-GENE tag/'O' tag
                <--------------------------------------------------------------------------------------------------------->
    '''
    def get_emission_param(self, word, tag):
        if self.words_count.get(word, 0) < 5:
            return self.emission_probs.get(tag).get('_RARE_')
        else:
            return self.emission_probs.get(tag).get(word, 0.0)
            
    def emission_param(self):
        # Get count of I-GENE tags [c(y)]
        tag_counts = self.count_tag()
  
        for elem in self.wordtag_counts:
            if elem[1] == "I-GENE":
                self.igene_emission_params[elem[0]] = float(self.wordtag_counts[elem]) / tag_counts[0]
            else:
                self.tagO_emission_params[elem[0]] = float(self.wordtag_counts[elem]) / tag_counts[1]
        
    '''
                <--------------------------------------------------------------------------------------------------------->
                                          Replace rare words (count < 5) with "_RARE_"
                <--------------------------------------------------------------------------------------------------------->
    '''
    def replace_rare_words(self):
        # replace rare words (count<5) with "_RARE_"
        self.wordtag_counts[("_RARE_", "I-GENE")] = 0
        self.wordtag_counts[("_RARE_", "O")] = 0
        for elem in self.words_count:
            if self.words_count[elem] < 5:
                if self.wordtag_counts.has_key((elem, "I-GENE")):  # if the word has been tagged as I-GENE, 
                                                                    # add its count to rare and remove it from dict
                    self.wordtag_counts[("_RARE_", "I-GENE")] += self.wordtag_counts[(elem, "I-GENE")]
                    del self.wordtag_counts[(elem, "I-GENE")]  
                                      
                if self.wordtag_counts.has_key((elem, "O")):  # if the word has been tagged as '0', 
                                                                    # add its count to rare and remove it from dict
                    self.wordtag_counts[("_RARE_", "O")] += self.wordtag_counts[(elem, "O")]
                    del self.wordtag_counts[(elem, "O")]    
    
    
    '''
                <--------------------------------------------------------------------------------------------------------->
                                            (1) Take gene.counts file as input and 
                                            (2) Store data in either N-Gram WordTag lists 
                <--------------------------------------------------------------------------------------------------------->
    '''
    
    def extract_TagCounts(self, count_file_loc):
        
        count_file = open(count_file_loc, "r")
        
        for line in count_file:
            
            # Append "WORDTAG" counts to Wordtag_counts list
            line_list = line.split()
            
            if line_list[1] == "WORDTAG":
                self.wordtag_counts[(line_list[3], line_list[2])] = int(line_list[0])
                self.words_count.setdefault(line_list[3], 0)
                self.words_count[line_list[3]] += int(line_list[0])
             
            # add NGram elements  to N-Gram List
            elif line_list[1] == "1-GRAM":
                self.ngram_counts.get(1)[tuple(line_list[2:])] = int(line_list[0])
            elif line_list[1] == "2-GRAM":
                self.ngram_counts.get(2)[tuple(line_list[2:])] = float(line_list[0])
            elif line_list[1] == "3-GRAM":
                self.ngram_counts.get(3)[tuple(line_list[2:])] = int(line_list[0])
        
     
    '''
                        <------------------------------------------------------------------------------------>
                                                        Start main method
                        <------------------------------------------------------------------------------------>
    '''

def main():
    tag = Tagger()
    
    tag.extract_TagCounts(count_file_loc)
    tag.replace_rare_words()
    tag.emission_param()
    # tag.tag_test_data(test_file_loc)            #Use this to tag data with the simple Unigram tagger
    tag.est_trigram_prob()
    tag.printer(test_file_loc)
    # print tag.get_emission_param('amino', 'O')
    

    
if __name__ == '__main__':
    
    if len(sys.argv) != 3:  # Expect exactly three arguments: the "1"-string-argument, gene.counts, gene.dev data file
        count_file_loc = "C:\\Users\\Gufran.Pathan\\NewWorkSpace\\Stanford NLP\\src\\nlangp1\\Resources\\gene.counts" 
        test_file_loc = "D:\\Coursera\\NLangP\\gene.dev"
    else:
        count_file_loc = sys.argv[1]
        test_file_loc = sys.argv[2]
    
    main()    
