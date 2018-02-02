import numpy as np
import csv
from copy import copy
import re
import sys

def complementary(seq):
    test_seq = seq.replace('G','c')
    test_seq = test_seq.replace('C','G')
    test_seq = test_seq.replace('c','C')
    test_seq = test_seq.replace('A','t')
    test_seq = test_seq.replace('T','A')
    test_seq = test_seq.replace('t','T')
    return test_seq

class Plasmid:
    def __init__(self, seq):
        self.sequence = seq
        self.comp = complementary(seq)
        self.blank_range = [0] * len(seq)
        self.oligo_list = {}
        self.complete_reads = []
        self.oligo_locs = {}
        self.num_iters = 0
        self.partial_reads = []
        
    def range_int(self, st, en):
    	if st < 1:
    		st = input("Start too low, pick new start. ")
    		self.range_int(st, en)
    	elif st > en:
    		st = input("Start larger than end, pick new start. ")
    		self.range_int(st, en)
    	elif en > len(self.sequence):
    		en = input("End out of range, pick new end. ")
    		self.range_int(st, en)

        self.target_seq = copy(self.sequence[st-1:en-1])
        self.target_range = [st-1, en] 
        
    def copy_pattern(self, oligos):
        #find each primers range
        #odd= start, even = end
        count = 1
        for i, o in enumerate(oligos):
            if self.blank_range[o[4][0]] == 0:
                self.blank_range[o[4][0]] = count
                self.oligo_list[count] = [o[0], o[3]]
                self.oligo_locs[count] = [o[4][0], o[4][1]]
                
                count += 1
                self.blank_range[o[4][1]] = count
                count += 1
            else:
                new_odd, new_even = self.adjuster(o[4][0], o[4][1])
                
                self.blank_range[new_odd] = count
                self.oligo_list[count] = [o[0], o[3]]
                self.oligo_locs[count] = [new_odd, new_even]
                
                count += 1
                self.blank_range[new_even] = count
                count += 1

    def adjuster(self, st, en):
        #if one of the original locations is occupied by another primer's range
        mag = 1
        vec = 1
        
        while True:
            new_st = st + (mag*vec)
            new_en = en + (mag*vec)
            
            if self.blank_range[new_st] == 0 and self.blank_range[new_en] == 0:
                return new_st, new_en
            elif vec > 0:
                vec *= -1
            else:
                vec *= -1
                mag += 1

    def seq_finder(self):
        self.find_full_seq(self.target_range[0], 1, [])
        
        if len(self.complete_reads)>=1:
            self.print_reads()
            self.replacer()
        else:
        	self.no_path()

    def replacer(self):
    	rep = input('Replace an oligo (input oligo number, 0 to exit)? ')

    	if rep != 0:
    		try:
    			plasmid.rem_oligo(rep)
    			plasmid.num_iters = 0
    			plasmid.complete_reads = []
    			plasmid.seq_finder()
    		except:
    			print "Oligo not found"
    			self.replacer()
    	else:
    		return

    def find_last_odd(self, curr, lowlim):
        all_rem_odds = [x for x in self.blank_range[lowlim+1:curr] if x%2==1]
        
        if all_rem_odds == []:
            return ['','']
        
        last_odd = all_rem_odds[-1]
        new_curr = plasmid.blank_range.index(last_odd)
        return [last_odd, new_curr]
              
    
    def find_full_seq(self, curr_ind, low_ind, path):
        
        self.num_iters += 1
        
        [next_odd, next_low] = self.find_last_odd(curr_ind, low_ind)

        if next_odd != '':

            path.append(next_odd)
            curr_ind = self.blank_range.index(next_odd + 1)
            
            if curr_ind >= self.target_range[1]:
                self.complete_reads.append(list(path))
                return 
            
            else:

                self.find_full_seq(curr_ind, next_low, path)
                
        else:
            self.partial_reads.append(list(path))
            return 
    
    def rem_oligo(self, oli):
        self.blank_range[self.blank_range.index(oli)] = 0
        self.blank_range[self.blank_range.index(oli+1)] = 0
        
    def print_reads(self):
        for read in self.complete_reads:
            print '--------------'
            for item in read:
                print str(item) + ": " + str(self.oligo_list[item])
            print '--------------'

    def find_empty_ranges(self):
        last_end = self.blank_range.index(self.partial_reads[-1][-1] + 1)
        next_start = self.find_next_odd(last_end)

        if next_start >= self.target_range[1] or next_start <= self.target_range[0]:
        	return

        self.find_full_seq(next_start+1, last_end, [])

        if len(self.complete_reads) == 1:
        	self.partial_reads.append(list(self.complete_reads[0]))
        	return
        else:
        	self.find_empty_ranges()
        
    def find_next_odd(self, curr):
    	future_odds = [x for x in self.blank_range[curr:] if x%2==1]
    	return self.blank_range.index(future_odds[0])

    def no_path(self):

    	# TODO			
        print 'No complete path found'
        self.find_empty_ranges()
        
        print self.partial_reads
        return 0
        #add functionality here
        #should describe missing range and exit function
            
class Oligos:
    def __init__(self, oliglist, br=50, rr=750):
        self.oligos = self.remove_size(oliglist)
        self.oligos = self.remove_name(oliglist)
        self.oligos = zip([item[0] for item in self.oligos], [item[1].upper() for item in self.oligos], [item[1].upper()[::-1] for item in self.oligos])
        
        self.binders = []
        self.buff_range = br
        self.read_range = rr
        
    def remove_size(self, oliglist):
        return [x for x in oliglist if len(x[1]) < 35 and len(x[1]) > 15]
    
    def remove_name(self, oliglist):
        return [x for x in oligolist if len(x[0]) <= 13]
        
    def find_bind(self):
    	#Would be nice to also find binding spots that are off by 1 nucleotide
        for i, o in enumerate(self.oligos):
            bind = [m.end() for m in re.finditer(o[1], plasmid.sequence)]
            bind_rev = [m.start() for m in re.finditer(o[2], plasmid.comp)]
        
            if len(bind) + len(bind_rev) == 1:
                direction, ran = self.pattern(bind, bind_rev)
                self.binders.append((o[0], o[1], o[2], direction, ran))
        
    def pattern(self, b, br):
        if b == []:
            return '-', [br[0]-self.read_range, br[0]-self.buff_range]
        elif br == []:
            return '+', [b[0]+self.buff_range, b[0]+self.read_range]





with open('Export_-_2018-01-17/Oligos_-_2018-01-17.csv', 'rb') as f:
    reader = csv.reader(f)
    sjoligos = list(reader)
    
with open('Export_-_2018-01-17/My_Oligos_-_2018-01-17.csv', 'rb') as f:
    reader = csv.reader(f)
    kltkoligos = list(reader)
    
with open('Export_-_2018-01-17/Elim_Primers_-_2018-01-17.csv', 'rb') as f:
    reader = csv.reader(f)
    elimoligos = list(reader)
    
with open('Export_-_2018-01-17/Sequetech_Primers_-_2018-01-17.csv', 'rb') as f:
    reader = csv.reader(f)
    seqoligos = list(reader)

elim = input('Use ELIM primers (0, 1)? ')
if elim == 1:
    alloligolist = sjoligos + kltkoligos + elimoligos + seqoligos
elif elim == 0:
    alloligolist = sjoligos + kltkoligos + seqoligos

names = [item[0] for item in alloligolist]
sequences = [item[1] for item in alloligolist]
oligolist = zip(names, sequences)

adj = input('Adjust primer ranges (0, 1)? ')
if adj == 1:
    buff = input('Buffer range (default = 50)? ')
    ran = input('Read range (default = 750)? ')
    
    o = Oligos(oligolist, buff, ran)
else:
    o = Oligos(oligolist)


with open('plasmid.txt', 'r') as f:
    plasmid_sequence = f.read()

plasmid = Plasmid(plasmid_sequence.upper())
range_start = input('Start of range? ')
range_end = input('End of range? ')
plasmid.range_int(range_start,range_end)
o.find_bind()
plasmid.copy_pattern(o.binders)
plasmid.seq_finder()