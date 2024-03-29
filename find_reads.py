import csv
import re
from tabulate import tabulate
from glob import glob
#import os
#will need it for window resizing
from os.path import join


def complementary(seq):

    test_seq = seq.replace('G', 'c')
    test_seq = test_seq.replace('C', 'G')
    test_seq = test_seq.replace('c', 'C')
    test_seq = test_seq.replace('A', 't')
    test_seq = test_seq.replace('T', 'A')
    test_seq = test_seq.replace('t', 'T')
    return test_seq

# ADD FUNCTIONALITY FOR TERMINAL WINDOW RESIZE

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
        self.read = []

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

        self.target_range = [st - 1, en]
        self.read_length = en - st

    def copy_pattern(self, oligos):
        # Find each primers range
        # Odd = start, even = end

        count = 1
        for i, o in enumerate(oligos):
            if self.blank_range[o[4][0]] == 0:
                self.blank_range[o[4][0]] = count
                self.oligo_list[count] = [o[0], o[3]]
                self.oligo_locs[count] = [o[4][0], o[4][1], o[5]]

                count += 1
                self.blank_range[o[4][1]] = count
                count += 1
            else:
                new_odd, new_even = self.adjuster(o[4][0], o[4][1])

                self.blank_range[new_odd] = count
                self.oligo_list[count] = [o[0], o[3]]
                self.oligo_locs[count] = [new_odd, new_even, o[5]]

                count += 1
                self.blank_range[new_even] = count
                count += 1

    def adjuster(self, st, en):
        # If one of the original locations is already occupied
        mag = 1
        vec = 1

        while True:
            new_st = st + (mag * vec)
            new_en = en + (mag * vec)

            if self.blank_range[new_st] == 0 and self.blank_range[new_en] == 0:
                return new_st, new_en
            elif vec > 0:
                vec *= -1
            else:
                vec *= -1
                mag += 1

    def seq_finder(self):

        self.find_full_seq(self.target_range[0], 1, [])

        if len(self.complete_reads) >= 1:
            print ""
            print "*** COMPLETE READ FOUND ***"

            basesper = self.read_length / len(self.complete_reads[0])
            print ""
            print "Read length: " + str(self.read_length)
            print "Avg bases per primer: " + str(basesper)
            print ""

            self.print_reads(self.complete_reads)
            self.replacer()
        else:
            self.no_path()

    def replacer(self):

        print ""
        rep = input('Replace an oligo (input Oligo ID, 0 to exit)? ')

        if rep == 0:
            return
        elif rep % 2 == 1:
            try:
                plasmid.rem_oligo(rep)
                plasmid.num_iters = 0
                plasmid.complete_reads = []
                plasmid.seq_finder()
            except:
                print "Oligo not found"
                self.replacer()
        elif rep % 2 == 0:
            print "Oligo not found"
            self.replacer()

    def find_last_odd(self, curr, lowlim):

        # NEED TO UPDATE TO INCLUDE PRIMERS THAT START BEFORE BASE 0 (without causing recursive error)
        # UPDATE TO INCLUDE PRIMER PREFERENCES

        all_rem_odds = [x for x in self.blank_range[lowlim+1:curr] if x % 2 == 1]

        if all_rem_odds == []:
            return ['', '']

        last_odd = all_rem_odds[-1]
        new_curr = plasmid.blank_range.index(last_odd)
        return [last_odd, new_curr]

    def find_full_seq(self, curr_ind, low_ind, path):

        self.num_iters += 1
        if self.num_iters == 50:
            return

        [next_odd, next_low] = self.find_last_odd(curr_ind, low_ind)

        if next_odd != '':

            path.append(next_odd)
            curr_ind = self.blank_range.index(next_odd + 1)

            if curr_ind >= self.target_range[1]:
                self.complete_reads.append(list(path))
                return

            elif curr_ind <= self.blank_range.index(next_odd):
                self.complete_reads.append(list(path))
                return

            else:
                self.find_full_seq(curr_ind, next_low, path)

        else:
            if path not in self.partial_reads:
                self.partial_reads.append(list(path))
            return

    def rem_oligo(self, oli):

        self.blank_range[self.blank_range.index(oli)] = 0
        self.blank_range[self.blank_range.index(oli + 1)] = 0

    def print_reads(self, readlist):

        for read in readlist:

            headers = ['Oligo ID', 'Oligo Name', 'Direction', 'Binding Location', 'Approx. Read Range', 'Overlap']

            table = []
            for item in read:
                oligoid = str(item)
                oligoname = self.oligo_list[item][0]
                direction = self.oligo_list[item][1]
                bindloc = self.oligo_locs[item][2][0]
                readran = [self.oligo_locs[item][0], self.oligo_locs[item][1]]

                if readran[0] < 0:
                    readran[0] = len(self.sequence) + readran[0]

                if len(table) >= 1:
                    overlap = table[-1][-2][-1] - readran[0]
                else:
                    overlap = '-'

                table.append([oligoid, oligoname, direction, bindloc, readran, overlap])

            print tabulate(table, headers=headers)


    def find_empty_ranges(self):

        # FIX OUT OF RANGE ERROR THAT I GET ON THIS INDEX
        last_end = self.blank_range.index(self.partial_reads[-1][-1] + 1)

        if last_end <= self.blank_range.index(self.partial_reads[-1][-1]):
            return

        next_start = self.find_next_odd(last_end)

        if next_start >= self.target_range[1] or next_start <= self.target_range[0]:
            return

        self.find_full_seq(next_start + 1, last_end, [])

        if len(self.complete_reads) == 1:
            self.partial_reads.append(list(self.complete_reads[0]))
            return
        else:
            self.find_empty_ranges()

    def find_next_odd(self, curr):

        future_odds = [x for x in self.blank_range[curr:] if x % 2 == 1]

        if future_odds == []:
            future_odds = [x for x in self.blank_range if x % 2 == 1]

        return self.blank_range.index(future_odds[0])

    def coverage(self):

        for p in self.partial_reads:
            st = self.blank_range.index(p[0])
            en = self.blank_range.index(p[-1] + 1)
            self.read.append(list([st, en]))

    def no_path(self):

        print ""
        print '*** NO COMPLETE READ FOUND ***'

        try:
            self.find_empty_ranges()
            self.coverage()
            print "These oligos will sequence part of the target range"
            self.print_reads(self.partial_reads)
            print " "
            print "The following ranges are covered by the above: "
            print " "
        except:
            print "No useful oligo binding locations found."

        for r in self.read:
            if r[1] < r[0]:
                print "Pos " + str(r[0]) + " to End"
            else:
                print "Pos " + str(r[0]) + " to Pos " + str(r[1])

        print ""
        return


class Oligos:

    def __init__(self, oliglist, br=36, rr=850):

        self.oligos = self.remove_size(oliglist)
        self.oligos = self.remove_name(self.oligos)

        #bad
        try:
            self.oligos = self.diluted_primers(self.oligos)
        except:
            pass

        self.oligos = zip([item[0] for item in self.oligos],
                          [item[1].upper() for item in self.oligos],
                          [item[1].upper()[::-1] for item in self.oligos])

        self.binders = []
        self.buff_range = br
        self.read_range = rr

    def remove_size(self, oliglist):

        return [x for x in oliglist if len(x[1]) < 35 and len(x[1]) >= 18]

    def remove_name(self, oliglist):

        return [x for x in oliglist if len(x[0]) <= 13]

    def find_bind(self):

        # Would be nice to also find binding spots that are off by 1 nucleotide
        for i, o in enumerate(self.oligos):
            bind = [m.end() for m in re.finditer(o[1], plasmid.sequence)]
            bind_rev = [m.start() for m in re.finditer(o[2], plasmid.comp)]

            if len(bind) + len(bind_rev) == 1:
                direction, ran = self.pattern(bind, bind_rev)
                self.binders.append((o[0], o[1], o[2], direction, ran, max(bind, bind_rev)))

    def pattern(self, b, br):
        # READ RANGE SHOULD BE UPDATED TO READ LENGTH

        if b == []:
            return '-', [br[0] - self.read_range, br[0] - self.buff_range]
        elif br == []:
            l = len(plasmid.sequence)

            if b[0] + self.buff_range >= l:
                return '+', [(b[0] + self.buff_range) - l, (b[0] + self.read_range) - l]
            elif b[0] + self.read_range >= l:
                return '+', [b[0] + self.buff_range, (b[0] + self.read_range) - l]
            else:
                return '+', [b[0] + self.buff_range, b[0] + self.read_range]

    def diluted_primers(self, oliglist):

        '''
        Compares list of oligo names to list of diluted primers in csv and stars the primers that already have dilutions made
        This is to make manual sorting easier, but in the future, would be nice if this also fed into the read finding algorithm
        '''

        with open('sequencing_primers.csv', 'rb') as sheet:
            reader = csv.reader(sheet)
            diluted_primers = list(reader)

        diluted_primers = [d for p in diluted_primers for d in p]

        new_list = []

        for x in oliglist:
            if x[0] in diluted_primers:
                star_name = x[0] + '*'
                new_list.append((star_name, x[1]))
            else:
                new_list.append(x)

        return new_list

files = glob(join('*Export*', '*'))

otheroligos = []
for f in files:

    with open(f,'rb') as sheet:
        reader = csv.reader(sheet)

        if 'Sequetech' in f:
            seqoligos = list(reader)
        elif 'Elim' in f:
            elimoligos = list(reader)
        elif 'Genewiz' in f:
            genewizoligos = list(reader)
        else:
            l = list(reader)
            otheroligos += l

# PRIMER REMOVAL BASED ON NAME AND LENGTH SHOULD HAPPEN HERE
print ""
# ADD QUINTARA PRIMERS
company = input('Comapny primers? (1: Sequetech, 2: ELIM, 3: Genewiz, 4: None) ')
if company == 3:
    alloligolist = otheroligos + genewizoligos
elif company == 1:
    alloligolist = otheroligos + seqoligos
elif company == 2:
    alloligolist = otheroligos + elimoligos
else:
    alloligolist = otheroligos


names = [item[0] for item in alloligolist]
sequences = [item[1] for item in alloligolist]
oligolist = zip(names, sequences)

adj = input('Adjust primer ranges (0, 1)? ')
if adj == 1:
    buff = input('Buffer range (default = 36)? ')
    ran = input('Read range (default = 950)? ')

    o = Oligos(oligolist, buff, ran)
else:
    o = Oligos(oligolist)

print ""
print "Total oligos found: " + str(len(oligolist))
print "Usable oligos found: " + str(len(o.oligos))


with open('plasmid.txt', 'r') as f:
    plasmid_sequence = f.read()
    plasmid_sequence = plasmid_sequence.strip()

plasmid = Plasmid(plasmid_sequence.upper())

print ""
print "Plasmid length: " + str(len(plasmid.sequence)) + " bp"
print ""

range_start = input('Start of range: ')
range_end = input('End of range: ')
plasmid.range_int(range_start, range_end)
o.find_bind()

plasmid.copy_pattern(o.binders)
plasmid.seq_finder()
