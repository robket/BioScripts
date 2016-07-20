from __future__ import division
import sys
import string
import math
import numpy as np


class Settings():
    def __init__(s):
        s.INDEL_LENGTH_MEAN = 21
        s.INDEL_LENGTH_STD = 6
        s.INDELS_MULTIPLES_OF_3 = True
        s.MINIMUM_REGION_LENGTH = 5
        s.ALLOW_INDEL_IN_IN = False

        # Simulation poperties
        s.NUMBER_OF_GENERATIONS = 4
        s.NUMBER_OF_OFFSPRING_PER_GENERATION = 3
        s.KEEP_PARENTS = True
        s.ALLOW_STOP_CODON_MUTATION = False
        s.ACTIVELY_DELETE_SEQUENCES_WITH_STOPS = True
        s.SHUFFLE_OUTPUT = True
        s.MAX_POOL_SIZE = 1024
        s.NUMBER_OF_REGIONS = 5
        s.CONSERVED_REGION_SIZE = 200
        s.VARIABLE_REGION_SIZE = 60
        s.CONSERVED_MUTATE = 0.001
        s.VARIABLE_MUTATE = 0.003
        s.CONSERVED_INDEL = 0.0
        s.VARIABLE_INDEL = 0.0005
        length = math.ceil(s.NUMBER_OF_REGIONS/2) * s.CONSERVED_REGION_SIZE + math.floor(s.NUMBER_OF_REGIONS/2) * s.VARIABLE_REGION_SIZE
        # Recombinations
        s.PERCENTAGE_OFFSPRING_THAT_ARE_RECOMBINATIONS = 0.1
        s.CONSERVED_RECO = 2.0/length
        s.VARIABLE_RECO = 0.0/length
        s.SWITCH_COUNT = 1

def true_or_false(probability):
    return np.random.rand(1) < probability

class Nucleotide:
    NUCLEOTIDES = set(["A", "C", "T", "G"])
    @staticmethod
    def mutate(nucleotide, mutation_rate):
        do_mutate =  true_or_false(mutation_rate)
        if do_mutate:
            return np.random.choice(list(Nucleotide.NUCLEOTIDES - set(nucleotide)))
        else:
            return nucleotide

class Strain:
    class Region:
        def __init__(self, start, length, mutation_rate, indel_rate, reco_rate):
            self.start = start
            self.length = length
            self.end = start + length
            self.mutation_rate = mutation_rate # The likelihood of individual nucleotide mutating to a different nucleotide
            self.indel_rate = indel_rate # The likelihood that an individual nucleotide will be the start of an insertion/deletion
            self.reco_rate = reco_rate

        def clone(self, start=None, length=None):
            if start is None:
                start = self.start
            if length is None:
                length = self.length
            return Strain.Region(start, length, self.mutation_rate, self.indel_rate, self.reco_rate)

    def __init__(self, sequence, regions, *parents):
        self.sequence = sequence
        self.contains_stop_codons = Strain.first_stop_codon_index(sequence) is not None
        self.regions = regions
        if len(parents) == 0:
           self.parents = None
        else:
           self.parents = parents

    @staticmethod
    def first_stop_codon_index(sequence, start_index=0):
        for index in range(start_index, len(sequence) - 2, 3):
            if sequence[index] == "T" and (sequence[index + 1] == "A" and (sequence[index + 2] == "A" or sequence[index + 2] == "G")
                                        or sequence[index + 1] == "G" and sequence[index + 2] == "A"):
                return index
        return None

    def fix_stop_codons(self):
        next_stop = Strain.first_stop_codon_index(self.sequence)
        while next_stop is not None:
            self.sequence[next_stop] = np.random.choice(list("ACG"))
            next_stop = Strain.first_stop_codon_index(self.sequence, next_stop + 3)
        self.contains_stop_codons = False

    def generate_offspring(self):
        mutated_strain = self._mutate()
        indeled_strain = mutated_strain._indel()
        return indeled_strain

    def _mutate(self):
        new_sequence = []
        new_regions = []
        for region in self.regions:
            for i in range(region.start, region.end):
                try:
                    new_sequence.append(Nucleotide.mutate(self.sequence[i], region.mutation_rate))
                except:
                    print("Error: %d"%(i))
            new_regions.append(region.clone())
        return Strain(new_sequence, new_regions, self)

    def _indel(self):
        new_sequence = []
        new_regions = []
        i = 0
        for region in self.regions:
            new_region_start = len(new_sequence)
            expected_region_length = region.length
            max_i = 0
            while i < region.length:
                max_i = max(i, max_i)
                if not s.ALLOW_INDEL_IN_IN and i < max_i:
                    is_indel_start = False # Do not start insertions/deletions while busy with an insertion/deletion
                else:
                    is_indel_start = true_or_false(region.indel_rate)
                if not is_indel_start:
                    new_sequence.append(self.sequence[max(0, region.start + i)])
                    i += 1
                    continue
                insertion_not_deletion = true_or_false(0.5)
                indel_length = abs(int(np.random.standard_normal() * s.INDEL_LENGTH_STD + s.INDEL_LENGTH_MEAN))
                if s.INDELS_MULTIPLES_OF_3:
                    indel_length = round(indel_length / 3) * 3
                if insertion_not_deletion:
                    # Go Back length letters, insert them again
                    i -= indel_length
                    i = max(i, 0 - region.start) # Cannot go back before the start of the sequence
                    expected_region_length += indel_length
                else:
                    # Check that this deletion won't bring the total region length <  MINIMUM_REGION_LONGTH
                    if expected_region_length - indel_length > s.MINIMUM_REGION_LENGTH:
                        i += indel_length
                        expected_region_length -= indel_length
                    else:
                        pass # Don't delete any
            i -= region.length # If there are still bases to be skipped (i > region.length), skip over the start of the next region
            new_regions.append(region.clone(new_region_start, len(new_sequence) - new_region_start))
        return Strain(new_sequence, new_regions, self)

    @staticmethod
    def get_random_recombine_positions(regions, N):
        region_probs = np.array([region.length * region.reco_rate for region in regions])
        region_probs /= sum(region_probs)
        positions = []
        for i in range(N):
            chosen_region = np.random.choice(regions, p=region_probs)
            position_in_region = np.random.randint(chosen_region.length)
            positions.append(chosen_region.start + position_in_region)
        positions.sort()
        return positions

    @staticmethod
    def recombine(strain1, strain2):
        active = 0
        active_strain = strain1
        switch_count = 0
        new_sequence = []
        new_regions = []
        if len(strain1.regions) != len(strain2.regions):
           raise ValueError("regions1 and regions2 should be similarly sized arrays")
        recombine_positions = Strain.get_random_recombine_positions(strain1.regions, s.SWITCH_COUNT)
        for region_index in range(len(strain1.regions)):
            active_region = active_strain.regions[region_index]
            i = 0
            new_region_start = len(new_sequence)
            while i < active_region.length:
                switch = switch_count < len(recombine_positions) and i + active_region.start == recombine_positions[switch_count]
                if switch:
                    print("Switched at %d"%(active_region.start + i), file=sys.stderr)
                    switch_count += 1
                    current_percent = i / active_region.length
                    active_strain = strain2 if active_strain == strain1 else strain1
                    active_region = active_strain.regions[region_index]
                    i = int(round(current_percent * active_region.length))
                new_sequence.append(active_strain.sequence[active_region.start + i])
                i += 1
            new_regions.append(active_region.clone(new_region_start, len(new_sequence) - new_region_start))
        print("Switched %d times"%(switch_count), file=sys.stderr)
        offspring = Strain(new_sequence, new_regions, strain1, strain2)
        return  offspring

def random_strain(number_of_regions, conserved_region_size, variable_region_size, conserved_mutate, variable_mutate, conserved_indel, variable_indel, conserved_reco, variable_reco, alphabet=list(Nucleotide.NUCLEOTIDES)):
    regions = []
    size = 0
    for i in range(number_of_regions):
        mutate_rate = conserved_mutate if i%2 == 0 else variable_mutate
        indel_rate = conserved_indel if i%2 == 0 else variable_indel
        reco_rate = conserved_reco if i%2 == 0 else variable_reco
        diff = conserved_region_size if i%2 == 0 else variable_region_size
        region = Strain.Region(size, diff, mutate_rate, indel_rate, reco_rate)
        size += diff
        regions.append(region)
    sequence = np.random.choice(alphabet, size=[size])
    return Strain(sequence, regions)

def simulation(starting_pool):
    pool = starting_pool
    for g in range(s.NUMBER_OF_GENERATIONS):
        print("Generation %d, population size %d"%(g, len(pool)) ,file=sys.stderr)
        offsprings = []
        for strain in pool:
            if s.KEEP_PARENTS:
                offsprings.append(strain)
            if not s.ALLOW_STOP_CODON_MUTATION and strain.contains_stop_codons:
                print("Strain %s has stop codons, no offspring for you." % (strain.id), file=sys.stderr)
                continue
            for i in range(s.NUMBER_OF_OFFSPRING_PER_GENERATION):
                should_recombine = true_or_false(s.PERCENTAGE_OFFSPRING_THAT_ARE_RECOMBINATIONS)
                if should_recombine:
                    other_strain = np.random.choice(pool)
                    offspring = Strain.recombine(strain, other_strain)
                    offspring.id = "-".join([s + o for s, o in zip(strain.id.split("-"), other_strain.id.split("-"))])
                else:
                    offspring = strain.generate_offspring()
                    offspring.id = strain.id + "-" + string.ascii_uppercase[i]
                if not s.ACTIVELY_DELETE_SEQUENCES_WITH_STOPS or not offspring.contains_stop_codons:
                    offsprings.append(offspring)
        # Choosing the next population
        if s.MAX_POOL_SIZE is not None and len(offsprings) > s.MAX_POOL_SIZE:
            indexes = list(np.random.choice(range(len(offsprings)), size=s.MAX_POOL_SIZE, replace=False))
            if not s.SHUFFLE_OUTPUT:
                indexes.sort()
        elif s.SHUFFLE_OUTPUT:
            indexes = list(range(len(offsprings)))
            np.random.shuffle(indexes)
        else:
            indexes = range(len(offsprings))

        pool = [offsprings[i] for i in indexes]
    # Check offsprings for stop codons
    return pool

def print_strain(strain, detailed_regions=False):
    if detailed_regions:
        print(" ".join(["%d -(%d)-> %d"%(r.start, r.length, r.end) for r in strain.regions]))
    print(" ".join(["".join(strain.sequence[r.start:r.end]) for r in strain.regions]))

def test_mutate():
    print("Test Mutate")
    starting_strain = random_strain(5, 40, 40, 0.05, 0.025, 6, 0.1, 0.1)
    print_strain(starting_strain)
    mutated_strain = starting_strain._mutate()
    print_strain(mutated_strain)

def test_indel():
    print("Test Indel")
    starting_strain = random_strain(5, 40, 40, 0.05, 0.025, 6, 0.1, 0.1, alphabet=list(string.ascii_letters + string.digits))
    print_strain(starting_strain, detailed_regions=True)
    mutated_strain = starting_strain._indel()
    print_strain(mutated_strain, detailed_regions=True)

def test_recombination():
    print("Test Recombination")
    strain1 = random_strain(5, 40, 40, 0.05, 0.025, 6, 0.1, 0.1)
    strain2 = Strain([s.lower() for s in strain1.sequence], strain1.regions)
    print_strain(strain1)
    print_strain(strain2)
    mutated_strain = Strain.recombine(strain1, strain2)
    print_strain(mutated_strain)

np.random.seed(0)

s = Settings()
# Initial random starting strain
#s.CONSERVED_MUTATE = 0.0005
#s.VARIABLE_MUTATE = 0.03
#s.CONSERVED_INDEL = 0.0
#s.VARIABLE_INDEL = 0.002

starting_strain = random_strain(s.NUMBER_OF_REGIONS, s.CONSERVED_REGION_SIZE, s.VARIABLE_REGION_SIZE, s.CONSERVED_MUTATE, s.VARIABLE_MUTATE, s.CONSERVED_INDEL, s.VARIABLE_INDEL, s.CONSERVED_RECO, s.VARIABLE_RECO, alphabet=list(Nucleotide.NUCLEOTIDES))
starting_strain.fix_stop_codons()
starting_strain.id = str(0)

s.KEEP_PARENTS = False
s.PERCENTAGE_OFFSPRING_THAT_ARE_RECOMBINATIONS = 0
s.MAX_POOL_SIZE = 1
s.NUMBER_OF_OFFSPRING_PER_GENERATION = 1
s.ALLOW_STOP_CODON_MUTATION = True
s.ACTIVELY_DELETE_SEQUENCES_WITH_STOPS = False
s.NUMBER_OF_GENERATIONS = 8
seeding_pool = [starting_strain]

for i in range(1, 4):
    starting_strain = simulation([starting_strain])[0]
    starting_strain.id = str(i)
    starting_strain.fix_stop_codons()
    seeding_pool.append(starting_strain)

s = Settings()
for i, strain in enumerate(seeding_pool):
    for i, region in enumerate(strain.regions):
        region.mutation_rate = s.CONSERVED_MUTATE if i%2 == 0 else s.VARIABLE_MUTATE
        region.indel_rate = s.CONSERVED_INDEL if i%2 == 0 else s.VARIABLE_INDEL
        region.reco_rate = s.CONSERVED_RECO if i%2 == 0 else s.VARIABLE_RECO

# Running the real simulation given the initial strains
pool = simulation(seeding_pool)

# Print output in fasta format
for i, strain in enumerate(pool):
    print("> Strain %d_%s"%(i, strain.id))
    print("".join(strain.sequence))
