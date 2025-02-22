# Phylogenetic Analysis using MEGA
The goal of this project is to understand the basic concepts of phylogenetic analysis. Phylogenetic analysis refers to the study of evolutionary relationships among organisms, genes, or proteins based on their genetic sequences. 
To identify the common ancestors and understand how species have evolved over time.

#Phylogenetic trees: 
A phylogenetic tree is constructed to represent evolutionary relationships.
There are two tree construction methods.
1) Distance-based methods (Neighbor-Joining method, Unweighted Pair Group Method with Arithmetic mean (UPGMA))
2) Character-based methods (Maximum Likelihood method, Maximum Parsimony method)

# Neighbor-Joining Method:
# UPGMA Method
It assumes constant rate of evolution (molecular clock) means all species in the tree have evolved for the same duration from their common ancestor.
For example, we have some sequences and we need to construct a phylogenetic tree using these sequences and UPGMA method.
Sequence A: ATCGTGGTACTG
Sequence B: CCGGAGAACTAG
Sequence C: AACGTGCTACTG
Sequence D: ATGGTGAAAGTG
Sequence E: CCGGAAAACTTG
Sequence F: TGGCCCTGTATC

Step 1: Compute the distance matrix (In matrix we calculate the paiwise genetic distances or number of mismatches between two sequences).
   A   B   C   D   E   F
A  0   9   2   4   9   11
B  0   0   9   6   2   11
C  0   0   0   5   9   11
D  0   0   0   0   6   10
E  0   0   0   0   0   10
F  0   0   0   0   0   0

Step 2: Identify the closest pair (pair with smallest distance)
In this example, A and C, B and E are the closest because number of mismatches are smallest for these.

Step 3: Merge A and C, B and E into a new cluster and calculate distances between remaining taxa and the new clusters using the before formula.

d(AC,X) = d(A,X)+d(C,X)/2


      AC   B   D   E   F
AC   0     9   4.5   9   11
B     0     0   6   2   11
D     0     0   0   6   10
E     0     0   0   0   10
F     0     0   0   0   0




# Link to install MEGA (https://www.megasoftware.net/releases/MEGA_12.0.8_win64_setup.exe)
