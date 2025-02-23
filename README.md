# Phylogenetic Analysis using MEGA
The goal of this project is to understand the basic concepts of phylogenetic analysis. Phylogenetic analysis refers to the study of evolutionary relationships among organisms, genes, or proteins based on their genetic sequences. 
To identify the common ancestors and understand how species have evolved over time.

#Phylogenetic trees: 
A phylogenetic tree is constructed to represent evolutionary relationships.
There are two tree construction methods.
1) Distance-based methods (Neighbor-Joining method, Unweighted Pair Group Method with Arithmetic mean (UPGMA))
2) Character-based methods (Maximum Likelihood method, Maximum Parsimony method)

# Neighbor-Joining Method:
Key points:
Starts from a distance matrix.
Pairs up the most closely related taxa.
Treat this pair as new taxa.
#Compute the distance matrix (In matrix we calculate the paiwise genetic distances or number of mismatches between two sequences).
d  A   B   C   D
A  0   17  21  27
B  0   0   12  18
C  0   0   0   14
D  0   0   0   0

#Round 1 
#Step 1: Compute r'i for each terminal node using the below formula (i = A/B/C/D).

r'i = Summation(Di,j)/(n-2)      # n = number of taxa, i and j are two taxa.
r'A = D(A,B)+D(A,C)+D(A,D)/(4-2) = 17+21+27 = 32.5
Similary compute r'i for each node.
d  A   B   C   D   r'
A  0   17  21  27  32.5
B  0   0   12  18  23.5
C  0   0   0   14  23.5
D  0   0   0   0   29.5

#Step 2: Compute d'i,j for each terminal node using the below formula.
d'ij = dij-r'i-r'j
d'AB = dAB-r'A-r'B = 17-32.5-23.5 = -39


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

Step 3: Merge A and C  into a new cluster and calculate distances between remaining taxa and the new cluster using the below formula.

d(AC,X) = d(A,X)+d(C,X)/2
For example, calculate distance between AC and B.
d(AC,B) = d(A,B)+d(C,B)/2 = 9+9/2 = 9

      AC   B   D   E   F
AC   0     9   4.5   9   11
B     0     0   6   2   11
D     0     0   0   6   10
E     0     0   0   0   10
F     0     0   0   0   0

Step 4: Merge B and E into a new cluster and calculate distance between remaining taxa and this cluster using the above formula.

   AC   BE   D   F
AC  0   9   4.5   11
BE  0   0   6     10.5
D   0   0   0     10
F   0   0   0     0

#Identify the taxon and the cluster which are more close to each other. For example, taxon D and cluster AC are more close to each other so merge them into a new cluster.

Step 5: Calculate distance between remaining taxon F and clusters ACD and BE using the above formula.

   ACD   BE   F
ACD  0   7.5  10.5
BE   0   0    10.5
F    0   0    0
#Taxon F is equidistant from cluster ACD and BE so it will be a new cluster. Cluster ACD and BE will be closer to each other than cluster F because distance between them is less than cluster F.

# Link to install MEGA (https://www.megasoftware.net/releases/MEGA_12.0.8_win64_setup.exe)
