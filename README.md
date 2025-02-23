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
Similarly compute d'ij for each node.
d'   A   B   C   D
A    0   -39 -35 -35
B    0   0   -35 -35
C    0   0   0   -39
D    0   0   0    0

#Select minimum d' (For e.g. minimum distance is -39 between A&B and C&D, we can select and join any of these two taxa).
First we consider joining of taxa C and D to a node u.

#Step 3: Calculate branch lengths.
Vi = 0.5(dij)+0.5(r'i-r'j)
Vj = 0.5(dji)+0.5(r'j-r'i)

VCu = 0.5(dCD)+0.5(r'C-r'D) = 0.5(14)+0.5(23.5-29.5) = 4
VDu = 0.5(dDC)+0.5(r'D-r'C) = 0.5(14)+0.5(29.5-23.5) = 10

Length of branches Cu and Du are 4 and 10 respectively.

#Round 2
#Step 1:

dij,k = (di,k+dj,k-di,j)/2   (ij = CD, k = A/B)
dCD,A = (dCA+dDA-dcD)/2 = (21+27-14)/2 = 17
dCD,B = (dCB+dDB-dcD)/2 = (12+18-14)/2 = 8

r'A = (dAB+dACD)/(n-2) = (17+17)/(3-2) = 34
Similary calculate r'B and r'CD.

d   A   B   CD    r'i
A   0   17  17    34
B   0   0   8     25
CD  0   0   0     25

#Step 2: Compute d'i,j for new matrix using above formula

d'   A   B   CD
A    0  -42 -42
B    0   0  -42
CD   0   0   0

#Select minimum d' (minimum distance is -42 so we can join any of three for e.g, A to B, A to CD, B to CD)
We join A to B to a node.
#Step 3: Calculate branch length using above formula.
VAV = 0.5(dAB)+0.5(r'A-r'B) = 0.5(17)+0.5(34-25) = 13
VBv = 0.5(dBA)+0.5(r'B-r'A) = 0.5(17)+0.5(25-34) = 4

Length of branches Av and Bv are 13 and 4 respectively.

#Round 3
#Step 1: Update distance matrix and calculate distance between AB and CD.
d   AB   CD
AB  0     4
CD  0     0


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
