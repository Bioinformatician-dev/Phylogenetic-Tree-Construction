from Bio import AlignIO
from Bio.Phylo import TreeConstruction
from Bio import Phylo
import matplotlib.pyplot as plt

# Load the genetic sequences
alignment = AlignIO.read("sequences.fasta", "fasta")

# Construct a distance matrix using the Neighbor-Joining method
calculator = TreeConstruction.DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

# Create a Neighbor-Joining tree
constructor = TreeConstruction.NeighborJoining()
nj_tree = constructor.build_tree(distance_matrix)

# Visualize the tree
Phylo.draw(nj_tree)
plt.title("Neighbor-Joining Phylogenetic Tree")
plt.show()
