from Bio import SeqIO
from Bio import AlignIO
from Bio.Phylo import DistanceTreeConstructor, TreeConstruction
from Bio import Phylo
import matplotlib.pyplot as plt

# Function to read DNA sequences from a FASTA file
def read_fasta(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

# Function to compute distance matrix
def compute_distance_matrix(sequences):
    # Create a MultipleSeqAlignment object
    alignment = AlignIO.MultipleSeqAlignment([SeqIO.SeqRecord(seq, id=seq_id) for seq_id, seq in sequences.items()])
    
    calculator = TreeConstruction.DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    return distance_matrix

# Function to construct a phylogenetic tree using Neighbor-Joining
def construct_tree(distance_matrix):
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)
    return tree

# Function to visualize the phylogenetic tree
def visualize_tree(tree):
    Phylo.draw(tree)
    plt.title("Phylogenetic Tree")
    plt.show()

# Main function to integrate all components
def main(file_path):
    sequences = read_fasta(file_path)
    print("Sequences read:")
    for seq_id in sequences:
        print(f"{seq_id}: {sequences[seq_id]}")

    # Compute the distance matrix
    distance_matrix = compute_distance_matrix(sequences)
    
    # Construct the phylogenetic tree
    tree = construct_tree(distance_matrix)
    
    # Visualize the tree
    visualize_tree(tree)

if __name__ == "__main__":
    fasta_file = "your_sequences.fasta"  # Replace with your FASTA file path
    main(fasta_file)
