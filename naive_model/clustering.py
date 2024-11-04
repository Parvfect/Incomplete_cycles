

import numpy as np
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq

class ConsensusSelector:
    def __init__(self, quality_threshold=30):
        self.quality_threshold = quality_threshold

    def select_final_consensus(self, cluster_sequences, quality_scores=None, method='weighted'):
        """
        Select final consensus sequence using specified method
        
        Parameters:
        cluster_sequences (list): List of DNA sequences in cluster
        quality_scores (list): Optional list of quality scores for each sequence
        method (str): Method to use ('majority', 'weighted', or 'quality_filtered')
        
        Returns:
        str: Final consensus sequence
        float: Confidence score
        """
        if method == 'majority':
            return self._majority_voting(cluster_sequences)
        elif method == 'weighted':
            return self._weighted_consensus(cluster_sequences, quality_scores)
        elif method == 'quality_filtered':
            return self._quality_filtered_consensus(cluster_sequences, quality_scores)
        else:
            raise ValueError("Unknown consensus method")

    def _majority_voting(self, sequences):
        """Simple majority voting at each position"""
        # Pad sequences to same length
        max_length = max(len(seq) for seq in sequences)
        padded_seqs = [seq.ljust(max_length, '-') for seq in sequences]
        
        consensus = ''
        confidence_scores = []
        
        # Vote at each position
        for i in range(max_length):
            bases = [seq[i] for seq in padded_seqs]
            base_counts = Counter(bases)
            
            # Remove gaps from consideration if possible
            if len(base_counts) > 1 and '-' in base_counts:
                del base_counts['-']
            
            # Get most common base and its frequency
            most_common = base_counts.most_common(1)[0]
            consensus += most_common[0]
            
            # Calculate confidence as fraction of votes for winning base
            confidence = most_common[1] / len(sequences)
            confidence_scores.append(confidence)
        
        return consensus.replace('-', ''), np.mean(confidence_scores)

    def _weighted_consensus(self, sequences, quality_scores):
        """Weight each sequence's vote by its quality score"""
        if not quality_scores:
            raise ValueError("Quality scores required for weighted consensus")
            
        # Convert quality scores to weights
        weights = np.array(quality_scores) / np.sum(quality_scores)
        max_length = max(len(seq) for seq in sequences)
        padded_seqs = [seq.ljust(max_length, '-') for seq in sequences]
        
        consensus = ''
        confidence_scores = []
        
        for i in range(max_length):
            bases = [seq[i] for seq in padded_seqs]
            base_weights = {}
            
            # Sum weights for each base
            for base, weight in zip(bases, weights):
                if base not in base_weights:
                    base_weights[base] = 0
                base_weights[base] += weight
            
            # Remove gaps if possible
            if len(base_weights) > 1 and '-' in base_weights:
                del base_weights['-']
            
            # Select base with highest weighted vote
            selected_base = max(base_weights.items(), key=lambda x: x[1])
            consensus += selected_base[0]
            confidence_scores.append(selected_base[1])
        
        return consensus.replace('-', ''), np.mean(confidence_scores)

    def _quality_filtered_consensus(self, sequences, quality_scores):
        """Filter low quality reads before consensus"""
        if not quality_scores:
            raise ValueError("Quality scores required for quality filtering")
            
        # Filter sequences based on quality threshold
        filtered_data = [
            (seq, score) for seq, score in zip(sequences, quality_scores)
            if score >= self.quality_threshold
        ]
        
        if not filtered_data:
            raise ValueError("No sequences pass quality threshold")
            
        filtered_sequences, filtered_scores = zip(*filtered_data)
        
        # Use weighted consensus on filtered sequences
        return self._weighted_consensus(filtered_sequences, filtered_scores)

    def evaluate_consensus(self, consensus, cluster_sequences, quality_scores=None):
        """
        Evaluate consensus sequence quality
        
        Returns:
        dict: Quality metrics including:
            - average_similarity: Mean similarity to consensus
            - coverage_uniformity: Uniformity of base coverage
            - confidence_score: Overall confidence in consensus
        """
        similarities = []
        for seq in cluster_sequences:
            similarity = self._sequence_similarity(consensus, seq)
            similarities.append(similarity)
        
        metrics = {
            'average_similarity': np.mean(similarities),
            'coverage_uniformity': self._calculate_coverage_uniformity(cluster_sequences),
            'confidence_score': np.mean(similarities) * self._calculate_coverage_uniformity(cluster_sequences)
        }
        
        return metrics

    def _sequence_similarity(self, seq1, seq2):
        """Calculate simple similarity score between sequences"""
        matches = sum(a == b for a, b in zip(seq1, seq2))
        return matches / max(len(seq1), len(seq2))

    def _calculate_coverage_uniformity(self, sequences):
        """Calculate how uniform the coverage is across positions"""
        max_length = max(len(seq) for seq in sequences)
        coverage = [0] * max_length
        
        for seq in sequences:
            for i in range(len(seq)):
                coverage[i] += 1
                
        # Calculate coefficient of variation
        coverage = np.array(coverage)
        return 1 - (np.std(coverage) / np.mean(coverage))