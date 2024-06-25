from math import log
from Bio import SeqIO, Entrez
from subprocess import run, call
import numpy as np
from hmmlearn import hmm
import pickle 

seq = 'AGCTCGCTCGCTGCGTATAAAATCGCATCGCGCGCAGC'


# Identificação de motifs específico em sequências
def find_motifs(seque: str):
    position = seq.find(seque)
    return position

# Exemplo livre de um motif
resideu_profile = {
'A':[ 61, 16,352, 3,354,268,360,222,155, 56, 83, 82, 82, 68, 77],
'C':[145, 46, 0, 10, 0, 0, 3, 2, 44,135,147,127,118,107,101],
'G':[152, 18, 2, 2, 5, 0, 10, 44,157,150,128,128,128,139,140],
'T':[ 31,309, 35,374, 30,121, 6,121, 33, 48, 31, 52, 61, 75, 71]}

def matchDnaProfile(seq, profile):
    """ Find the best-matching position and score when comparing a DNA
    sequence with a DNA sequence profile """
    bestScore = 0
    bestPosition = None # Just to start with
    width = len(profile['A'])
    for i in range(len(seq)-width):
        score = 0
        for j in range(width):
            letter = seq[i+j]
            score += profile[letter][j]
        if score > bestScore:
            bestScore = score
            bestPosition = i
    return bestScore, bestPosition

#Contando o conteúdo GC

def calcGcContent(seq, winSize=10):
    gcValues = []
    for i in range(len(seq)-winSize):
        subSeq = seq[i:i+winSize]
        numGc = subSeq.count('G') + subSeq.count('C')
        value = numGc/float(winSize)
        gcValues.append(value)
        
    return gcValues

#Contando a repetitividade das sequências
 
def calcRelativeEntropy(seq, resCodes):
    """Calculate a relative entropy value for the residues in a
    sequence compared to a uniform null hypothesis.
    """
    N = float(len(seq))
    base = 1.0/len(resCodes)
    prop = {}
    for r in resCodes:
        prop[r] = 0
    for r in seq:
        prop[r] += 1

    for r in resCodes:
        prop[r] /= N
    H = 0
    for r in resCodes:
        if prop[r] != 0.0:
            h = prop[r]* log(prop[r]/base, 2.0)
            H += h
    H /= log(base, 2.0)
    return H

name = ""
fileObj = open(name, "rU")


#Coletando dadosd do NCBI
from Bio import Entrez
Entrez.email = 'rodolfo.eli.jrles@gmail.com'
handle = Entrez.efetch(db="nucleotide", id="AY851612", rettype="gb", retmode="text")
print(handle.readline().strip())


Entrez.email = 'mickey@disney.com'
socketObj = Entrez.efetch(db="protein", rettype="fasta", id="71066805")
dnaObj = SeqIO.read(socketObj, "fasta")
socketObj.close()
print(dnaObj.description)
print(dnaObj.seq)

def SequenceIdentity(seqA, seqB):
    numPlaces = min(len(seqA), len(seqB))
    score = 0.0
    for nucleotide in range(numPlaces):
        if seqA[nucleotide] == seqB[nucleotide]:
            score+=1.0
    return 100.0 * score/numPlaces

def seqSimilarity(seqA,seqB,simMatrix):
    numPlaces = min(len(seqA), len(seqB))
    totalScore=0.0
    for i in range(numPlaces):
        
        residueA = seqA[i]
        residueB = seqB[i]
        
        totalScore =+ simMatrix[residueA][residueB]
        
        return totalScore

def pairAlignScore(alignA, alignB, simMatrix, insert=8, extend=4):
    totalScore = 0.0
    n = min(len(alignA), len(alignB))
    for i in range(n):
        residueA = alignA[i]
        residueB = alignB[i]
        
        if '-' not in (residueA, residueB):
            simScore = simMatrix[residueA][residueB]
        elif (i>0) and ('-' in (alignA[i-1],alignB[i-1])):
            simScore = -extend
        else:
            simScore = -insert
            
        totalScore += simScore
        
    return totalScore


def transition_matrix(states, n_states):
    trans_counts = np.zeros((n_states, n_states))
    
    for state_seq in states:
        for (s1, s2) in zip(state_seq[:-1], state_seq[1:]):
            trans_counts[s1, s2] += 1
    
    transmat_ = trans_counts / trans_counts.sum(axis=1, keepdims=True)
    return transmat_

def emission_matrix(states, observations, n_states, n_observations):
    emiss_counts = np.zeros((n_states, n_observations))
    
    for state_seq, obs_seq in zip(states, observations):
        for s, o in zip(state_seq, obs_seq):
            emiss_counts[s, o] += 1
    
    emissionprob_ = emiss_counts / emiss_counts.sum(axis=1, keepdims=True)
    return emissionprob_

def startprob_matrix(states, n_states): 
    start_counts = np.zeros(n_states)
    
    for state_seq in states:
        start_counts[state_seq[0]] += 1
    
    startprob_ = start_counts / start_counts.sum()
    return startprob_

def trainning_markov_model(HiddenStates, observation, n_HiddenStates:int, n_observations:int):
    
    emissions_array = np.concatenate([np.array(seq).reshape(-1, 1) for seq in observation])
    lengths = [len(seq) for seq in observation]
    
    model = hmm.MultinomialHMM(n_components=n_HiddenStates, n_iter=100, tol=0.01)
    
    model.fit(emissions_array, lengths=lengths)
    
    startprob_ = startprob_matrix(HiddenStates, n_HiddenStates)
    transmat_ = transition_matrix(HiddenStates, n_HiddenStates)
    emissionprob_ = emission_matrix(HiddenStates, observation, n_HiddenStates, n_observations)
    
    model.startprob_ = startprob_
    model.transmat_ = transmat_
    model.emissionprob_ = emissionprob_
    
    with open("algorithms/model.pkl", "wb") as arquivo:
	    pickle.dump(model, arquivo)
    
def run_markov_model(new_observations):
    
    model = pickle.load(open("model.pkl", "rb"))
    
    state_sequence = model.decode(new_observations, algorithm="viterbi")
    return state_sequence