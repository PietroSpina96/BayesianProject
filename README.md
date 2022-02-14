# Distance-based probabilistic clustering for functional data
Somatosensory Evoked Potential (SEP) are electrical activity occurring in the central nervous system following a stimulus, detected within 2 seconds after the stimulus over 1600 instants.
SEP reflects the processing of information from the peripheral levels down to the cortical structures, and are perhaps useful in assessing the integrity of nerve conduction pathways and the different areas involved in the reception and processing of sensory stimuli.

The patients under study, in a state of coma, underwent a neuropsychological evaluation within 72 hours after surgery. We consider a dataset composed of 26 multivariate functional observation
with 4 components, each of which represents a Somatosensory Evoked Potential, detected within 2 seconds after the stimulus over 1600 time instants.

The aim of the project is to implement a Bayesian functional clustering algorithm in a Generalized Bayes Product Partition Model framework (GB-PPM), considering a generalization of the Mahalanobis distance in a functional setting.
At the beginning we will consider an uniform prior over the cluster partition with the number of clusters fixed a priori, then we will switch our focus over an EPPF Pitman-Yor process prior. Both algorithms will first be tested on simulated data and then on clinical ones.
