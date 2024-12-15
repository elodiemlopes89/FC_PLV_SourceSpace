# FC_PLV_SourceSpace
Function created under the framework of the PhD thesis Novel Contributions to Personalized Brain Stimulation Biomarkers for Better Management of Neurological Disorders" - Doctoral Program in Biomedical Engineering (FEUP), supervised by João P. Cunha (INESC TEC, Porto, Portugal).

Function to estimate electroencephalogram (EEG)-based functional connectivity (FC) networks in the source space, using the Phase Locking Value (PLV) metric.

## PLV_ud_source.m
This function computes the phase locking value (PLV) between EEG sources in different brain regions for a given patient. It performs preprocessing, segmentation, source analysis, and computes undirected (PLVu) and directed (PLVd) network measures. Additionally, it outputs global connectivity metrics for the given EEG segments. Estimation of Brain Functional Connectivity (FC) Networks in the source space Fieldtrip is used to construct the Head Model

Headmodel: brainstorm_ICBM152_mni

### Syntax:
out_seg = PLV_ud_source(PatID, N, Stim, FR)

INPUTS
* PatID: A string representing the patient ID (e.g., '01').
* N: A string or number representing the segment number (e.g., '1').
* Stim: A character ('a' or 'b') indicating whether the data is from after ('a') or before ('b') stimulation.
* FR: A vector of two values, [fmin, fmax], specifying the frequency range of interest for filtering (e.g., [1, 30] Hz).

OUTPUTS

out_seg: A structure containing the following fields:
* eeg: EEG data for the segment.
* sources: The source-reconstructed EEG signals.
* net_u: Undirected network connectivity matrix.
* net_d: Directed network connectivity matrix.
* GC: Global connectivity metrics (Clustering coefficient, Node degree, etc.)
* ROIs: Labels for the regions of interest.

### Packages:
* Fieldtrip - 20191119 (@Free Software Foundation, Inc, 1991)
* BrainNet Viewer (https://www.mathworks.com/matlabcentral/fileexchange/68881-brainnet-viewer)
* Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/)
* EEG Preprocessing custom-made package (https://github.com/elodiemlopes89/Electrophysiological-Data-Preprocessing-Visualization)


### Pipeline

1. PREPROCESSING DATA
* load data
* remove irrelevant channels
* set frequency range for filtering

2. ELECTRODE INFO
* Load electrode information from the standard 1005 electrode layout (standard_1005.elc)

3. DEFINE ELECTRODE POSITIONS FOR THE PATIENT
* In this example, we also defined additional positions for non-standard elecrodes: F11, F12, P11, P12

4. SEGEMENTATION OF EEG DATA
* Segment epochs into 20s to estimate FC networks, as recomended in the literature.

5. CONVERT DATA INTO FIELDTRIP STRUCUTRE
* Using the data2fieldtrip.m function

6. PREPARE HEAD MODEL
* Using the template brainstorm_ICBM152_mni.

7. PREPARE ELECTRODE MODEL

8. PREPARE SOURCE GRID

9. PREPARE LEADFIELD.

