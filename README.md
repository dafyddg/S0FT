# S0FT
S0FT (Simple F0 Tracker): Simple online illustration of the S0FT F0 estimation module used in the CRAFT, CR0FT and SoundOps teaching packages.

S0FT is the name of the simple parametrised F0 estimator used in the CRAFT, S0FT and SoundOps package. It was created for teaching purposes, in particular to create awareness of the algorithmic factors affecting F0 visualisations.
The main features are:
1. Preprocessing:
  - centre-clipping
  - low and high pass filtering
2. F0 estimation:
  - FFT peak determination
  - Rising zero-crossings
  - Peak-picking (zero-crossings of first derivative)
3. Postprocessing:
  - Median filtering
 
For comparison, a visualisation using the RAPT 'gold standard' is included, as well as a spectrogram, in which the shape of the harmonics can be seen to parallel the shape of the F0 estimations. A figure of the output is included.
