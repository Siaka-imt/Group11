# Group11
STAT4060J - Final project

Matrix factorization is fundamental to modern data science, yet classical deterministic algorithms
often struggle with computational efficiency on massive datasets or convergence stability in non-convex
problems. This project implements and compares Deterministic (Normal) versus Randomized approaches for 
three decomposition techniques: Singular Value Decomposition (SVD), Non-negative Matrix Factorization (NMF),
and CUR Decomposition. Using a custom-built R Shiny dashboard and real-world datasets (MovieLens and Volcano
Topography), we demonstrate that randomized initialization strategies significantly improve NMF convergence 
on sparse data, while Randomized SVD robustly captures natural low-rank structures in dense topographic data.
