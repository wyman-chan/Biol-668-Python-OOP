# Biol-668-Python-OOP

# Overview
This project demonstrates the use of Object-Oriented Programming (OOP) in Python to model biological sequences, including DNA, RNA, and proteins. Through a series of classes, the code replicates basic bioinformatics functionality such as sequence formatting, k-mer generation, codon translation, reverse complementation, and computation of protein properties—all without using specialized libraries like Biopython.

# Contents and Project Structure
CHAN_OOP_FinalProject_2023.py
Core Python script implementing the class hierarchy:

- Seq: Base class for generic sequences

- DNA: Inherits from Seq and adds DNA-specific utilities

- RNA: Inherits from DNA and adds transcription/translation logic

- Protein: Inherits from Seq and calculates properties like hydrophobicity and molecular weight

CHAN_PythonProject_2025.ipynb
A Jupyter Notebook showcasing usage examples and test cases for the implemented classes, with formatted output and explanatory cells. This uses CHAN_OOP_FinalProject_2023.py.

CHAN_04A_Biopython_Tutorial_LSH-1.ipynb
A Jupyter notebook exploring fundamental features of Biopython, focusing on biological sequence handling and Locality-Sensitive Hashing (LSH) techniques for efficient sequence similarity searches. This also uses the KRAS dependencies.

# Example
Examples were performed with the test files provided. 


CHAN_PythonProject_2025.ipynb
![7CC4EB4C-9F01-4EC5-9F80-441F1B01C145](https://github.com/user-attachments/assets/fcfcce33-5872-4576-ba70-dfd67014337e)

CHAN_04A_Biopython_Tutorial_LSH-1.ipynb
![CD9360CF-1F43-4BE2-AEC8-F74B0242A4FB](https://github.com/user-attachments/assets/86bcdffd-420a-40bb-8253-ac81ba2c079e)
