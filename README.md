# sc-type-py

How to Use:

This GitHub repository is designed to facilitate the use of sctype in Python. Follow these steps to get started:

**Step 1**: Cleaning Gene Sets
Execute "Step 1": Run the original gene_sets_prepare from sctype. This step is necessary to clean the gene set file. Since there isn't a straightforward and efficient way to validate gene symbols in Python, we'll stick to the R version for this task. The "Step 1.R" file will convert the input file into a cleaned and corrected HGNC symbol XLSX file.

**Step 2**: Calculating Scytype Score
Next, load your scaled data into Python and run the "Step 2.PY" file. This file will compute the scytype score, and the results will be saved in a .txt file for future reference. This file will serve as an annotation x cell file.

**Step 3**: Continuing Analysis
Take the annotation x cell file generated in Step 2 and incorporate it into your R workflow as usual. Proceed with the analysis using the "Step 3.R" file.


I'll prob be able to work on this later or something to make a checksymbols function that works a bit better/faster in python
