Starting with scBaseCamp


Creating a comprehensive catalog of all single-cell studies out there is Herculean effort - let alone getting access to all the data.
Therefore, we start with scBaseCamp, which currently has uniformly preprocessed counts matrices for all 10x Genomics single-cell data on SRA.


## scBaseCamp is not collectively exhaustive

It currently ignores datasets that are:
* Under lock and key on platforms like dbGap and EGA
* published only on an academic lab website
* non-10x Genomics modalities

That said, it is to be expected that Arc Institute will continue to close these gaps in the coming months and years.


## scBaseCamp does not provide detailed sample-level metadata

It only scrapes a few fields from the SRA study page.
This ignores much information from the associated publication.
Even finding the publication, let alone the supplementary files, that correspond to a SRP ID can be a challenging tasks.
On a sunny day, the GSE page contains a neat link to the publication, and it's open access.
On a rainy day, there is no link to a paper (even though a paper may exist), and we need to do e.g. a Google search to find it.
Even then, the paper might not be open-access.

