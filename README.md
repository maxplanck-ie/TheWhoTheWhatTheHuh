This will be our new bcl to fastq pipeline. Features will include:

  * Do we ever have the same sample name split across libraries in the same project/flow cell? The current implementation will overwrite the fastqc results for that. A better option might be to make subdirectories within each FASTQC_projectID directory
  * It must run as a background process (ideally accepting HUP signals)
  * copy unmatched indexes (and other important information) to /sequencing_data/final
  * Add a quality index to the reporting and status email
  * Combine flowcells/lanes from submission (with our flowcell focus we still have issues if submissions are spread over many flow cells)
  * Compile an automated project-report (pdf), to go along with each submission
  * Functions divided into more meaningful subfiles in a module.
  * Written explicitly for python3, just to future proof things a bit.
