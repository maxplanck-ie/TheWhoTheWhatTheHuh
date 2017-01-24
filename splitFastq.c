#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>

int isDuplicate(char *line) {
    char *p;
    size_t l = strlen(line) - 11;

    p = line + l;
    if(strcmp(p, " duplicate\n") == 0) {
        //Strip " duplicate" from the end, fixing the line ending
        *(p++) = '\n';
        *p = 0;
        return 1;
    }
    return 0;
}

int processSingle(FILE *ifile, char *bname, uint64_t *total, uint64_t *dupes, char* pigz_threads) {
    char *cmd = NULL;
    FILE *o1=NULL, *o1d=NULL;
    FILE *out;
    char *line = malloc(1024);
    if(!line) return 1;

    //Non-duplicate output
    cmd = calloc(strlen(bname) + strlen("_R1.fastq.gz pigz -p >  ") + strlen(pigz_threads), sizeof(char));
    if(!cmd) return 1; 
    sprintf(cmd, "pigz -p %s > %s_R1.fastq.gz", pigz_threads, bname);
    o1 = popen(cmd, "w");
    free(cmd);

    //duplicate output
    cmd = calloc(strlen(bname) + strlen("_R1.optical_duplicates.fastq.gz pigz -p >  ") + strlen(pigz_threads), sizeof(char));
    if(!cmd) return 1;
    sprintf(cmd, "pigz -p %s > %s_R1_optical_duplicates.fastq.gz", pigz_threads, bname);
    o1d = popen(cmd, "w");
    free(cmd);

    *total = 0;
    *dupes = 0;
    //Iterate through the lines, writing as appropriate
    while((line = fgets(line, 1024, ifile))) {
        // Check if this is a duplicate
        if(isDuplicate(line)) {
             *dupes += 1;
             out = o1d;
        } else {
             out = o1;
        }
        //Read 1
        fputs(line, out);
        if(!(line = fgets(line, 1024, ifile))) return 1;
        fputs(line, out);
        if(!(line = fgets(line, 1024, ifile))) return 1;
        fputs(line, out);
        if(!(line = fgets(line, 1024, ifile))) return 1;
        fputs(line, out);
        *total += 1;
    }
    free(line);

    return 0;
}
    
int processPaired(FILE *ifile, char *bname, uint64_t *total, uint64_t *dupes, char *pigz_threads) {
    char *cmd = NULL;
    FILE *o1=NULL, *o2=NULL, *o1d=NULL, *o2d=NULL;
    FILE *out1, *out2;
    char *line = malloc(1024);
    if(!line) return 1;

    //Non-duplicate output
    cmd = calloc(strlen(bname) + strlen("_R1.fastq.gz pigz -p >  ") + strlen(pigz_threads), sizeof(char));
    if(!cmd) return 1; 
    sprintf(cmd, "pigz -p %s > %s_R1.fastq.gz", pigz_threads, bname);
    o1 = popen(cmd, "w");
    sprintf(cmd, "pigz -p %s > %s_R2.fastq.gz", pigz_threads, bname);
    o2 = popen(cmd, "w");
    free(cmd);

    //duplicate output
    cmd = calloc(strlen(bname) + strlen("_R1.optical_duplicates.fastq.gz pigz -p >  ") + strlen(pigz_threads), sizeof(char));
    if(!cmd) return 1; 
    sprintf(cmd, "pigz -p %s > %s_R1_optical_duplicates.fastq.gz", pigz_threads, bname);
    o1d = popen(cmd, "w");
    sprintf(cmd, "pigz -p %s > %s_R2_optical_duplicates.fastq.gz", pigz_threads, bname);
    o2d = popen(cmd, "w");
    free(cmd);

    *total = 0;
    *dupes = 0;
    //Iterate through the lines, writing as appropriate
    while((line = fgets(line, 1024, ifile)) != NULL) {
        // Check if this is a duplicate
        if(isDuplicate(line)) {
             *dupes += 1;
             out1 = o1d;
             out2 = o2d;
        } else {
             out1 = o1;
             out2 = o2;
        }
        //Read 1
        fputs(line, out1);
        if(!(line = fgets(line, 1024, ifile))) return 1;
        fputs(line, out1);
        if(!(line = fgets(line, 1024, ifile))) return 1;
        fputs(line, out1);
        if(!(line = fgets(line, 1024, ifile))) return 1;
        fputs(line, out1);
        //Read 2
        if(!(line = fgets(line, 1024, ifile))) return 1;
        fputs(line, out2);
        if(!(line = fgets(line, 1024, ifile))) return 1;
        fputs(line, out2);
        if(!(line = fgets(line, 1024, ifile))) return 1;
        fputs(line, out2);
        if(!(line = fgets(line, 1024, ifile))) return 1;
        fputs(line, out2);
        *total += 1;
    }
    free(line);

    return 0;
}
    

int main(int argc, char *argv[]) {
    int PE, rv;
    char *bname, *cmd=NULL;
    uint64_t total=0, dupes=0;
    char *pigz_threads = "1";
    FILE *ifile;

    if(argc < 4) {
        fprintf(stderr, "Usage %s input.fastq.gz paired basename [pigzThreads]\n", argv[0]);
        fprintf(stderr, "Splits an interleaved single or paired-end file into 2/4 output files according\n"
"to whether 'duplicate' is in the read name. This is useful to split the output\n"
"of clumpify into duplicates and non-duplicates. The metrics are then written to\n"
"basename.deduplicate.txt and files to basename_R1.fastq.gz,\n"
"basename_R1_optical_duplicates.fastq.gz and so on.\n\n");
        fprintf(stderr, "paired:      1 is the input is PE, 0 otherwise\n");
        fprintf(stderr, "basename:    The basename for the output files\n");
        fprintf(stderr, "pigzThreads: (optional) Number of threads used by each pigz process. Note that there can be 4 of these.\n");
        return 0;
    }

    PE = atoi(argv[2]);
    bname = argv[3];
    if(argc == 5) pigz_threads = argv[4];

    cmd = calloc(strlen(argv[1]) + strlen("zcat  "), sizeof(char));
    assert(cmd);
    sprintf(cmd, "zcat %s", argv[1]);
    ifile = popen(cmd, "r");

    if(PE==1) rv = processPaired(ifile, bname, &total, &dupes, pigz_threads);
    else rv = processSingle(ifile, bname, &total, &dupes, pigz_threads);
    pclose(ifile);

    if(rv) {
        fprintf(stderr, "We encountered some sort of error!\n");
        return 1;
    }

    //Output the metrics to a text file
    cmd = calloc(strlen(bname) + strlen(".duplicate.txt "), sizeof(char));
    sprintf(cmd, "%s.duplicate.txt", bname);
    ifile = fopen(cmd, "w");
    free(cmd);
    fprintf(ifile, "%"PRIu64"\t%"PRIu64"\n", dupes, total);
    fclose(ifile);

    return 0;
}
